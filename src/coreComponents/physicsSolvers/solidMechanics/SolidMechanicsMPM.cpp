/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsMPM.cpp
 */

#include "SolidMechanicsMPM.hpp"
#include "codingUtilities/Utilities.hpp"

#include "kernels/ImplicitSmallStrainNewmark.hpp"
#include "kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "kernels/ExplicitSmallStrain.hpp"
#include "kernels/ExplicitFiniteStrain.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "common/TypeDispatch.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/Kinematics.h"
#include "LvArray/src/output.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/TractionBoundaryCondition.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "events/mpmEvents/MPMEvents.hpp"
#include "physicsSolvers/solidMechanics/LogLevelsInfo.hpp"

#include "chrono"
#include "thread"
#include <random>

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

// Flattened combinations function to avoid performance hit from recursive calls
GEOS_HOST_DEVICE
array2d< int > generateCombinations(array1d< array1d< int > > sets)
{
  int numSets = sets.size();
  int numCombinations = 1;
  array1d< int > m(numSets);
  for( int s=0; s < numSets; s++)
  {
    numCombinations *= sets[s].size();
    m[s] = numCombinations;
  }

  array2d< int > combinations( numCombinations, numSets);
  for( int c = 0; c < numCombinations; c++ )
  {
    for( int s = 0; s < numSets; s++ )
    {
      int i = 0;
      if ( s == 0 )
      {
        i = c % m[s];
      }
      else
      {
        i = ( c % m[s] ) / m[s-1];
      }
      combinations[c][s] = sets[s][i];
    }
  }

  return combinations;
}

// // Swap the values at indices i and j with one another
// template< typename T >
// GEOS_HOST_DEVICE
// void swap( T * vec, localIndex i, localIndex j )
// {
//   T temp = vec[i];
//   vec[i] = vec[j];
//   vec[j] = temp;
// }

// template< typename T >
// GEOS_HOST_DEVICE
// localIndex partition( T * vec, localIndex * indices,  localIndex low, localIndex high )
// {
//     // Selecting last element as the pivot
//     T pivot = vec[high];

//     // Index of elemment just before the last element
//     // It is used for swapping
//     localIndex i = (low - 1);

//     for (localIndex j = low; j <= high - 1; ++j) {

//         // If current element is smaller than or
//         // equal to pivot
//         if (vec[j] <= pivot) {
//             ++i;
//             swap(vec, i, j);
//             swap(indices, i, j);
//         }
//     }

//     // Put pivot to its position
//     swap(vec, i+1, high);
//     swap(indices, i+1, high);

//     // Return the point of partition
//     return (i + 1);
// }

// template< typename T >
// GEOS_HOST_DEVICE
// void quickSortRecursive( T * vec, localIndex * indices, localIndex low, localIndex high )
// {
//   // Base case: This part will be executed till the starting
//   // index low is lesser than the ending index high
//   if (low < high)
//   {
//       // arr[p] is now at right place
//       localIndex partitionIndex = partition( vec, indices, low, high );

//       // Separately sort elements before and after the partition Index pi
//       quickSortRecursive( vec, indices, low,  partitionIndex - 1 );
//       quickSortRecursive( vec, indices, partitionIndex + 1, high );
//   }
// }

// template< typename T >
// GEOS_HOST_DEVICE
// void quickSort( T * vec, localIndex * indices )
// {
//   // Should do a size check, or at least pass in a size since we won't know from the pointers
//   for( localIndex i = 0; i < 64; ++i)
//   {
//     indices[i] = i;
//   }
//   quickSortRecursive( vec, indices, 0, 63 );
// }

GEOS_HOST_DEVICE
void swap(localIndex* a, localIndex* b)
{
  localIndex temp =*a;
  *a = *b;
  *b = temp;
}

// GEOS_HOST_DEVICE
// localIndex partition(localIndex (& arr)[64], localIndex (& arr2)[64], localIndex low, localIndex high)
// {
//   localIndex p  = arr[low];
//   localIndex i = low;
//   localIndex j = high;

//   while(i<j)
//   {
//     while(arr[i] <= p && i <= high - 1 )
//     {
//       ++i;
//     }

//     while(arr[j] > p && j >= low + 1 )
//     {
//       j--;
//     }

//     if(i < j)
//     {
//       swap(&arr[i], &arr[j]);
//       swap(&arr2[i], &arr2[j]);
//     }
//   }
//   swap(&arr[low], &arr[j]);
//   swap(&arr2[low], &arr2[j]);
//   return j;
// }

// GEOS_HOST_DEVICE
// void quickSort(localIndex (& arr)[64], localIndex (& arr2)[64], localIndex low, localIndex high)
// {
//   if(low < high)
//   {
//     localIndex pi = partition(arr, arr2, low, high);

//     quickSort(arr, arr2, low, pi-1);
//     quickSort(arr, arr2, pi+1, high);
//   }
// }

GEOS_HOST_DEVICE
localIndex partitionIterative(localIndex (& arr)[64], localIndex (& arr2)[64], int l, int h)
{
    localIndex x = arr[h];
    localIndex i = (l - 1);

    for (localIndex j = l; j <= h - 1; ++j) {
        if (arr[j] <= x) {
            ++i;
            swap(&arr[i], &arr[j]);
            swap(&arr2[i], &arr2[j]);
        }
    }
    swap(&arr[i + 1], &arr[h]);
    swap(&arr2[i + 1], &arr2[h]);
    return (i + 1);
}

GEOS_HOST_DEVICE
void quickSortIterative(localIndex (& arr)[64], localIndex (& arr2)[64], localIndex l, localIndex h)
{
    // Create an auxiliary stack
    localIndex stack[64]; //[h - l + 1];

    // initialize top of stack
    localIndex top = -1;

    // push initial values of l and h to stack
    stack[++top] = l;
    stack[++top] = h;

    // Keep popping from stack while is not empty
    while (top >= 0) {
        // Pop h and l
        h = stack[top--];
        l = stack[top--];

        // Set pivot element at its correct position
        // in sorted array
        localIndex p = partitionIterative(arr, arr2, l, h);

        // If there are elements on left side of pivot,
        // then push left side to stack
        if (p - 1 > l) {
            stack[++top] = l;
            stack[++top] = p - 1;
        }

        // If there are elements on right side of pivot,
        // then push right side to stack
        if (p + 1 < h) {
            stack[++top] = p + 1;
            stack[++top] = h;
        }
    }
}

SolidMechanicsMPM::SolidMechanicsMPM( const string & name,
                                      Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_solverProfiling( 0 ),
  m_logStartCycle( 0 ),
  m_plottableFields(),
  m_plottableFieldsSorted(),
  m_plotGridFields( 0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_updateMethod( UpdateMethodOption::FLIP ),
  m_updateOrder( 2 ),
  m_prescribedBcTable( 0 ),
  m_boundaryConditionTypes(),
  m_boundaryFaceCoefficientsOfRestitution(),
  m_boundaryFaceFrictionCoefficients(),
  m_bcTable(),
  m_prescribedFTable( 0 ),
  m_prescribedBoundaryFTable( 0 ),
  m_fTableInterpType( SolidMechanicsMPM::InterpolationOption::Linear ),
  m_fTable(),
  m_domainF(),
  m_domainL(),
  m_enablePrescribedBoundaryTransverseVelocities(),
  m_prescribedBoundaryTransverseVelocities(),
  m_globalFaceReactions(),
  m_bodyForce(),
  m_enableSurfaceTension( 0 ),
  m_surfaceTensionCoefficient( 0.0 ),
  m_enableBoreholePressure( 0 ),
  m_boreholeStress(),
  m_boreholeRadius( 0.0 ),
  m_stressControl(),
  m_stressTableInterpType( SolidMechanicsMPM::InterpolationOption::Linear ),
  m_stressTable(),
  m_stressControlKp( 0.1 ),
  m_stressControlKi( 0.0 ),
  m_stressControlKd( 0.0 ),
  m_domainStress(),
  m_stressControlLastError(),
  m_stressControlITerm(),
  m_boxAverageHistory( 0 ),
  m_boxAverageWriteInterval( 0.0 ),
  m_nextBoxAverageWriteTime( 0.0 ),
  m_boxAverageMin(),
  m_boxAverageMax(),
  m_reactionHistory( 0 ),
  m_reactionWriteInterval( 0.0 ),
  m_nextReactionWriteTime( 0.0 ),
  m_writeParticleData( 0 ),
  m_particleDataWriteInterval( 0.0 ),
  m_nextParticleDataWriteTime( 0.0 ),
  m_explicitSurfaceNormalInfluence( 0.0 ),
  m_computeParticleSurfaceNormalsAndPositions( 0 ),
  m_normalAndPositionMethod( NormalsAndPositionsMethodOption::LR ),
  m_overwriteExistingNormalsAndPositions( 0 ),
  m_gpuScheme( GPUSchemeOption::Atomics ),
  m_referenceCohesiveZone( 0 ),
  m_enableCohesiveLaws( 0 ),
  m_cohesiveLaw( CohesiveLawOption::NeedlemanXu ),
  m_enableCohesiveFailure( 0 ),
  m_preventCZInterpentration( 0 ),
  m_normalForceConstant( 0.0 ),
  m_shearForceConstant( 0.0 ),
  m_reinitializeCohesiveZones( 0 ),
  m_totalBinderVolume( 0.0 ),
  m_numSurfaceIntegrationPoints( 200 ),
  m_maxCohesiveNormalStress( 0 ),
  m_maxCohesiveShearStress( 0 ),
  m_characteristicNormalDisplacement( 1 ),
  m_characteristicTangentialDisplacement( 1 ),
  m_maxCohesiveNormalDisplacement( 0 ),
  m_maxCohesiveTangentialDisplacement( 0 ),
  m_polymerCZThickness( 1.0 ),
  m_polymerCZBulkModulus( 1.0 ),
  m_polymerCZShearModulus( 1.0 ),
  m_polymerCZYieldStrength0( 0.0 ),
  m_polymerCZR0( 0.0 ),
  m_polymerCZR1( 0.0 ),
  m_polymerCZR2( 1.0 ),
  m_polymerCZGr( 0.0 ),
  m_polymerCZMaxStretch( 0.0 ),
  m_needsNeighborList( 0 ),
  m_needsNodalNeighborList( 0 ),
  m_logisticRegression( 0 ),
  m_neighborRadius( -1.0 ),
  m_binSizeMultiplier( 1 ),
  m_maxNumNeighbors( 0 ),
  m_thinFeatureDFGThreshold( DBL_MAX ),
  m_FSubcycles( 1 ),
  m_LBar( 0 ),
  m_LBarScale( 0.0 ),
  m_exactJIntegration( 0 ),
  m_maxParticleVelocity( 1e6 ), // Floating point exception if this is set to DBL_MAX when squared
  m_maxParticleVelocitySquared( DBL_MAX ),
  m_minParticleJacobian( 0.1 ),
  m_maxParticleJacobian( 10.0 ),
  m_overlapCorrection( OverlapCorrectionOption::Off ),
  m_overlapThreshold1( 1.00 ),
  m_overlapThreshold2( 1.10 ),
  m_computeSPHJacobian( 0 ),
  m_shockHeating( 0 ),
  m_computeInternalEnergyAndTemperature( 0 ),
  m_useArtificialViscosity( 0 ),
  m_artificialViscosityQ0( 0.0 ),
  m_artificialViscosityQ1( 0.0 ),
  m_cpdiDomainScaling( 0 ),
  m_subdivideParticles( 0 ),
  m_disableSurfaceNormalsAndPositionsOnCPDIScaling( 0 ),
  m_disableSurfaceNormalsAndPositionsOnDamage( 0 ),
  m_smallMass( DBL_MAX ),
  m_numContactGroups( 0 ),
  m_numContactFlags( 0 ),
  m_numVelocityFields( 0 ),
  m_hasContact( 0 ),
  m_computeNodalArea( 0 ),
  m_useNodePosForArea( 0 ),
  m_separabilityMinDamage( 0.5 ),
  m_treatFullyDamagedAsSingleField( 0 ),
  m_surfaceDetection( 0 ),
  m_damageFieldPartitioning( 0 ),
  m_useSurfacePositionForContact( 0 ),
  m_contactNormalType( ContactNormalTypeOption::MassWeighted ),
  m_contactNormalExponent( 1.0 ),
  m_contactGapCorrection( ContactGapCorrectionOption::Simple ),
  m_directionalOverlapCorrection( 0 ),
  m_resetDefGradForFullyDamagedParticles( 0 ),
  m_useReferenceVectorsForParticleUpdate( 0 ),
  m_plotUnscaledParticles( 0 ),
  m_frictionCoefficient( -1.0 ),
  m_frictionCoefficientTable(),
  m_planeStrain( 0 ),
  m_numDims( 3 ),
  m_generalizedVortexMMS( 0 ),
  m_ijkMap(),
  m_numberOfSubRegions( 0 ),
  m_useEvents( 0 ),
  m_mpmEventManager( nullptr ),
  m_surfaceHealing( 0 ),
  m_computeXProfile( 0 ),
  m_xProfileWriteInterval( 0.0 ),
  m_nextXProfileWriteTime( 0.0 ),
  m_xProfileVx0( 0.0 ),
  m_implicitContinuumFluidPressure( 0.0 )
{
  GEOS_LOG_RANK("MPM Constructor");
  // setInputFlags( InputFlags::OPTIONAL );

  registerWrapper( "solverProfiling", &m_solverProfiling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for timing subroutines in the solver" );

  registerWrapper( "logStartCycle", &m_logStartCycle ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Cycle to start logging output for debugging" );

  registerWrapper( "plottableFields", &m_plottableFields ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "If specified, determines which fields are written to vtk. Helpful for reducing the vtk output file sizes" );

  registerWrapper( "plottableFieldsSorted", &m_plottableFieldsSorted ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Sorted array of plottable fields" );

  registerWrapper( "plotGridFields", &m_plotGridFields ).
    setApplyDefaultValue( m_plotGridFields ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable plotting of grid fields in vtk output" );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( "updateMethod", &m_updateMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_updateMethod ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Update method. Options are:\n* " + EnumStrings< UpdateMethodOption >::concat( "\n* ") );

  registerWrapper( "updateOrder", &m_updateOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_updateOrder ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription("Order for update method, only applies to XPIC and FMPM");

  registerWrapper( "prescribedBcTable", &m_prescribedBcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_prescribedBcTable ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to have time-dependent boundary condition types" );

  registerWrapper( "boxAverageHistory", &m_boxAverageHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_boxAverageHistory ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to output box average history" );

  registerWrapper( "boxAverageWriteInterval", &m_boxAverageWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_boxAverageWriteInterval ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval between writing box averages to files" );

  registerWrapper( "boxAverageMin", &m_boxAverageMin).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Minimum corner position of box average" );

  registerWrapper( "boxAverageMax", &m_boxAverageMax).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Maximum corner position of box average" );

  registerWrapper( "nextBoxAverageWriteTime", &m_nextBoxAverageWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next time to write box averages" );

  registerWrapper( "reactionHistory", &m_reactionHistory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to output face reaction history" );

  registerWrapper( "reactionWriteInterval", &m_reactionWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_reactionWriteInterval ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval between writing reactions to files" );

  registerWrapper( "nextReactionWriteTime", &m_nextReactionWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next time to write reactions" );

  registerWrapper( "writeParticleData", &m_writeParticleData ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_writeParticleData ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable writing particle data to file" );

  registerWrapper( "particleDataWriteInterval", &m_particleDataWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_particleDataWriteInterval ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval to write particle data to file" );

  registerWrapper( "nextParticleDataWriteTime", &m_nextParticleDataWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_nextParticleDataWriteTime ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Next time to write particle data" );

  registerWrapper( "explicitSurfaceNormalInfluence", &m_explicitSurfaceNormalInfluence ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_explicitSurfaceNormalInfluence ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Determines the relative weighting for implicit and explicit surface normals" );

  registerWrapper( "computeParticleSurfaceNormalsAndPositions", &m_computeParticleSurfaceNormalsAndPositions ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_computeParticleSurfaceNormalsAndPositions ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to compute particle surface normals and positions" );

  registerWrapper( "normalAndPositionMethod", &m_normalAndPositionMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_normalAndPositionMethod ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Method for computing particle surface normals and positions" );

  registerWrapper( "overwriteExistingNormalsAndPositions", &m_overwriteExistingNormalsAndPositions ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_overwriteExistingNormalsAndPositions ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to that determines if existing particle normals and positions are overriden" );

  registerWrapper( "boundaryConditionTypes", &m_boundaryConditionTypes ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Boundary conditions on x-, x+, y-, y+, z- and z+ faces. Options are:\n* " + EnumStrings< BoundaryConditionOption >::concat( "\n* " ) );

 registerWrapper( "boundaryFaceCoefficientsOfRestitution", &m_boundaryFaceCoefficientsOfRestitution ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Boundary face coefficients of restitution for BC type = 3 on x-, x+, y-, y+, z- and z+ faces." );

 registerWrapper( "boundaryFaceFrictionCoefficients", &m_boundaryFaceFrictionCoefficients ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Boundary face fricition coefficients for BC type = 3 on x-, x+, y-, y+, z- and z+ faces." );

  registerWrapper( "bodyForce", &m_bodyForce ).
    setInputFlag(InputFlags::OPTIONAL).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores uniform body force" );

  registerWrapper( "enableSurfaceTension", &m_enableSurfaceTension ).
    setApplyDefaultValue( m_enableSurfaceTension ).
    setInputFlag(InputFlags::OPTIONAL).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable surface tension forces on particles" );

  registerWrapper( "surfaceTensionCoefficient", &m_surfaceTensionCoefficient ).
    setApplyDefaultValue( m_surfaceTensionCoefficient ).
    setInputFlag(InputFlags::OPTIONAL).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Stores the value of the surface tension coefficient" );

  registerWrapper( "enableBoreholePressure", &m_enableBoreholePressure ).
    setInputFlag(InputFlags::FALSE).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable borehole pressure" );

  registerWrapper( "boreholeStress", &m_boreholeStress ).
    setInputFlag(InputFlags::FALSE).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores borehole stress state" );

  registerWrapper( "boreholeRadius", &m_boreholeRadius ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Radius of borehole" );

  registerWrapper( "bcTable", &m_bcTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores time-dependent bc types on x-, x+, y-, y+, z- and z+ faces." );

  registerWrapper( "prescribedFTable", &m_prescribedFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to have time-dependent superimposed velocity gradient for triply periodic simulations" );

  registerWrapper( "prescribedBoundaryFTable", &m_prescribedBoundaryFTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to have time-dependent boundary conditions described by a global background grid F" );

  registerWrapper( "fTableInterpType", &m_fTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setApplyDefaultValue( m_fTableInterpType ).
    setDescription( "The type of F table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "fTable", &m_fTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores time-dependent grid-aligned stretches interpreted as a gloabl background grid F read from the XML file." );

  registerWrapper( "enablePrescribedBoundaryTransverseVelocities", &m_enablePrescribedBoundaryTransverseVelocities ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array storing a flag to enable transverse velocity boundary conditions for each face of domain" );

  registerWrapper( "prescribedBoundaryTransverseVelocities", &m_prescribedBoundaryTransverseVelocities ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array storing the transverse velocity boundary conditions for each face of domain" );

  registerWrapper( "stressControl" , &m_stressControl).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether stress control using box averages is enabled" );

  registerWrapper( "stressTableInterpType", &m_stressTableInterpType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setApplyDefaultValue( m_stressTableInterpType ).
    setDescription( "The type of stress table interpolation. Options are 0 (linear), 1 (cosine), 2 (quintic polynomial)." );

  registerWrapper( "stressTable", &m_stressTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Array that stores the time-depended grid aligned stresses" );

  registerWrapper( "stressControlKp", &m_stressControlKp ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Proportional gain of stress PID controller" );

  registerWrapper( "stressControlKi", &m_stressControlKi ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Integral gain of stress PID controller" );

   registerWrapper( "stressControlKd", &m_stressControlKd ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Derivative gain of stress PID controller" );

  registerWrapper( "domainStress", &m_domainStress ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stores current target domain stress as driven by stress control" );

  registerWrapper( "stressControlLastError", &m_stressControlLastError ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stress control error from previous timestep" );

  registerWrapper( "stressControlITerm", &m_stressControlITerm ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stress control integral term" );

  registerWrapper( "referenceCohesiveZone", &m_referenceCohesiveZone ).
    setInputFlag( InputFlags::FALSE).
    setApplyDefaultValue( m_referenceCohesiveZone ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to reference cohesive zone" );

  registerWrapper( "enableCohesiveLaws", &m_enableCohesiveLaws ).
    setInputFlag( InputFlags::FALSE).
    setApplyDefaultValue( m_enableCohesiveLaws ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag to enable cohesive laws" );

  registerWrapper( "cohesiveLaw", &m_cohesiveLaw ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_cohesiveLaw ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Type of cohesive law" );

  registerWrapper( "enableCohesiveFailure", &m_enableCohesiveFailure ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_enableCohesiveFailure ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable cohesive failure" );

  registerWrapper( "preventCZInterpentration", &m_preventCZInterpentration ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_preventCZInterpentration ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to use contact normal forces instead of CZ forces in compression" );

  registerWrapper( "normalForceConstant", &m_normalForceConstant ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_normalForceConstant ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Normal force constant for uncoupled cohesive law" );

  registerWrapper( "shearForceConstant", &m_shearForceConstant ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_shearForceConstant ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Shear force constant for uncoupled cohesive law" );

  registerWrapper( "reinitializeCohesiveZones", &m_reinitializeCohesiveZones ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_reinitializeCohesiveZones ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Triggers reinitialization of cohesive zones" );

  registerWrapper( "totalBinderVolume", &m_totalBinderVolume ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_totalBinderVolume ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Total volume of binder" );

  registerWrapper( "gpuScheme", &m_gpuScheme ).
    setInputFlag( InputFlags::OPTIONAL).
    setApplyDefaultValue( m_gpuScheme ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "GPU Scheme (Atomics, No atomics, Random shuffling, Minimal atomics)" );

  registerWrapper( "numSurfaceIntegrationPoints", &m_numSurfaceIntegrationPoints ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_numSurfaceIntegrationPoints ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Number of surface integration points for cohesive grid nodes" );

  registerWrapper( "maxCohesiveNormalStress", &m_maxCohesiveNormalStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Max cohesive normal stress value" );

  registerWrapper( "maxCohesiveShearStress", &m_maxCohesiveShearStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Max cohesive shear stress value" );

  registerWrapper( "characteristicNormalDisplacement", &m_characteristicNormalDisplacement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Characteristic normal displacement value" );

  registerWrapper( "characteristicTangentialDisplacement", &m_characteristicTangentialDisplacement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Characteristic tangential displacement value" );

  registerWrapper( "maxCohesiveNormalDisplacement", &m_maxCohesiveNormalDisplacement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Max cohesive normal displacement value" );

  registerWrapper( "maxCohesiveTangentialDisplacement", &m_maxCohesiveTangentialDisplacement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Max cohesive tangential displacement value" );

  registerWrapper( "polymerCZThickness", &m_polymerCZThickness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law thickness" );

  registerWrapper( "polymerCZBulkModulus", &m_polymerCZBulkModulus ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law bulk modulus" );

  registerWrapper( "polymerCZShearModulus", &m_polymerCZShearModulus ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law shear modulus" );

  registerWrapper( "polymerCZYieldStrength0", &m_polymerCZYieldStrength0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law initial yield strength" );

  registerWrapper( "polymerCZR0", &m_polymerCZR0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law r0 parameter" );

  registerWrapper( "polymerCZR1", &m_polymerCZR1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law r1 parameter" );

  registerWrapper( "polymerCZR2", &m_polymerCZR2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law r2 parameter" );

  registerWrapper( "polymerCZGr", &m_polymerCZGr ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law Gr parameter" );

  registerWrapper( "polymerCZMaxStretch", &m_polymerCZMaxStretch ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Polymer cohesive law max stretch" );

  registerWrapper( "cohesiveNodeGlobalIndices", &m_cohesiveNodeGlobalIndices ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Array of the global indices for cohesive grid nodes" );

  registerWrapper( "referenceCohesiveGridNodePartitioningSurfaceNormals", &m_referenceCohesiveGridNodePartitioningSurfaceNormals ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference partitioning surface normal for cohesive grid nodes" );

  registerWrapper( "referenceCohesiveGridNodeAreas", &m_referenceCohesiveGridNodeAreas ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node areas" );

  registerWrapper( "referenceCohesiveGridNodePositions", &m_referenceCohesiveGridNodePositions ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node positions" );

  registerWrapper( "cohesiveGridNodeDamages", &m_cohesiveGridNodeDamages ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Cohesive grid node damages" );

  registerWrapper( "referenceCohesiveGridNodeSurfaceNormals", &m_referenceCohesiveGridNodeSurfaceNormals ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Reference cohesive grid node surface normals" );

  registerWrapper( "maxCohesiveGridNodeNormalDisplacement", &m_maxCohesiveGridNodeNormalDisplacement ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Max cohesive normal displacement for each cohesive grid node" );

  registerWrapper( "maxCohesiveGridNodeTangentialDisplacement", &m_maxCohesiveGridNodeTangentialDisplacement ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Max cohesive tangential displacement for each cohesive grid node" );

  registerWrapper( "needsNeighborList", &m_needsNeighborList ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_needsNeighborList ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to construct neighbor list" );

  registerWrapper( "needsNodalNeighborList", &m_needsNodalNeighborList ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_needsNodalNeighborList ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to construct nodal neighbor list" );

  registerWrapper( "logisticRegression", &m_logisticRegression ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_logisticRegression ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for whether to utilize logistic regression for contact" );

  registerWrapper( "neighborRadius", &m_neighborRadius ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_neighborRadius ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Neighbor radius for SPH-type calculations" );

  registerWrapper( "binSizeMultiplier", &m_binSizeMultiplier ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_binSizeMultiplier ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Multiplier for setting bin size, used to speed up particle neighbor sorting" );

  registerWrapper( "thinFeatureDFGThreshold", &m_thinFeatureDFGThreshold ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_thinFeatureDFGThreshold ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Threshold to treat relatively thin ( compared to grid spacing ) damaged material to avoid spurious slip surfaces" );

  registerWrapper( "useDamageAsSurfaceFlag", &m_useDamageAsSurfaceFlag ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Indicates whether particle damage at the beginning of the simulation should be interpreted as a surface flag" );

  registerWrapper( "FSubcycles", &m_FSubcycles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_FSubcycles ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Number of sub cycles to more accurately integrate the deformation gradient" );

  registerWrapper( "LBar", &m_LBar ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_LBar ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Used to choose the Lbar method. 0 = do nothing, 1 = run trD through the SPH kernel, 2 = straight-up volume average" );

  registerWrapper( "LBarScale", &m_LBarScale ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_LBarScale ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Scales the contribution of volume-averaged trD to particle L, alleviates checkerboarding" );

  registerWrapper( "exactJIntegration", &m_exactJIntegration ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_exactJIntegration ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Will force integration of F to have an exact integral of J." );

  registerWrapper( "maxParticleVelocity", &m_maxParticleVelocity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_maxParticleVelocity ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Velocity above which particles are deleted" );

  registerWrapper( "maxParticleVelocitySquared", &m_maxParticleVelocitySquared).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Square of the max particle velocity" );

  registerWrapper( "minParticleJacobian", &m_minParticleJacobian ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_minParticleJacobian ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Jacobian value below which particles are deleted" );

  registerWrapper( "maxParticleJacobian", &m_maxParticleJacobian ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_maxParticleJacobian ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Jacobian value above which particles are deleted" );

  registerWrapper( "overlapCorrection", &m_overlapCorrection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_overlapCorrection ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction flag (0=off, 1=increase normal force to remove gap. (not fully developed), 2=compute the SPH Jacobian and use it to scale particle density to improve the overlap. )" );

  registerWrapper( "overlapThreshold1", &m_overlapThreshold1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_overlapThreshold1 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction threshold 1" );

  registerWrapper( "overlapThreshold2", &m_overlapThreshold2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_overlapThreshold2 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction threshold 2" );

  registerWrapper( "computeSPHJacobian", &m_computeSPHJacobian ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_computeSPHJacobian ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Overlap correction flag (0=off, 1=increase normal force to remove gap. (not fully developed), 2=compute the SPH Jacobian and use it to scale particle density to improve the overlap. )" );

  registerWrapper( "shockHeating", &m_shockHeating ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_shockHeating ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable shock heating" );

  registerWrapper( "computeInternalEnergyAndTemperature", &m_computeInternalEnergyAndTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_computeInternalEnergyAndTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable computing changes to internal energy and temperature" );

  registerWrapper( "useArtificialViscosity", &m_useArtificialViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_useArtificialViscosity ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable artificial viscosity" );

  registerWrapper( "artificialViscosityQ0", &m_artificialViscosityQ0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_artificialViscosityQ0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Artificial viscosity q0 parameter" );

  registerWrapper( "cpdiDomainScaling", &m_cpdiDomainScaling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_cpdiDomainScaling ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for CPDI domain scaling" );

  registerWrapper( "subdivideParticles", &m_subdivideParticles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_subdivideParticles ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for splitting particles when they span more than a grid cell (prevents numerical fracture)" );

 registerWrapper( "computeNodalArea", &m_computeNodalArea ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_computeNodalArea ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option to compute nodal area from grid surface normals and positions (results may be inaccurate if implicit surface normals and positiosn are used)" );

 registerWrapper( "useNodePosForArea", &m_useNodePosForArea ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_useNodePosForArea ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option to compute nodal area from grid surface normals and positions (results may be inaccurate if implicit surface normals and positiosn are used)" );

  registerWrapper( "disableSurfaceNormalsAndPositionsOnCPDIScaling", &m_disableSurfaceNormalsAndPositionsOnCPDIScaling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_disableSurfaceNormalsAndPositionsOnCPDIScaling ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for disabling explicit surface normals and positions when a particle has severely deformed" );

  registerWrapper( "disableSurfaceNormalsAndPositionsOnDamage", &m_disableSurfaceNormalsAndPositionsOnDamage ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_disableSurfaceNormalsAndPositionsOnDamage ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for disabling explicit surface normals and positions when a particle has been severely damaged" );

  registerWrapper( "generalizedVortexMMS", &m_generalizedVortexMMS ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_generalizedVortexMMS ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Option for CPDI domain scaling" );

  registerWrapper( "smallMass", &m_smallMass ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_smallMass ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The small mass threshold for ignoring extremely low-mass nodes." );

  registerWrapper( "numContactGroups", &m_numContactGroups ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( m_numContactGroups ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Number of prescribed contact groups" );

  registerWrapper( "numContactFlags", &m_numContactFlags ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Number of contact flags that can appear due to damage" );

  registerWrapper( "numVelocityFields", &m_numVelocityFields ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Number of velocity fields" );

  registerWrapper( "separabilityMinDamage", &m_separabilityMinDamage ).
    setApplyDefaultValue( m_separabilityMinDamage ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Damage threshold for field separability" );

  registerWrapper( "treatFullyDamagedAsSingleField", &m_treatFullyDamagedAsSingleField ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_treatFullyDamagedAsSingleField ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Whether to consolidate fully damaged fields into a single field. Nice for modeling damaged mush." );

  registerWrapper( "surfaceDetection", &m_surfaceDetection ).
    setApplyDefaultValue( m_surfaceDetection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for automatic surface detection on the 1st cycle" );

  registerWrapper( "useSurfacePositionForContact", &m_useSurfacePositionForContact ).
    setApplyDefaultValue( m_useSurfacePositionForContact ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for to use particle mapped surface positions in contact calculations" );

  registerWrapper( "damageFieldPartitioning", &m_damageFieldPartitioning ).
    setApplyDefaultValue( m_damageFieldPartitioning ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for using the gradient of the particle damage field to partition material into separate velocity fields" );

  registerWrapper( "contactNormalType", &m_contactNormalType ).
    setApplyDefaultValue( m_contactNormalType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for contact normal type" );

  registerWrapper( "contactNormalExponent", &m_contactNormalExponent ).
    setApplyDefaultValue( m_contactNormalExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Exponent value for surface normal weights in contact normal type aligned" );

  registerWrapper( "contactGapCorrection", &m_contactGapCorrection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_contactGapCorrection ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for mitigating contact gaps" );

  registerWrapper( "resetDefGradForFullyDamagedParticles", &m_resetDefGradForFullyDamagedParticles ).
    setApplyDefaultValue( m_resetDefGradForFullyDamagedParticles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for resetting deformation gradient of fully damaged particles" );

  registerWrapper( "useReferenceVectorsForParticleUpdate", &m_useReferenceVectorsForParticleUpdate ).
    setApplyDefaultValue( m_useReferenceVectorsForParticleUpdate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for determining the integration scheme of particle vector attributes" );

  registerWrapper( "plotUnscaledParticles", &m_plotUnscaledParticles ).
    setApplyDefaultValue( m_plotUnscaledParticles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for plotting particles without CPDI domain scaling" );

  registerWrapper( "directionalOverlapCorrection", &m_directionalOverlapCorrection ).
    setApplyDefaultValue( m_directionalOverlapCorrection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for mitigating pile-up of particles at contact interfaces" );

  registerWrapper( "frictionCoefficient", &m_frictionCoefficient ).
    setApplyDefaultValue( m_frictionCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Coefficient of friction, currently assumed to be the same everywhere" );

  registerWrapper( "frictionCoefficientTable", &m_frictionCoefficientTable ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Friction coefficient table for different groups" );

  registerWrapper( "planeStrain", &m_planeStrain ).
    setApplyDefaultValue( m_planeStrain ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Flag for performing plane strain calculations" );

  registerWrapper( "numDims", &m_numDims ).
    setApplyDefaultValue( m_numDims ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "The number of active spatial dimensions, 2 for plane strain, 3 otherwise" );

  registerWrapper( "m_ijkMap", &m_ijkMap ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Map from indices in each spatial dimension to local node ID" );

  registerWrapper("numberOfSubRegions", &m_numberOfSubRegions).
    setInputFlag( InputFlags::FALSE ).
    setDefaultValue( m_numberOfSubRegions ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stores the total number of subregions" );

  registerWrapper("useEvents", &m_useEvents).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_useEvents ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Enable events" );

  registerWrapper( "computeXProfile", &m_computeXProfile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_computeXProfile ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag to enable x profiling" );

  registerWrapper( "xProfileWriteInterval", &m_xProfileWriteInterval ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_xProfileWriteInterval ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Interval between writing x profiling data" );

  registerWrapper( "nextXProfileWriteTime", &m_nextXProfileWriteTime ).
    setInputFlag( InputFlags::FALSE ).
    setDefaultValue( m_nextXProfileWriteTime ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Next x profile writing time" );

  registerWrapper( "xProfileVx0", &m_xProfileVx0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_xProfileVx0 ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Velocity to subtract from x profile" );

  registerWrapper( "elementSize", &m_hEl ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimum element size in x, y and z" );

  registerWrapper( "localMinimum", &m_xLocalMin ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local minimum" );

  registerWrapper( "localMaximum", &m_xLocalMax ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local maximum" );

    registerWrapper( "localMinimumNoGhost", &m_xLocalMinNoGhost ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local minimum without ghost cells" );

  registerWrapper( "localMaximumNoGhost", &m_xLocalMaxNoGhost ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "local maximum without ghost cells" );

  registerWrapper( "globalMinimum", &m_xGlobalMin ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "global minimum" );

  registerWrapper( "globalMaximum", &m_xGlobalMax ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "global maximum" );

  registerWrapper( "partitionExtent", &m_partitionExtent ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "parititon extent" );

  registerWrapper( "domainExtent", &m_domainExtent ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "domain extent" );

  registerWrapper( "domainF", &m_domainF ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "domain deformation gradient" );

  registerWrapper( "domainL", &m_domainL).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "domain L" );

  registerWrapper( "globalFaceReactions", &m_globalFaceReactions).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Stores the global face reactions for stress control between timesteps" );

  registerWrapper( "numElements", &m_nEl ).
    setInputFlag( InputFlags::FALSE).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "number of elements along partition directions" );

  registerWrapper( "implicitContinuumFluidPressure", &m_implicitContinuumFluidPressure ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_implicitContinuumFluidPressure ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Pressure of implicit continuum fluid" );

  m_mpmEventManager = &registerGroup< MPMEventManager >( groupKeys.mpmEventManager );
}

void SolidMechanicsMPM::postRestartInitialization()
{
  PhysicsSolverBase::postRestartInitialization();

  // Initialize friction coefficient table
  if(m_numContactGroups > 0 )
  {
    initializeFrictionCoefficients();
  }

  // Initialize box average history file and write its header
  if( m_boxAverageHistory == 1 )
  {
    // Initialize box average extent
    if(  m_boxAverageMin.size() == 0  )
    {
      m_boxAverageMin.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMin, m_xGlobalMin );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMin.size() != 3, "Box average min must have 3 elements" );
    }

    if(  m_boxAverageMax.size() == 0  )
    {
      m_boxAverageMax.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMax, m_xGlobalMax );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMax.size() != 3, "Box average min must have 3 elements" );
    }

    GEOS_ERROR_IF( ( m_boxAverageMin[0] > m_boxAverageMax[0] ) || ( m_boxAverageMin[1] > m_boxAverageMax[1] ) || ( m_boxAverageMin[2] > m_boxAverageMax[2] ) , "Box minimums must be less than box maximums");
  }
}

void SolidMechanicsMPM::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

  if( m_overlapCorrection == OverlapCorrectionOption::SPH )
  {
    GEOS_LOG_RANK_0_IF( m_computeSPHJacobian != 1, "Warning! overlapCorrection=SPH sets computeSPHJacobian=1" );
    m_computeSPHJacobian = 1;
  }

  // Activate neighbor list if necessary
  if( m_damageFieldPartitioning == 1 || 
      m_surfaceDetection == 1 || 
      m_computeSPHJacobian > 0 || 
      m_LBar > 0 || 
      m_directionalOverlapCorrection > 0 ||
      m_enableSurfaceTension == 1 )
  {
    m_needsNeighborList = 1;
  }

  // Activate nodal neighbor list if necessary
  if( m_logisticRegression == 1 )
  {
    m_needsNodalNeighborList = 1;
  }

  // Set number of active dimensions based on m_planeStrain
  m_numDims = m_planeStrain ? 2 : 3;

  // Initialize boundary condition types if they're not specified by the user
  if( m_boundaryConditionTypes.size() == 0 )
  {
    m_boundaryConditionTypes.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_boundaryConditionTypes, 0 );
  }
  else
  {
    // Throw error if boundary conditions are incorrectly specified
    GEOS_ERROR_IF( m_boundaryConditionTypes.size() != 6,
                  "boundaryConditionTypes must be of length 6. "
                  "The 6 entries correspond to BCs on the x-, x+, y-, y+, z- and z+ faces." );

    // Check elements to ensure they correspond to defined boundary condition
    // Should probably use enum for boundary conditions instead of integers (or at least unsigned, negative values shouldn't mean anything)
    GEOS_ERROR_IF( std::any_of( m_boundaryConditionTypes.begin(), m_boundaryConditionTypes.end(), []( int & bc ){ return bc < 0 || bc > 3; } ),
                  "Unknown boundary condition type specified, possible values are 0 (Outflow), 1 (Symmetry), 2 (Moving), and 3 (Contact)." );
  }

  // Initialize boundary condition types if they're not specified by the user
  if( m_enablePrescribedBoundaryTransverseVelocities.size() == 0 )
  {
    m_enablePrescribedBoundaryTransverseVelocities.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_enablePrescribedBoundaryTransverseVelocities, 0 );
  }
  else
  {
    // Throw error if boundary conditions are incorrectly specified
    GEOS_ERROR_IF( m_enablePrescribedBoundaryTransverseVelocities.size() != 6,
                   "enablePrescribedBoundaryTransverseVelocities must be of length 6. "
                   "The 6 entries correspond to transverse velocity BCs on the x-, x+, y-, y+, z- and z+ faces." );
  }

  // Initialize boundary condition types if they're not specified by the user
  if( m_boundaryFaceCoefficientsOfRestitution.size() == 0 )
  {
    m_boundaryFaceCoefficientsOfRestitution.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_boundaryFaceCoefficientsOfRestitution, 1.0 );
  }
  else
  {
    // Throw error if boundary face coefficients of restitution are incorrectly specified
    GEOS_ERROR_IF( m_boundaryFaceCoefficientsOfRestitution.size() != 6,
                  "boundaryFaceCoefficientsOfRestitution must be of length 6. "
                  "The 6 entries correspond to coefficients of restitution on the x-, x+, y-, y+, z- and z+ faces." );
  }

  // Initialize boundary condition types if they're not specified by the user
  if( m_boundaryFaceFrictionCoefficients.size() == 0 )
  {
    m_boundaryFaceFrictionCoefficients.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_boundaryFaceFrictionCoefficients, 0.0 );
  }
  else
  {
    // Throw error if boundary face friction coefficients are incorrectly specified
    GEOS_ERROR_IF( m_boundaryFaceFrictionCoefficients.size() != 6,
                  "boundaryFaceFrictionCoefficients must be of length 6. "
                  "The 6 entries correspond to friction coefficients on the x-, x+, y-, y+, z- and z+ faces." );

    GEOS_ERROR_IF( std::any_of( m_boundaryFaceFrictionCoefficients.begin(), m_boundaryFaceFrictionCoefficients.end(), []( real64 & mu ){ return mu < 0.0; } ),
                   "Boundary face coefficients of friction must be greater than 0." );
  }

    // Read and distribute BC table
  if( m_prescribedBcTable == 1 )
  {
    // Reads the FTable directly from the xml
    int numRows = m_bcTable.size( 0 );

    GEOS_ERROR_IF(numRows == 0, "Prescribed boundary conditions are enabled but no bcTable was specified.");
    for(int i = 0; i < numRows; ++i){
      GEOS_ERROR_IF(m_bcTable[i].size() != 7, "BCtable row " << i+1 << " must have 7 elements.");

      GEOS_ERROR_IF(m_bcTable[i][0] < 0, "BCTable times must be positive.");

      GEOS_ERROR_IF( !isEqual( LvArray::math::round(m_bcTable[i][1]), m_bcTable[i][1] ) ||
                     !isEqual( LvArray::math::round(m_bcTable[i][2]), m_bcTable[i][2] ) ||
                     !isEqual( LvArray::math::round(m_bcTable[i][3]), m_bcTable[i][3] ) ||
                     !isEqual( LvArray::math::round(m_bcTable[i][4]), m_bcTable[i][4] ) ||
                     !isEqual( LvArray::math::round(m_bcTable[i][5]), m_bcTable[i][5] ) ||
                     !isEqual( LvArray::math::round(m_bcTable[i][6]), m_bcTable[i][6] ), "Only integer boundary condition types are permitted.");

      GEOS_ERROR_IF( ( m_bcTable[i][1] < 0 ||
                       m_bcTable[i][2] < 0 ||
                       m_bcTable[i][3] < 0 ||
                       m_bcTable[i][4] < 0 ||
                       m_bcTable[i][5] < 0 ||
                       m_bcTable[i][6] < 0 ) &&
                     ( m_bcTable[i][1] > 3 ||
                       m_bcTable[i][2] > 3 ||
                       m_bcTable[i][3] > 3 ||
                       m_bcTable[i][4] > 3 ||
                       m_bcTable[i][5] > 3 ||
                       m_bcTable[i][6] > 3 ), "Boundary types must be 0, 1, 2 or 3");
    }
  }

  // Initialize domain F and L, then read and distribute F table
  if( m_domainF.size() == 0 )
  {
    m_domainF.resize( 3 );
    LvArray::tensorOps::fill< 3 >(m_domainF, 1.0);
  }

  if( m_domainL.size() == 0 )
  {
    m_domainL.resize( 3 );
    LvArray::tensorOps::fill< 3 >(m_domainL, 0.0);
  }

  if( m_prescribedBoundaryTransverseVelocities.size(0) == 0 )
  {
    m_prescribedBoundaryTransverseVelocities.resize(6, 2);
    LvArray::tensorOps::fill< 6, 2 >( m_prescribedBoundaryTransverseVelocities, 0.0 );
  }
  else
  {
    GEOS_ERROR_IF(!( m_prescribedBoundaryTransverseVelocities.size(0) == 6 && m_prescribedBoundaryTransverseVelocities.size(1) == 2 ), "Check dimensions of prescribedBoundaryTransverseVelocities!" );
  }

  if( m_prescribedBoundaryFTable == 1 && m_prescribedFTable == 1 )
  {
    // Reads the FTable directly from the xml
    int numRows = m_fTable.size( 0 );
    GEOS_ERROR_IF(numRows == 0, "Prescribed boundary deformation is enabled but no F table was specified.");
    for(int i = 0; i < numRows; ++i){
      GEOS_ERROR_IF(m_fTable[i].size() != 4, "F table row " << i+1 << " must have 4 elements.");
      GEOS_ERROR_IF(m_fTable[i][0] < 0, "F table times must be positive.");

      if(i == 0)
      {
        GEOS_ERROR_IF( isEqual( m_fTable[i][1], 1.0 ) ||
                       isEqual( m_fTable[i][2], 1.0 ) ||
                       isEqual( m_fTable[i][3], 1.0 ) , "Deformation of first row of F table must be 1." );
      }
      else
      {
        GEOS_ERROR_IF( m_fTable[i][1] < 0 ||
                       m_fTable[i][2] < 0 ||
                       m_fTable[i][3] < 0, "F table entries must be positive." );
      }
    }
  }

  // Throw error if boundary conditions are incorrectly specified
  GEOS_ERROR_IF( m_bodyForce.size() != 3 && m_bodyForce.size() > 0,
                 "bodyForce must be of length 3. ");

  //Initialize body force if they're not specified by the user
  if(m_bodyForce.size() == 0)
  {
    m_bodyForce.resize(3);
    LvArray::tensorOps::fill< 3 >(m_bodyForce, 0.0 );
  }
  else
  {
    GEOS_ERROR_IF( m_bodyForce.size() != 3, "Body force input must have size 3." );
  }

  // Initialze reactions for each face of box to 0
  if( m_globalFaceReactions.size() == 0 )
  {
    m_globalFaceReactions.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_globalFaceReactions, 0.0 );
  }

 // Check stress control
  if( m_boreholeStress.size() == 0 )
  {
    m_boreholeStress.resize( 6 );
    LvArray::tensorOps::fill< 6 >( m_boreholeStress, 0 );
  }
  else
  {
    GEOS_ERROR_IF( m_stressControl.size() != 6, "Borehole stress must have size 6.");
  }

  // Check stress control
  if( m_stressControl.size() == 0 )
  {
    m_stressControl.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_stressControl, 0 );
  }
  else
  {
    GEOS_ERROR_IF( m_stressControl.size() != 3, "Stress control input must have size 3.");
  }

  if( m_domainStress.size() == 0 )
  {
    m_domainStress.resize( 3 );
    LvArray::tensorOps::fill< 3 >( m_domainStress, 0.0 );
  }

  if( m_stressControlLastError.size() == 0 )
  {
    m_stressControlLastError.resize(3);
    LvArray::tensorOps::fill< 3 >( m_stressControlLastError, 0.0 );
  }

  if( m_stressControlITerm.size() == 0 )
  {
    m_stressControlITerm.resize(3);
    LvArray::tensorOps::fill< 3 >( m_stressControlITerm, 0.0 );
  }

  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    GEOS_ERROR_IF(m_stressTable.size(0) == 0, "Stress table cannot be empty if stress control is enabled");
    GEOS_ERROR_IF(m_stressTable.size(1) == 0, "Stress table must have 4 columns");

    for(int i = 0; i < m_stressTable.size(0) ; ++i)
    {
      GEOS_ERROR_IF(m_stressTable[i][0] < 0.0, "Stress table times must be positive");
    }
  }

  // Sort the grid indices and move any duplicates to the end.
  std::ptrdiff_t const numUniqueValues = LvArray::sortedArrayManipulation::makeSortedUnique( m_plottableFields.begin(),
                                                                                             m_plottableFields.end() );
  // Move unique global grid node indices to member variable
  m_plottableFieldsSorted.insert( m_plottableFields.begin(), m_plottableFields.begin() + numUniqueValues );
  GEOS_LOG_RANK_0( "Plotting fields: " << m_plottableFieldsSorted );
}

void SolidMechanicsMPM::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::mpm;

  ExecutableGroup::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                                    MeshLevel & meshLevel,
                                                    string_array const & regionNames )
  {
    ParticleManager & particleManager = meshLevel.getParticleManager();

    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

    // Set constitutive names on particles
    if( meshBody.hasParticles() )
    {
      particleManager.forParticleSubRegions< ParticleSubRegionBase >( regionNames,
                                                                      [&]( localIndex const,
                                                                           ParticleSubRegionBase & subRegion )
      {
        setParticlesConstitutiveNames( subRegion );
      } );
    }

  } );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                                    MeshLevel & meshLevel,
                                                    string_array const & regionNames )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
    if( meshBody.hasParticles() ) // Particle field registration? TODO: What goes here?
    {
      ParticleManager & particleManager = meshLevel.getParticleManager();

      particleManager.forParticleSubRegions< ParticleSubRegion >( regionNames,
                                                                  [&]( localIndex const,
                                                                       ParticleSubRegion & subRegion )
      {
        // Registration automatically allocates the 1st dimension (i.e. particle index) of the associated arrays.
        // Vector/tensor fields need to have their other dimensions specified.

        string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

        // Single-indexed fields (scalars)
        subRegion.registerField< particleMaterialType >( getName() );
        subRegion.registerField< particleMass >( getName() );
        subRegion.registerField< particleDensity >( getName() );
        subRegion.registerField< particleReferenceVolume >( getName() );
        subRegion.registerField< particleReferencePorosity >( getName() );
        subRegion.registerField< particleCrystalHealFlag >( getName() );
        subRegion.registerField< particleDeleteFlag >( getName() );
        subRegion.registerField< particleWavespeed >( getName() );
        subRegion.registerField< particleKineticEnergy >( getName() );
        subRegion.registerField< particleCohesiveZoneFlag >( getName() );
        subRegion.registerField< particleColor >( getName() );
        
        // Double-indexed fields (vectors and symmetric tensors stored in Voigt notation)
        subRegion.registerField< particleBodyForce >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleStress >( getName() ).setDimLabels( 1, voightLabels ).reference().resizeDimension< 1 >( 6 );
        subRegion.registerField< particlePlasticStrain >( getName() ).setDimLabels( 1, voightLabels ).reference().resizeDimension< 1 >( 6 );
        subRegion.registerField< particleDamageGradient >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleReferencePosition >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleReferenceMaterialDirection >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleReferenceSurfaceNormal >( getName() ).reference().resizeDimension< 1 >( 3 );

        subRegion.registerField< particleReferenceSurfacePosition >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleReferenceSurfaceTraction >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleCohesiveForce >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleReferenceMappedNodes >( getName() ).reference().resizeDimension< 1 >( 8 * subRegion.numberOfVerticesPerParticle() );
        subRegion.registerField< particleReferenceShapeFunctionValues >( getName() ).reference().resizeDimension< 1 >( 8 * subRegion.numberOfVerticesPerParticle() );
        subRegion.registerField< particleCohesiveFieldMapping >( getName() ).reference().resizeDimension< 1 >( 8 * subRegion.numberOfVerticesPerParticle() );
        
        subRegion.registerField< particleCohesiveReferenceSurfaceNormal >( getName() ).reference().resizeDimension< 1 >( 3 );
        subRegion.registerField< particleCohesiveReferenceSurfacePosition >( getName() ).reference().resizeDimension< 1 >( 3 );

        // Triple-indexed fields (vectors of vectors, non-symmetric tensors)
        subRegion.registerField< particleCohesiveReferenceDeformationGradient >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleReferenceRVectors >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleDeformationGradient >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleFDot >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleVelocityGradient >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        subRegion.registerField< particleReferenceShapeFunctionGradientValues >( getName() ).reference().resizeDimension< 1, 2 >( 8 * subRegion.numberOfVerticesPerParticle(), 3 );
        
        // Development
        subRegion.registerField< particleEstimatedSurfacePosition >( getName() ).reference().resizeDimension< 1 >( 3 );

        // Conditional registrations based on flags and other criteria, reduced memory overhead in problems that do not require these fields
        subRegion.registerField< particleReferenceTemperature >( getName() );

        if( m_cpdiDomainScaling == 1 )
        {
          subRegion.registerField< particleDomainScaledFlag >( getName() );
        }
        
        if( m_computeSPHJacobian > 0 )
        {
          subRegion.registerField< particleSPHJacobian >( getName() );
        }
        
        if( m_directionalOverlapCorrection > 0 )
        {
          subRegion.registerField< particleOverlap >( getName() );
          subRegion.registerField< particleSPHF >( getName() ).reference().resizeDimension< 1, 2 >( 3, 3 );
        }
        
        if( m_subdivideParticles == 1 )
        {
          subRegion.registerField< particleSubdivideFlag >( getName() );
          subRegion.registerField< particleCopyFlag >( getName() );
        }
        
        if( m_computeInternalEnergyAndTemperature == 1 )
        {
          subRegion.registerField< particleHeatCapacity >( getName() );
          subRegion.registerField< particleInternalEnergy >( getName() );
        }

        if( m_shockHeating == 1 || m_useArtificialViscosity == 1 || m_computeInternalEnergyAndTemperature == 1)
        {
          subRegion.registerField< particleArtificialViscosity >( getName() );
        }

        if( m_enableSurfaceTension == 1)
        {
          subRegion.registerField< particleSurfaceCurvature >( getName() );
        }

        // Adjust plotting levels for particle fields to match those specified in plottable fields
        if( m_plottableFieldsSorted.size() > 0 )
        {
          stdVector< string > const wrapperNames  = subRegion.getWrappersNames();
          for( const string & wrapperName : wrapperNames )
          {
            WrapperBase & wrapper = subRegion.getWrapperBase( wrapperName );
            if( !(m_plottableFieldsSorted.contains( wrapperName )) )
            {
              wrapper.setPlotLevel( PlotLevel::NOPLOT );
            }
            else
            {
              wrapper.setPlotLevel( PlotLevel::LEVEL_0 );
            }
          }
        }
      } );
    }
    else // Background grid field registration
    {
      NodeManager & nodeManager = meshLevel.getNodeManager();

      PlotLevel gridFieldPlotLevel = m_plotGridFields == 1 ? PlotLevel::LEVEL_0 : PlotLevel::NOPLOT;

      nodeManager.registerWrapper< array1d< int > >( viewKeyStruct::gridCohesiveNodeString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the flag for whether node is part of a cohesive zone" );

      // Begin debugging fields, currently using to check mappings of surface normals*******************************

      nodeManager.registerWrapper< array1d< real64 > >( viewKeyStruct::gridSurfaceMassString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the mass for partitioning cohesive particles" );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the mass of surface particles for each field" );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridExplicitSurfaceNormalString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface normal for partitioning cohesive particles" );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridPrincipalExplicitSurfaceNormalString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface normal corresponding to the particle with the largest ID" );

      nodeManager.registerWrapper< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds flag for whether a cohesive particle was mapped to the field of a grid node" );

      nodeManager.registerWrapper< array1d< globalIndex > >( viewKeyStruct::gridMaxMappedParticleIDString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "Holds the max global ID of the particles that mapped to the grid node" );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridCohesiveAreaString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the cohesive area for each field" );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridCohesiveForceString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the cohesive force for each field" );

      //End debugging fields *********************************************************************************************

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridReferenceSurfacePositionString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the  surface position" );

      nodeManager.registerWrapper< array1d< real64 > >( viewKeyStruct::gridReferenceMaterialVolumeString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the reference material volume" );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridReferenceAreaVectorString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the reference area vector" );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridDisplacementString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that hold the displacement on the nodes for each field");

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that hold the center of volume on the nodes for each field");

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridParticleMappedSurfaceNormalString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface normals interpoalted from particles for each field");

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridMassString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the mass on the nodes." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the material volume on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridVelocityString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current velocity on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridDVelocityString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the change in velocity from enforcing contact." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridMomentumString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current momentum on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridAccelerationString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the current acceleration on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridExternalForceString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                        " conditions as well as coupling forces such as hydraulic forces." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridSurfaceTensionForceString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the forces due to surface tension on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridInternalForceString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the internal forces on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridContactForceString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the contact force on the nodes." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridDamageString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle damage to the nodes." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle damage gradients to the nodes." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the maximum damage of any particle mapping to a given node." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightsString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds weights of surface normals at nodes." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightNormalizationString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the normalization factor for grid surface normal weights of each field." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the contact surface normals on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the contact surface positions on the nodes." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridSurfaceAreaString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface area on the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle positions to the nodes." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridNormalStressString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle normal stresses to the nodes for x profiling." );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of mapping particle mass weighted damage to the nodes for x profiling." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridVPlusString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of each XPIC and FMPM order iteration" );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridDVPlusString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the result of each XPIC and FMPM order iteration for multifield contact" );

      nodeManager.registerWrapper< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the borehole stress state for grid nodes." );

      // Development fields -------------------------------------------------------------------------------------------------

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridLRSurfaceNormalString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface normal calculated from logistic regression." );

      nodeManager.registerWrapper< array3d< real64 > >( viewKeyStruct::gridLRSurfacePositionString() ).
        setPlotLevel( gridFieldPlotLevel ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the surface position calculated from logistic regression." );

      // End development fields -------------------------------------------------------------------------------------------

      Group & nodeSets = nodeManager.sets();

      nodeSets.registerWrapper< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::WRITE_AND_READ );

      nodeSets.registerWrapper< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::WRITE_AND_READ );

      // Adjust plotting levels for grid fields to match those specified in plottable fields
      if( m_plottableFieldsSorted.size() > 0 )
      {
        stdVector< string > const wrapperNames  = nodeManager.getWrappersNames();
        for( const std::string & wrapperName : wrapperNames )
        {
          WrapperBase & wrapper = nodeManager.getWrapperBase( wrapperName );
          if( !m_plottableFieldsSorted.contains( wrapperName ) )
          {
            wrapper.setPlotLevel( PlotLevel::NOPLOT );
          }
        }
      }
    }
  } );
}

void SolidMechanicsMPM::initializePreSubGroups()
{
  PhysicsSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  Group & meshBodies = domain.getMeshBodies();

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const & meshBodyName,
                                                    MeshLevel & meshLevel,
                                                    string_array const & regionNames )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );

    if( meshBody.hasParticles() ) // Only particle regions will hold actual materials. Background grid currently holds a null material so
                                  // that the input file parser doesn't complain, but we don't need to actually do anything with it.
    {
      ParticleManager & particleManager = meshLevel.getParticleManager();
      particleManager.forParticleSubRegions< ParticleSubRegion >( regionNames, [&]( localIndex const,
                                                                                    ParticleSubRegion & subRegion )
      {
        string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
        solidMaterialName = PhysicsSolverBase::getConstitutiveName< ContinuumBase >( subRegion );
      } );
    }
  } );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = feDiscretizationManager.getGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOS_UNUSED_VAR( feDiscretization );
}


bool SolidMechanicsMPM::execute( real64 const time_n,
                                 real64 const dt,
                                 integer const cycleNumber,
                                 integer const GEOS_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  m_nextDt = solverStep( time_n,
                         dt,
                         cycleNumber,
                         domain );

  return false;
}

real64 SolidMechanicsMPM::solverStep( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  real64 dtReturn = dt;

  PhysicsSolverBase * const surfaceGenerator = this->getParent().getGroupPointer< PhysicsSolverBase >( "SurfaceGen" );

  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );

    if( surfaceGenerator!=nullptr )
    {
      surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain );
    }
  }
  else
  {
    GEOS_ERROR( "MPM solver only currently supports explicit time stepping!" );
  }

  return dtReturn;
}


void SolidMechanicsMPM::initialize( NodeManager & nodeManager,
                                    ParticleManager & particleManager,
                                    SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;
  
  // Initialize neighbor lists
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    neighborList.setParticleManager( particleManager );

    ++m_numberOfSubRegions;
  } );

  // GEOS_LOG_RANK_0( "Initialized neighbor list");

  // Get nodal position
  arrayView1d< int const > const periodic = partition.getPeriodic();
  int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridPosition = nodeManager.referencePosition();

  for(int i =0; i < 3; ++i)
  {
    if( periodic[i] && (partition.getCoords()[i] == 0 || partition.getCoords()[i] == partition.getPartitions()[i]-1) )
    {
      real64 xExtent = partition.getGlobalMax()[i] - partition.getGlobalMin()[i];
      for(int g=0; g<nodeManager.size(); ++g)
      {
        // if (gridPosition[g][i] < partition.getLocalMin()[i] && gridPosition[g][i] > partition.getLocalMax()[i] ){
          //Partition is on positive face
          if( partition.getCoords()[i] == partition.getPartitions()[i]-1) // CC: Does this need to be toleranced?
          {
            if(gridPosition[g][i] < partition.getLocalMin()[i] - xExtent/2)
            {
              gridPosition[g][i] += xExtent; //Do I nee to subtract two cells that are ghost? Shouldn't have those if periodic boundaries are on
            }
          }

          //Partition is on negative face
          if( partition.getCoords()[i] == 0){
            if(gridPosition[g][i] > partition.getLocalMax()[i] + xExtent/2) // CC: Does this need to be toleranced?
            {
              gridPosition[g][i] -= xExtent; //Do I nee to subtract two cells that are ghost? Shouldn't have those if periodic boundaries are on
            }
          }
        // }

      }
    }
  }

  // GEOS_LOG_RANK_0( "Fix periodic nodes");

  // Get local domain extent
  m_xLocalMin.resize(3);
  LvArray::tensorOps::fill< 3 >(m_xLocalMin, DBL_MAX);
  m_xLocalMax.resize(3);
  LvArray::tensorOps::fill< 3 >(m_xLocalMax, -DBL_MAX);
  for( int g=0; g<numNodes; ++g )
  {
    for( int i=0; i<3; ++i )
    {
      m_xLocalMin[i] = LvArray::math::min( m_xLocalMin[i], gridPosition[g][i] );
      m_xLocalMax[i] = LvArray::math::max( m_xLocalMax[i], gridPosition[g][i] );
    }
  }

  // m_xLocalMin0.resize( 3 );
  // LvArray::tensorOps::copy< 3 >( m_xLocalMin0, m_xLocalMin);
  // m_xLocalMax0.resize( 3 );
  // LvArray::tensorOps::copy< 3 >( m_xLocalMax0, m_xLocalMax);

  m_xLocalMinNoGhost.resize(3);
  m_xLocalMaxNoGhost.resize(3);
  m_partitionExtent.resize(3);
  for( int i=0; i<3; ++i )
  {
    m_xLocalMinNoGhost[i] = partition.getLocalMin()[i];
    m_xLocalMaxNoGhost[i] = partition.getLocalMax()[i];
    m_partitionExtent[i] = m_xLocalMax[i] - m_xLocalMin[i];
  }

  // CC: why not compute element size directly from domain extent and number of cpps across direction?
  // Get element size
  m_hEl.resize(3);
  LvArray::tensorOps::fill< 3 >(m_hEl, DBL_MAX);
  for( int g=0; g<numNodes; ++g )
  {
    for( int i=0; i<3; ++i )
    {
      real64 test = gridPosition[g][i] - m_xLocalMin[i]; // By definition, this should always be positive. Furthermore, the gridPosition
                                                         // should only be those on the local partition
      if( test > 0.0 ) // We're looking for the smallest nonzero distance from the "min" node. TODO: Could be vulnerable to a finite
                       // precision bug.
      {
        m_hEl[i] = LvArray::math::min( test, m_hEl[i] );
      }
    }
  }

  // m_hEl0.resize( 0 );
  // LvArray::tensorOps::copy< 3 >( m_hEl0, m_hEl);

  // Set SPH neighbor radius if necessary
  if( m_neighborRadius <= 0.0 )
  {
    if( m_planeStrain == 1 )
    {
      m_neighborRadius *= -1.0 * sqrt( m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] );
    }
    else
    {
      m_neighborRadius *= -1.0 * sqrt( m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] + m_hEl[2] * m_hEl[2] );
    }
  }

  // Get global domain extent excluding buffer nodes
  m_xGlobalMin.resize(3);
  m_xGlobalMax.resize(3);
  m_domainExtent.resize(3);
  for( int i=0; i<3; ++i )
  {
    m_xGlobalMin[i] = partition.getGlobalMin()[i];
    m_xGlobalMax[i] = partition.getGlobalMax()[i];
    if( !periodic[i] )
    {
      m_xGlobalMin[i] += m_hEl[i];
      m_xGlobalMax[i] -= m_hEl[i];
    }

    m_domainExtent[i] = m_xGlobalMax[i] - m_xGlobalMin[i];
  }

  // Get number of elements in each direction
  m_nEl.resize( 3 );
  for( int i=0; i<3; ++i )
  {
    m_nEl[i] = LvArray::math::round( m_partitionExtent[i] / m_hEl[i] );
  }

  // Create element map
  m_ijkMap.resize( m_nEl[0] + 1, m_nEl[1] + 1, m_nEl[2] + 1 );
  for( localIndex g=0; g<numNodes; ++g )
  {
    int i = LvArray::math::round( ( gridPosition[g][0] - m_xLocalMin[0] ) / m_hEl[0] );
    int j = LvArray::math::round( ( gridPosition[g][1] - m_xLocalMin[1] ) / m_hEl[1] );
    int k = LvArray::math::round( ( gridPosition[g][2] - m_xLocalMin[2] ) / m_hEl[2] );
    m_ijkMap[i][j][k] = g;
  }

  // GEOS_LOG_RANK_0( "Get grid information from partition");

  // Identify node sets for applying boundary conditions. We need boundary nodes and buffer nodes.
  Group & nodeSets = nodeManager.sets();
  array1d< SortedArray< localIndex > > & boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );
  boundaryNodes.resize( 6 );
  bufferNodes.resize( 6 );

  for( int face=0; face<6; face++ )
  {
    std::set< localIndex > tmpBoundaryNodes;
    std::set< localIndex > tmpBufferNodes;

    int dir0 = face / 2; // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)

    int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1
    int sign = 2 * positiveNormal - 1;
    real64 minOrMax = positiveNormal == 0 ? m_xGlobalMin[dir0] : m_xGlobalMax[dir0];

    for( localIndex g=0; g<numNodes; ++g )
    {
      real64 positionRelativeToBoundary = gridPosition[g][dir0] - minOrMax;
      real64 tolerance = 1.0e-12 * m_hEl[dir0]; // small multiple of element dimension

      if( LvArray::math::abs( positionRelativeToBoundary ) < tolerance )
      {
        tmpBoundaryNodes.insert( g );
      }

      if( sign * positionRelativeToBoundary > 0 && LvArray::math::abs( positionRelativeToBoundary ) > tolerance ) // basically a dot product with the
                                                                                                         // face normal
      {
        tmpBufferNodes.insert( g );
      }
    }

    boundaryNodes[face].insert( tmpBoundaryNodes.begin(), tmpBoundaryNodes.end() );
    bufferNodes[face].insert( tmpBufferNodes.begin(), tmpBufferNodes.end() );
  }

  initializeParticleFields( particleManager );

  
  // Initialize reaction force history file and write its header
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 && m_reactionHistory == 1 )
  {
    std::ofstream file;
    file.open( "reactionHistory.csv", std::ios::out); // | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    file << "time, F00, F11, F22, length_x, length_y, length_z, Rx-, Rx+, Ry-, Ry+, Rz-, Rz+, L00, L11, L22" << std::endl;
    file << std::setprecision( std::numeric_limits< long double >::digits10 )
         << 0.0 << ","
         << 1.0 << "," << 1.0 << "," << 1.0 << ","
         << m_domainExtent[0] << "," << m_domainExtent[1] << "," << m_domainExtent[2] << ","
         << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ","
         << 0.0 << "," << 0.0 << "," << 0.0
         << std::endl;
  }

  // Initialize box average history file and write its header
  if( m_boxAverageHistory == 1 )
  {
    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
    {
      std::ofstream file;
      file.open( "boxAverageHistory.csv", std::ios::out );
      if( file.fail() )
      {
        throw std::ios_base::failure( std::strerror( errno ) );
      }
      file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
      file << "Time, Sxx, Syy, Szz, Syz, Sxz, Sxy, Density, Damage, Internal Energy, Kinetic Energy, epxx, epyy, epzz, epyz, epxz, epxy, volume" << std::endl;
    }
    MpiWrapper::barrier( MPI_COMM_GEOS ); // wait for the header to be written

    // Initialize box average extent
    if(  m_boxAverageMin.size() == 0  )
    {
      m_boxAverageMin.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMin, m_xGlobalMin );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMin.size() != 3, "Box average min must have 3 elements" );
    }

    if(  m_boxAverageMax.size() == 0  )
    {
      m_boxAverageMax.resize( 3 );
      LvArray::tensorOps::copy< 3 >( m_boxAverageMax, m_xGlobalMax );
    }
    else
    {
      GEOS_ERROR_IF( m_boxAverageMax.size() != 3, "Box average min must have 3 elements" );
    }

    GEOS_ERROR_IF( ( m_boxAverageMin[0] > m_boxAverageMax[0] ) || ( m_boxAverageMin[1] > m_boxAverageMax[1] ) || ( m_boxAverageMin[2] > m_boxAverageMax[2] ) , "Box minimums must be less than box maximums");

    computeAndWriteBoxAverage( 0.0, 0.0, particleManager );
  }

  // GEOS_LOG_RANK_0( "Initialized files");

  // Resize grid arrays according to number of velocity fields
  localIndex maxLocalGroupNumber = 0; // Maximum contact group number on this partition.
  localIndex maxGlobalGroupNumber; // Maximum contact group number on global domain.
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< localIndex > const particleGroup = subRegion.getParticleGroup();
    for( localIndex p=0; p<subRegion.size(); ++p )
    {
      if( particleGroup[p] > maxLocalGroupNumber )
      {
        maxLocalGroupNumber = particleGroup[p];
      }
    }
  } );
  MPI_Allreduce( &maxLocalGroupNumber,
                 &maxGlobalGroupNumber,
                 1,
                 MPI_INT,
                 MPI_MAX,
                 MPI_COMM_GEOS );

  // Number of contact groups
  m_numContactGroups = maxGlobalGroupNumber + 1;

  // Specified number of damage flags.
  m_numContactFlags = m_damageFieldPartitioning == 1 ? 2 : 1;

  // Total number of velocity fields:
  m_numVelocityFields = m_numContactGroups * m_numContactFlags;

  // Resize grid field arrays
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridReferenceAreaVectorString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridReferenceSurfacePositionString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array1d< real64 > >( viewKeyStruct::gridReferenceMaterialVolumeString() ).resize( numNodes );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDisplacementString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridParticleMappedSurfaceNormalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() ).resize( numNodes, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVelocityString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceTensionForceString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridContactForceString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightsString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightNormalizationString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceAreaString() ).resize( numNodes, m_numVelocityFields);
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridPrincipalExplicitSurfaceNormalString() ).resize( numNodes, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVPlusString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVPlusString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array1d< real64 > >( viewKeyStruct::gridSurfaceMassString() ).resize( numNodes );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridExplicitSurfaceNormalString() ).resize( numNodes, 3 );
  nodeManager.getReference< array1d< int > >( viewKeyStruct::gridCohesiveNodeString() ).resize( numNodes );
  nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() ).resize( numNodes, m_numVelocityFields );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCohesiveAreaString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCohesiveForceString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array1d< globalIndex > >( viewKeyStruct::gridMaxMappedParticleIDString() ).resize( numNodes );
  nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() ).resize( numNodes, 6 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfaceNormalString() ).resize( numNodes, m_numVelocityFields, 3 );
  nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfacePositionString() ).resize( numNodes, m_numVelocityFields, 3 );

  // GEOS_LOG_RANK_0( "Initialized grid fields");

  // Precompute the square to save on computation later
  m_maxParticleVelocitySquared = m_maxParticleVelocity * m_maxParticleVelocity;

  initializeConstitutiveModelDependencies( particleManager );

  // GEOS_LOG_RANK_0( "Initialized constitutive model dependencies");

  initializeFrictionCoefficients();

  // GEOS_LOG_RANK_0( "Finished initialization");
}

localIndex SolidMechanicsMPM::partitionField( int numContactGroups,
                                              int damageFieldPartitioning,
                                              localIndex particleGroup,
                                              arraySlice1d< real64 const > const particleDamageGradient,
                                              arraySlice1d< real64 const > const particleSurfaceNormal,
                                              arraySlice1d< real64 const > const gridDamageGradient )
{
  real64 mappingNormal[3] = {};

  // By definition particle surface normals have magnitude 1
  if( isZero(LvArray::tensorOps::l2Norm< 3 >( particleSurfaceNormal ), 1e-12) )
  // if( LvArray::tensorOps::l2Norm< 3 >( particleDamageGradient ) > m_explicitSurfaceNormalInfluence / m_neighborRadius )
  {
    LvArray::tensorOps::copy< 3 >( mappingNormal, particleDamageGradient );
  }
  else
  {
    LvArray::tensorOps::copy< 3 >( mappingNormal, particleSurfaceNormal );
  }
  localIndex nodeFlag = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient, mappingNormal ) < 0.0 ) ? 1 : 0; // 0 undamaged for "A" field, 1 for "B" field
  return nodeFlag * numContactGroups + particleGroup; // This ranges from 0 to nMatFields-1
}

void SolidMechanicsMPM::initializeParticleFields( ParticleManager & particleManager )
{
  // Initialize particle fields that weren't intialized by reading the particle input file
  int const computeSPHJacobian = m_computeSPHJacobian;
  int const directionalOverlapCorrection = m_directionalOverlapCorrection;
  int const subdivideParticles = m_subdivideParticles;
  int const cpdiDomainScaling = m_cpdiDomainScaling;
  int const computeInternalEnergyAndTemperature = m_computeInternalEnergyAndTemperature;
  int const shockHeating = m_shockHeating;
  int const useArtificialViscosity = m_useArtificialViscosity;
  int const enableSurfaceTension = m_enableSurfaceTension;

  RAJA::ReduceMin< parallelDeviceReduce, real64 > minMassLocal( LvArray::NumericLimits< real64 >::max );
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );
    arrayView2d< real64 const > const constitutiveDensity = constitutiveModel.getDensity();

    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particlePorosity = subRegion.getParticlePorosity();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView2d< real64 const > const particleMaterialDirection = subRegion.getParticleMaterialDirection();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();
    arrayView1d< real64 const > const particleTemperature = subRegion.getParticleTemperature();

    arrayView1d< int > const particleMaterialType = subRegion.getField< fields::mpm::particleMaterialType >();
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 > const particleReferenceVolume = subRegion.getField< fields::mpm::particleReferenceVolume >();
    arrayView1d< real64 > const particleReferencePorosity = subRegion.getField< fields::mpm::particleReferencePorosity >();
    
    // Are these fields automatically set on initiailization?
    arrayView1d< int > const particleDeleteFlag = subRegion.getField< fields::mpm::particleDeleteFlag >();
    arrayView1d< real64 > const particleKineticEnergy = subRegion.getField< fields::mpm::particleKineticEnergy >();
    
    arrayView1d< int > const particleCrystalHealFlag = subRegion.getField< fields::mpm::particleCrystalHealFlag >();
    arrayView1d< int > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    arrayView1d< int > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    
    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView3d< real64 > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    
    arrayView2d< real64 > const particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView2d< int > const particleCohesiveFieldMapping = subRegion.getField< fields::mpm::particleCohesiveFieldMapping >();
    arrayView2d< real64 > const particleReferenceMaterialDirection = subRegion.getField< fields::mpm::particleReferenceMaterialDirection >();
    arrayView3d< real64 > const particleCohesiveReferenceDeformationGradient = subRegion.getField< fields::mpm::particleCohesiveReferenceDeformationGradient >();
    arrayView2d< real64 > const particleReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
    arrayView2d< real64 > const particleReferenceSurfacePosition = subRegion.getField< fields::mpm::particleReferenceSurfacePosition >();
    arrayView2d< real64 > const particleReferenceSurfaceTraction = subRegion.getField< fields::mpm::particleReferenceSurfaceTraction >();
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleReferenceRVectors = subRegion.getField< fields::mpm::particleReferenceRVectors >();
    
    // Development
    arrayView2d< real64 > const particleEstimatedSurfacePosition = subRegion.getField< fields::mpm::particleEstimatedSurfacePosition >();

    // Based on conditional registration of fields, requires checks of flags
    arrayView1d< int > particleDomainScaledFlag;
    if( cpdiDomainScaling == 1 )
    {
      particleDomainScaledFlag = subRegion.getField< fields::mpm::particleDomainScaledFlag >();
    }
    
    arrayView1d< real64 > particleSPHJacobian;
    if( computeSPHJacobian > 0 )
    {
      particleSPHJacobian = subRegion.getField< fields::mpm::particleSPHJacobian >();
    }
    
    arrayView3d< real64 > particleSPHF;
    if( directionalOverlapCorrection > 0 )
    {
      particleSPHF = subRegion.getField< fields::mpm::particleSPHF >();
    }
    
    arrayView1d< int > particleSubdivideFlag;
    arrayView1d< int > particleCopyFlag;
    if( subdivideParticles == 1 )
    {
      particleSubdivideFlag = subRegion.getField< fields::mpm::particleSubdivideFlag >();
      particleCopyFlag = subRegion.getField< fields::mpm::particleCopyFlag >();
    }
        
    arrayView1d< real64 > particleHeatCapacity;
    arrayView1d< real64 > particleReferenceTemperature;
    if( computeInternalEnergyAndTemperature == 1 )
    {
      particleHeatCapacity = subRegion.getField< fields::mpm::particleHeatCapacity >();
      particleReferenceTemperature = subRegion.getField< fields::mpm::particleReferenceTemperature >();
    }

    arrayView1d< real64 > particleInternalEnergy;
    if( computeInternalEnergyAndTemperature == 1 )
    {
      particleInternalEnergy = subRegion.getField< fields::mpm::particleInternalEnergy >();
    }

    arrayView1d< real64 > particleArtificialViscosity;
    if( shockHeating == 1 || useArtificialViscosity == 1 || computeInternalEnergyAndTemperature == 1 )
    {
      particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    }    

    arrayView1d< real64 > particleSurfaceCurvature;
    if( enableSurfaceTension == 1 )
    {
      particleSurfaceCurvature = subRegion.getField< fields::mpm::particleSurfaceCurvature >();
    }

    ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
    localIndex regionIndexOfSubRegion = region.getIndexInParent();
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
    {      
      particleMaterialType[p] = regionIndexOfSubRegion;

      // Initialize field from constitutive model
      particleDensity[p] = constitutiveDensity[p][0];

      particleMass[p] = particleDensity[p] * particleVolume[p];
      minMassLocal.min( particleMass[p] );

      // Initialize reference fields
      particleReferenceVolume[p] = particleVolume[p];
      particleReferencePorosity[p] = particlePorosity[p];

      // Initialize particle flags
      particleDeleteFlag[p] = 0;
      particleCrystalHealFlag[p] = 0;

      // Initiale particle fields
      particleKineticEnergy[p] = 0.0;

      if( computeInternalEnergyAndTemperature == 1 )
      {
        particleReferenceTemperature[p] = particleTemperature[p];
        particleHeatCapacity[p] = DBL_MAX; // CC: TODO Need to get this from constitutive models
      }

      // Should already be initialized by default value from DECLARE_FIELD, these might be redundant
      if( computeInternalEnergyAndTemperature == 1 )
      {
        particleInternalEnergy[p] = 0.0;
      }

      if( shockHeating == 1 || useArtificialViscosity == 1 || computeInternalEnergyAndTemperature == 1 )
      {
        particleArtificialViscosity[p] = 0.0;
      }

      if( computeSPHJacobian > 1 )
      {
        particleSPHJacobian[p] = 1.0;
      }
      
      particleCohesiveZoneFlag[p] = 0;

      if( subdivideParticles == 1 )
      {
        particleSubdivideFlag[p] = 0;
        particleCopyFlag[p] = -1;
      }

      if( cpdiDomainScaling == 1 )
      {
        particleDomainScaledFlag[p] = 0;
      }

      if( enableSurfaceTension == 1 )
      {
        particleSurfaceCurvature[p] = 0.0;
      }

      particleColor[p] = 0;

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        particleCohesiveFieldMapping[p][g] = particleGroup[p];
      }

      LvArray::tensorOps::fill< 3 >( particleBodyForce[p], 0.0 );
      LvArray::tensorOps::fill< 3 >( particleCohesiveForce[p], 0.0);
      LvArray::tensorOps::fill< 3 >( particleEstimatedSurfacePosition[p], 0.0 );

      LvArray::tensorOps::copy< 3 >( particleReferencePosition[p], particlePosition[p] );
      LvArray::tensorOps::copy< 3 >( particleReferenceMaterialDirection[p], particleMaterialDirection[p] );
      LvArray::tensorOps::copy< 3, 3 >( particleReferenceRVectors[p], particleRVectors[p]);

      LvArray::tensorOps::copy< 3 >( particleReferenceSurfaceNormal[p], particleSurfaceNormal[p] );
      LvArray::tensorOps::copy< 3 >( particleReferenceSurfacePosition[p], particleSurfacePosition[p]);
      LvArray::tensorOps::copy< 3 >( particleReferenceSurfaceTraction[p], particleSurfaceTraction[p]);

      // Initialize deformation gradient and velocity gradient
      LvArray::tensorOps::fill< 3, 3 >( particleFDot[p], 0.0 );
      LvArray::tensorOps::fill< 3, 3 >( particleVelocityGradient[p], 0.0 );

      for( localIndex i = 0; i < 3; ++i )
      {
        for( localIndex j = 0; j < 3; ++j )
        {
          particleDeformationGradient[p][i][j] = i == j ? 1.0 : 0.0;
          particleCohesiveReferenceDeformationGradient[p][i][j] = i == j ? 1.0 : 0.0; // Need to check if this is still needed
          if( directionalOverlapCorrection > 0 )
          {
            particleSPHF[p][i][j] = i == j ? 1.0 : 0.0;
          }
        }
      }

      // LvArray::tensorOps::fill< 6 >( particlePlasticStrain[p], 0.0 );
    } );
  } );

  // Set small mass threshold
  m_smallMass = fmin( MpiWrapper::min( minMassLocal.get() ) * 1.0e-12, m_smallMass );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.setActiveParticleIndices(); // Needed for computeAndWriteBoxAverage().
  } );
}

// Load any particle dependent properties into constitutive model such as strength scale which is determiend per particle not by material model input
void SolidMechanicsMPM::initializeConstitutiveModelDependencies( ParticleManager & particleManager )
{
  // Load strength scale into constitutive model (for ceramic)
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );
    if( constitutiveModel.hasWrapper( "strengthScale" ) )
    {
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      arrayView1d< real64 > const particleStrengthScale = subRegion.getParticleStrengthScale();
      arrayView1d< real64 > const constitutiveStrengthScale = constitutiveModel.getReference< array1d< real64 > >( "strengthScale" );
      forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveStrengthScale[p] = particleStrengthScale[p];
      } );
    }
  } );
}

real64 SolidMechanicsMPM::explicitStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Get spatial partition, get node and particle managers." );
  solverProfiling( "Get spatial partition, get node and particle managers." );

  // SpatialPartition & partition = dynamic_cast< SpatialPartition & >( domain.getReference< PartitionBase >( keys::partitionManager ) );
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >( domain.getGroup( domain.groupKeys.partitionManager ) );
  arrayView1d< int const > const periodic = partition.getPeriodic();

  // ***** We assume that there are exactly two mesh bodies, and that one has particles and one does not. *****
  Group & meshBodies = domain.getMeshBodies();

  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >( 0 );
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >( 1 );
  MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
  MeshBody & grid = !meshBody1.hasParticles() ? meshBody1 : meshBody2;

  ParticleManager & particleManager = particles.getBaseDiscretization().getParticleManager();
  MeshLevel & mesh = grid.getBaseDiscretization();
  NodeManager & nodeManager = mesh.getNodeManager();
  //#######################################################################################

  //#######################################################################################
  if( cycleNumber == 0 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "At time step zero, perform initialization calculations" );
    solverProfiling( "At time step zero, perform initialization calculations");
    initialize( nodeManager, particleManager, partition );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Set grid multi-field labels to avoid a VTK output bug" );
  solverProfiling( "Set grid multi-field labels to avoid a VTK output bug" );
  // Must be done every time step despite grid fields being registered
  // TODO: Only doing this on the 1st cycle breaks restarts, can we do better? Why aren't labels part of the restart data?
  setGridFieldLabels( nodeManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Set grid fields to zero" );
  solverProfiling( "Set grid fields to zero" );
  initializeGridFields( nodeManager );
  //#######################################################################################


  //#######################################################################################
  if( m_useEvents == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Check if event has been triggered" );
    solverProfiling( "Check if event has been triggered");
    triggerEvents( dt,
                   time_n,
                   particleManager,
                   partition );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_subdivideParticles == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Subdivide highly deformed particles" );
    solverProfiling( "Subdivide highly deformed particles");
    subdivideParticles( particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update global-to-local map" );
  solverProfiling( "Update global-to-local map" );
  particleManager.updateMaps();
  //#######################################################################################


  //#######################################################################################
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) > 1 && m_needsNeighborList == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines,  "Perform particle ghosting" );
    solverProfiling( "Perform particle ghosting" );

    // Move everything into host memory space
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      subRegion.forWrappers( [&]( WrapperBase & wrapper )
      {
        wrapper.move( LvArray::MemorySpace::host, true );
      } );
    } );

    // Perform ghosting
    partition.getGhostParticlesFromNeighboringPartitions( domain, m_neighborRadius );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Get indices of non-ghost particles" );
  solverProfiling( "Get indices of non-ghost particles" );
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.setActiveParticleIndices();
  } );
  //#######################################################################################


  //#######################################################################################
  if( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Correct ghost particle centers across periodic boundaries" );
    solverProfiling( "Correct ghost particle centers across periodic boundaries" );
    correctGhostParticleCentersAcrossPeriodicBoundaries(particleManager, partition);
  }
  //#######################################################################################


  //#######################################################################################
  if( m_needsNeighborList == 1 )
  { // This optimization compares neighbor list creation time for different
    // bin sizes, and finds an optimum.  It can be non deterministic and
    // is currently disabled to ensure integrated tests pass.
    // Optimize
    // if( cycleNumber == 0 )
    // {
    //   optimizeBinSort( particleManager );
    // }
    // else
    // {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Construct neighbor list" );
    solverProfiling( "Construct neighbor list" );
    (void) computeNeighborList( particleManager );
    // }
  }
  //#######################################################################################


  //#######################################################################################
  if( m_surfaceDetection > 0 && ( cycleNumber == 0 || m_surfaceHealing == 1 ) )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute surface flags" );
    solverProfiling( "Compute surface flags" );
    m_surfaceHealing = false;
    computeSurfaceFlags( particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  // GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Perform 1D logistic regression to find normals and positions" );
  // solverProfiling( "Perform 1D logistic regression to find normals and positions" );
  // logisticRegression1D( particleManager,
  //                      nodeManager );
  //#######################################################################################


  //#######################################################################################
  if( m_prescribedBcTable == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update BCs based on bcTable" );
    solverProfiling( "Update BCs based on bcTable" );
    boundaryConditionUpdate( dt, time_n );
  }

  // Update other boundary condition dependent solver flags
  m_hasContact = m_numVelocityFields > 1 || std::any_of(m_boundaryConditionTypes.begin(), m_boundaryConditionTypes.end(),[]( int & bc ){ return bc == 3; });
  //#######################################################################################


  //#######################################################################################
  if( m_disableSurfaceNormalsAndPositionsOnDamage == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Disable explicit normals and positions of fully damaged particles" );
    solverProfiling( "Disable explicit normals and positions of fully damaged particles" );

    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      arrayView1d< real64 > const particleDamage = subRegion.getParticleDamage();
      arrayView2d< real64 > const particleReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
      arrayView2d< real64 > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
      arrayView2d< real64 > const particleReferenceSurfacePosition = subRegion.getField< fields::mpm::particleReferenceSurfacePosition >();
      arrayView2d< real64 > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
      // arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        if( particleDamage[p] >= 0.9999 ) // Should this threshold be an optional user input
        {
          LvArray::tensorOps::fill< 3 >( particleSurfaceNormal[p], 0.0 );
          LvArray::tensorOps::fill< 3 >( particleSurfacePosition[p], 0.0 );

          LvArray::tensorOps::fill<3>( particleReferenceSurfaceNormal[p], 0.0 );
          LvArray::tensorOps::fill<3>( particleReferenceSurfacePosition[p], 0.0 );

          // if( particleSurfaceFlag[p] > 1 )
          // {
          //   particleSurfaceFlag[p] = 1;
          // }
        }
      } );
    } );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_cpdiDomainScaling == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Perform r-vector scaling (CPDI domain scaling)" );
    solverProfiling( "Perform r-vector scaling (CPDI domain scaling)" );
    cpdiDomainScaling( particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Resize and populate mapping arrays" );
  solverProfiling( "Resize and populate mapping arrays" );
  resizeMappingArrays( particleManager );
  populateMappingArrays( particleManager, nodeManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Create nodal neighborlist" );
  solverProfiling( "Create nodal neighborlist" );
  if( m_needsNodalNeighborList == 1 )
  {
    generateNodalNeighborList( particleManager,
                               nodeManager );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_damageFieldPartitioning == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute damage field gradient" );
    solverProfiling( "Compute damage field gradient" );
    computeDamageFieldGradient( particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update surface flag overload" );
  solverProfiling( "Update surface flag overload" );
  // We now read in surface flags so don't need to update them based on damage.
  // Keeping this here as it will eventually be used to support other surface
  // flag usages.
  // updateSurfaceFlagOverload( particleManager );
  //#######################################################################################


  //#######################################################################################
  if( m_computeSPHJacobian > 0 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute SPH Jacobian for overlap correction" );
    solverProfiling( "Compute SPH Jacobian for overlap correction" );
    computeSPHJacobian( particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_damageFieldPartitioning == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Project damage field gradient to the grid and then sync" );
    solverProfiling( "Project damage field gradient to the grid and then sync" );
    projectDamageFieldGradientToGrid( domain,
                                      particleManager,
                                      nodeManager,
                                      mesh );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute particle field mappings" );
  solverProfiling( "Compute particle field mappings" );
  computeParticleFieldMappings( particleManager,
                                nodeManager );
  //#######################################################################################


  //#######################################################################################
  if( m_computeParticleSurfaceNormalsAndPositions == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute particle surface normals and positions" );
    solverProfiling( "Compute particle surface normals and positions" );
    m_computeParticleSurfaceNormalsAndPositions = 0;
    computeParticleSurfaceNormalsAndPositions( particleManager,
                                               nodeManager );
  }
  //#######################################################################################
  

  //#######################################################################################
  if( m_reinitializeCohesiveZones == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Reinitialize cohesive zones" );
    solverProfiling( "Reinitialize cohesive zones" );


  }
  //#######################################################################################


  //#######################################################################################
  if( m_enableCohesiveLaws == 1 )
  {
    if( m_referenceCohesiveZone == 1 )
    {
      GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Initialize cohesive zones" );
      solverProfiling( "Initialize cohesive zones" );
      projectParticleSurfaceNormalsToGrid( domain, particleManager, nodeManager, mesh );

      initializeCohesiveReferenceConfiguration( domain,
                                                particleManager,
                                                nodeManager,
                                                mesh );

      m_referenceCohesiveZone = 0;
    }

    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute reactions due to cohesive laws" );
    solverProfiling( "Compute reactions due to cohesive laws" );
    enforceCohesiveLaw( particleManager,
                        nodeManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute particle body forces" );
  solverProfiling( "Compute particle body forces" );
  computeBodyForce( particleManager);
  //#######################################################################################


  //#######################################################################################
  if( m_enableSurfaceTension == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute sph surface curvature" );
    solverProfiling( "Compute sph surface curvature" );
    computeSPHSurfaceCurvature( particleManager );
  }
  //#######################################################################################

  //#######################################################################################
  if( m_generalizedVortexMMS == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute particle body forces" );
    solverProfiling( "Compute particle body forces" );
    computeGeneralizedVortexMMSBodyForce( time_n,
                                          particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_enableBoreholePressure == 1 ){
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update borehole stress state on grid nodes" );
    solverProfiling( "Update borehole stress state on grid nodes" );
    updateGridBoreholeStress( nodeManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Particle color sort" );
  solverProfiling( "Particle color sort" );
  particleColorSort( particleManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Particle-to-grid interpolation" );
  solverProfiling( "Particle-to-grid interpolation" );
  // Other schemes than Atomics are deprecated
  switch( m_gpuScheme )
  {
    case GPUSchemeOption::Atomics:
      particleToGrid( time_n,
                      cycleNumber,
                      particleManager,
                      nodeManager );
      break;
    case GPUSchemeOption::NoAtomics:
      particleToGrid_noAtomics( time_n,
                                cycleNumber,
                                particleManager,
                                nodeManager );
      break;
    case GPUSchemeOption::RandomMix:
      particleToGrid_randomMix( time_n,
                                cycleNumber,
                                particleManager,
                                nodeManager );
      break;
    case GPUSchemeOption::MinimalAtomics:
      particleToGrid_minimalAtomics( time_n,
                                    cycleNumber,
                                    particleManager,
                                    nodeManager );
      break;
    case GPUSchemeOption::Reduction:
      particleToGrid_reduction( time_n,
                                cycleNumber,
                                particleManager,
                                nodeManager );
      break;
    case GPUSchemeOption::Colors:
      particleToGrid_colors( time_n,
                             cycleNumber,
                             particleManager,
                             nodeManager );
      break;
  }

  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Grid MPI operations" );
  solverProfiling( "Grid MPI operations" );
  stdVector< std::string > fieldNames1 = { viewKeyStruct::gridMassString(),
                                             viewKeyStruct::gridDamageString(),
                                             viewKeyStruct::gridMaterialVolumeString(),
                                             viewKeyStruct::gridMomentumString(),
                                             viewKeyStruct::gridCenterOfMassString(),
                                             viewKeyStruct::gridInternalForceString(),
                                             viewKeyStruct::gridExternalForceString() };
  if( m_hasContact == 1 )
  {
    fieldNames1.push_back( viewKeyStruct::gridSurfaceFieldMassString() );
    fieldNames1.push_back( viewKeyStruct::gridSurfaceNormalString() );
    fieldNames1.push_back( viewKeyStruct::gridSurfacePositionString() );
  }
  syncGridFields( fieldNames1, domain, nodeManager, mesh, MPI_SUM );
  stdVector< std::string > fieldNames2 = { viewKeyStruct::gridMaxDamageString() };
  syncGridFields( fieldNames2, domain, nodeManager, mesh, MPI_MAX );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Determine trial momenta and velocities based on acceleration due to internal and external forces, but before contact enforcement" );
  solverProfiling( "Determine trial momenta and velocities based on acceleration due to internal and external forces, but before contact enforcement" );
  gridTrialUpdate( dt, nodeManager );
  //#######################################################################################


  //#######################################################################################
  if( m_hasContact == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Contact enforcement" );
    solverProfiling( "Contact enforcement" );
    enforceContact( dt, 
                    // domain,
                    // particleManager,
                    nodeManager
                    // mesh
                  );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_computeNodalArea == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute nodal surface area" );
    solverProfiling( "Compute nodal surface area" );
    
    computeNodalAreas( nodeManager );
  }

  //#######################################################################################


  //#######################################################################################
  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Interpolate stress table" );
    solverProfiling( "Interpolate stress table" );
    interpolateStressTable( dt, time_n );
  }
  //#######################################################################################


  //#######################################################################################
  if( !( m_stressControl[0] == 1 && m_stressControl[1] == 1 && m_stressControl[2] == 1 ) && ( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 ) )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Interpolate F table" );
    solverProfiling( "Interpolate F table" );
    interpolateFTable( dt, time_n );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Apply essential boundary conditions" );
  solverProfiling( "Apply essential boundary conditions" );
  applyEssentialBCs( dt, time_n, nodeManager );
  //#######################################################################################


  //#######################################################################################
  if( m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) ) )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute local averages to generate profile in x-direction: write to file." );
    solverProfiling( "Compute local averages to generate profile in x-direction: write to file." );
    computeXProfile( cycleNumber,
                     time_n,
                     dt,
                     nodeManager,
                     partition );
    m_nextXProfileWriteTime += m_xProfileWriteInterval;
  }
  //#######################################################################################


  // //#######################################################################################
  // TODO: Figure out syncing on ghost particles and move this to after the F update
  // if( m_directionalOverlapCorrection > 0 )
  // {
  //   GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Directional overlap correction" );
  //   solverProfiling( "Directional overlap correction" );
  //   directionalOverlapCorrection( dt, particleManager );
  // }
  // //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Grid-to-particle interpolation" );
  solverProfiling( "Grid-to-particle interpolation" );
  gridToParticle( dt, particleManager, nodeManager, domain, mesh );
  //#######################################################################################


  //#######################################################################################
  if( m_prescribedFTable == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update particle positions according to prescribed F Table" );
    solverProfiling( "Update particle positions according to prescribed F Table" );
    applySuperimposedVelocityGradient( dt,
                                       particleManager,
                                       partition );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update deformation gradient" );
  solverProfiling( "Update deformation gradient" );
  updateDeformationGradient( dt, particleManager );
  //#######################################################################################


  //#######################################################################################
  // Scale Jacobian to prevent overdensification if overlap correction type 2 is used.
  if( m_overlapCorrection == OverlapCorrectionOption::SPH )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Scale F based on SPH-J" );
    solverProfiling( "Scale F based on SPH-J" );
    sphOverlapCorrection( dt,
                       particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update particle geometry (e.g. volume, r-vectors) and density" );
  solverProfiling( "Update particle geometry (e.g. volume, r-vectors) and density" );
  particleKinematicUpdate( dt,
                           particleManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute kinetic energy" );
  solverProfiling( "Compute kinetic energy" );
  // 1/2*m*v^2 calculation, TODO add micro kinetic energy.
  computeKineticEnergy( particleManager );
  //#######################################################################################


  //#######################################################################################
  if( m_useArtificialViscosity == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute artificial viscosity" );
    solverProfiling( "Compute artificial viscosity" );
    computeArtificialViscosity( particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_computeInternalEnergyAndTemperature == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Increment internalEnergy with old stress and Fdot" );
    solverProfiling( "Increment internalEnergy with old stress and Fdot" );
    computeInternalEnergyAndTemperature( dt,
                                        particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update constitutive model dependencies" );
  solverProfiling( "Update constitutive model dependencies" );
  updateConstitutiveModelDependencies( particleManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update stress" );
  solverProfiling( "Update stress" );
  updateStress( dt, particleManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update solver dependencies" );
  solverProfiling( "Update solver dependencies" );
  updateSolverDependencies( particleManager );
  //#######################################################################################


  //#######################################################################################
  if( m_computeInternalEnergyAndTemperature == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Increment internalEnergy with new stress and Fdot" );
    solverProfiling( "Increment internalEnergy with new stress and Fdot" );
    computeInternalEnergyAndTemperature( dt,
                                        particleManager );
  }
  //#######################################################################################


  //#######################################################################################
  if( m_boxAverageHistory == 1 && time_n + dt >= m_nextBoxAverageWriteTime )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Compute and write box averages" );
    solverProfiling( "Compute and write box averages" );
    computeAndWriteBoxAverage( time_n, dt, particleManager );
    m_nextBoxAverageWriteTime += m_boxAverageWriteInterval;
  }
  //#######################################################################################


  //#######################################################################################
  if( m_writeParticleData == 1 && time_n + dt >= m_nextParticleDataWriteTime )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Write particle data to file" );
    solverProfiling( "Write particle data to file" );
    writeParticleData( time_n, particleManager );
    m_nextParticleDataWriteTime += m_particleDataWriteInterval;
  }
  //#######################################################################################


  //#######################################################################################
  // stress control reads a principal stress table.  we compute the
  // difference between most recent box sum and the prescribed value and increment
  // the domain strain rate accordingly.  This is a simple control loop.
  if( m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Stress control" );
    solverProfiling( "Stress control" );
    stressControl( dt,
                   particleManager,
                   partition );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Calculate stable time step" );
  solverProfiling( "Calculate stable time step" );
  real64 dtReturn = getStableTimeStep( particleManager );
  //#######################################################################################


  // TODO check if this should be moved in front of flag out of range particles
  //#######################################################################################
  if( m_prescribedBoundaryFTable == 1 || m_prescribedFTable == 1 || m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Resize grid based on F-table" );
    solverProfiling( "Resize grid based on F-table" );
    resizeGrid( partition, nodeManager, dt );
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Update global-to-local map" );
  solverProfiling( "Update global-to-local map" );
  flagOutOfRangeParticles( particleManager );
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Delete bad particles" );
  solverProfiling( "Delete bad particles" );
  deleteBadParticles( particleManager );
  //#######################################################################################


  //#######################################################################################
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) > 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Particle repartitioning" );
    solverProfiling( "Particle repartitioning" );

    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      subRegion.forWrappers( [&]( WrapperBase & wrapper )
      {
        wrapper.move( LvArray::MemorySpace::host, true );
      } );
      partition.repartitionMasterParticles( domain, subRegion );

      subRegion.setActiveParticleIndices(); // Can't we just loop over all the particles when correcting across periodic boundaries? Do we really need to reset the active particle indices again?
    } );

    // Correct particle centers across periodic boundaries
    if( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 )
    {
        correctParticleCentersAcrossPeriodicBoundaries(particleManager, partition);
    }
  }
  //#######################################################################################


  //#######################################################################################
  // Option to set F for fully damaged particles to J^(1/3)*[I], in 3D, so we don't get
  // negative J for super sheared particles with finite precision F Update.  This
  // Should only be used when all materials in the domain have a hypo-elastic
  // deviatoric update.
  if( m_resetDefGradForFullyDamagedParticles == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Set F to scaled value for damaged particles, maintain J" );
    solverProfiling( "Set F to scaled value for damaged particles, maintain J" );
  }
  resetDeformationGradient( particleManager );
  //#######################################################################################


  //#######################################################################################
  if ( m_plotUnscaledParticles == 1 )
  {
    GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "Reset CPDI R-vectors for plotting if option is selected to plot unscaled domains" );
    solverProfiling( "Reset CPDI R-vectors for plotting if option is selected to plot unscaled domains" );
    unscaleCPDIVectors( particleManager);
  }
  //#######################################################################################


  //#######################################################################################
  GEOS_LOG_LEVEL_BY_RANK( logInfo::MPMSubroutines, "End of explicitStep");
  solverProfiling( "End of explicitStep" );
  if( m_solverProfiling >= 1 )
  {
    printProfilingResults();
  }
  //#######################################################################################


  // Return stable time step
  return dtReturn;
}


void SolidMechanicsMPM::triggerEvents( const real64 dt,
                                       const real64 time_n,
                                       ParticleManager & particleManager,
                                       SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  int const planeStrain = m_planeStrain;

  //Iterate over every MPM Events to check if conditions are met and if so perform event
  m_mpmEventManager->forSubGroups< MPMEventBase >( [&]( MPMEventBase & event )
  {
    real64 eventTime = event.getTime();
    real64 eventInterval = event.getInterval();

    if( ( eventTime - dt / 2 <= time_n && time_n <= eventTime + eventInterval + dt/2 ) && !event.isComplete() )
    {
      if( event.getName() == "MaterialSwap" )
      {
        MaterialSwapMPMEvent & materialSwap = dynamicCast< MaterialSwapMPMEvent & >( event );
        GEOS_LOG_RANK_0("Performing material swap");
        performMaterialSwap( particleManager, materialSwap.getSourceRegion(), materialSwap.getDestinationRegion() );
        event.setIsComplete( 1 );
      }

      if( event.getName() == "BoreholePressure" )
      {
        BoreholePressureMPMEvent & boreholePressureEvent = dynamicCast< BoreholePressureMPMEvent & >( event );
        GEOS_LOG_RANK_0("Setting borehole pressure");

        m_enableBoreholePressure = 1;
        m_boreholeRadius = boreholePressureEvent.getBoreholeRadius();
        real64 startPressure = boreholePressureEvent.getStartPressure();
        real64 endPressure = boreholePressureEvent.getEndPressure();
        SolidMechanicsMPM::InterpolationOption interpolationType = static_cast< InterpolationOption >( boreholePressureEvent.getInterpType() ); // Interpolation type should be available in events

        // Set boreholePressure to interpolated value.  The default is 0, but at the end of the event
        // it won't be reset.

        real64 boreholePressure = 0.0;
        interpolateValueInRange( time_n,
                                eventTime,
                                eventTime + eventInterval,
                                startPressure,
                                endPressure,
                                boreholePressure, // output, overwritten from interpolaiton.
                                interpolationType );

        m_boreholeStress[0] = -boreholePressure;
        m_boreholeStress[1] = -boreholePressure;
        m_boreholeStress[2] = -boreholePressure;
        m_boreholeStress[3] = 0.0;
        m_boreholeStress[4] = 0.0;
        m_boreholeStress[5] = 0.0;

        //sevent.setIsComplete( 1 );
      }

      if( event.getName() == "Anneal" )
      {
        AnnealMPMEvent & anneal = dynamicCast< AnnealMPMEvent & >( event );

        particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region )
        {
          if( region.getName() == anneal.getTargetRegion() || anneal.getTargetRegion() == "all" )
          {
            auto & targetSubRegions = region.getSubRegions();
            for( int r=0; r < targetSubRegions.size(); ++r)
            {
              ParticleSubRegion & targetSubRegion = dynamicCast< ParticleSubRegion & >( *targetSubRegions[r] );

              // Get constitutive model reference
              string const & solidMaterialName = targetSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
              ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( targetSubRegion, solidMaterialName );

              GEOS_ERROR_IF( !constitutiveModel.hasWrapper( "oldStress" ), "Cannot anneal constitutive model that does not have oldStress wrapper!");
              arrayView3d< real64 > const constitutiveOldStress = constitutiveModel.getReference< array3d< real64 > >( "oldStress" );

              // SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
              forAll< serialPolicy >( constitutiveOldStress.size(0), [=] GEOS_HOST_DEVICE ( localIndex const p )
              {
                  real64 knockdown = 0.0; // At the end of the annealing process, strongly enforce zero deviatoric stress.
                  if( !( time_n - dt / 2 <= eventTime + eventInterval && time_n + dt / 2 > eventTime + eventInterval ) )
                  {
                    // Otherwise, scale it down by a number that gradually goes from 1 to 0
                    knockdown = fmax( 0.0, 1.0 - dt * 20.0 * ( time_n - eventTime ) / ( eventInterval * eventInterval ) ); // Gaussian decay
                  }

                  // Smoothly knock down the deviatoric stress to simulate annealing
                  real64 stress[6];
                  LvArray::tensorOps::copy< 6 >( stress, constitutiveOldStress[p][0]);

                  real64 trialP;
                  real64 trialQ;
                  real64 deviator[6];
                  twoInvariant::stressDecomposition( stress,
                                                    trialP,
                                                    trialQ,
                                                    deviator );

                  LvArray::tensorOps::copy< 6 >( constitutiveOldStress[p][0], stress );
                  twoInvariant::stressRecomposition( trialP,
                                                    knockdown * trialQ,
                                                    deviator,
                                                    stress );
                  LvArray::tensorOps::copy< 6 >( constitutiveOldStress[p][0], stress );
              } );
            }
          }
        });
      }

      if( event.getName() == "Heal" )
      {
        // HealMPMEvent & heal = dynamicCast< HealMPMEvent & >( event );

        m_surfaceHealing = true;
        particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
        {
          arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();

          SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
          {
            localIndex const p = activeParticleIndices[pp];
            if( particleSurfaceFlag[p] < 3 ) //  3 and 4 denote cohesive surface flags (still a binder interface)
            {
              particleSurfaceFlag[p] = 0;
            }
          });
        });
        event.setIsComplete( 1 );
      }

      if( event.getName() == "InsertPeriodicContactSurfaces" )
      {
          // InsertPeriodicContactSurfacesMPMEvent & insertPeriodicContactSurfaces = dynamicCast< InsertPeriodicContactSurfacesMPMEvent & >( event );

          real64 hEl[3] = {};
          LvArray::tensorOps::copy< 3 >(hEl, m_hEl);
          real64 xGlobalMin[3] = {};
          LvArray::tensorOps::copy< 3 >(xGlobalMin, m_xGlobalMin);
          real64 xGlobalMax[3] = {};
          LvArray::tensorOps::copy< 3 >(xGlobalMax, m_xGlobalMax);

          // Perhaps I need to find the largest particle size during initialize
          auto periodic = partition.getPeriodic();

          particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
          {
            arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
            arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();

            SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
            forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
            {
              localIndex const p = activeParticleIndices[pp];
              for(int i =0; i < 3; ++i)
              {
                if( periodic[i] == 1 && ( ( particlePosition[p][i] - xGlobalMin[i] < hEl[i] ) || (  xGlobalMax[i] - particlePosition[p][i] < hEl[i] ) ) )
                {
                  particleSurfaceFlag[p] = 2;
                  break;
                }
              }
            });
          });
          event.setIsComplete( 1 );
      }

      if( event.getName() == "CrystalHeal" ){
        CrystalHealMPMEvent & crystalHeal = dynamicCast< CrystalHealMPMEvent & >( event );

        int healType = crystalHeal.getHealType();
        particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion & region )
        {
          if( region.getName() == crystalHeal.getTargetRegion() || crystalHeal.getTargetRegion() == "all" )
          {
            // Copy particle data from source sub region to destination sub region
            subGroupMap & targetSubRegions = region.getSubRegions();
            for( int r=0; r < targetSubRegions.size(); ++r)
            {
              ParticleSubRegion & targetSubRegion = dynamicCast< ParticleSubRegion & >( *targetSubRegions[r] );

              arrayView1d< real64 > const particleDamage = targetSubRegion.getParticleDamage();
              arrayView1d< int > const particleCrystalHealFlag = targetSubRegion.getField< fields::mpm::particleCrystalHealFlag >();

              SortedArrayView< localIndex const > const activeParticleIndices = targetSubRegion.activeParticleIndices();
              string const & solidMaterialName = targetSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
              ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( targetSubRegion, solidMaterialName );
              if( !crystalHeal.getMarkedParticlesToHeal() )
              {
                arrayView2d< real64 > const particleStress = targetSubRegion.getField< fields::mpm::particleStress >();
                arrayView3d< real64 > const particleRVectors = targetSubRegion.getParticleRVectors();
                arrayView3d< real64 > const particleDeformationGradient = targetSubRegion.getField< fields::mpm::particleDeformationGradient >();
                arrayView1d< real64 > const particleReferenceVolume = targetSubRegion.getField< fields::mpm::particleReferenceVolume >();

                arrayView2d< real64 const > const particleDamageGradient = targetSubRegion.getField< fields::mpm::particleDamageGradient >();

                // CC: TODO Need to add this for Mike
                // arrayView1d< real64 const > const particleReferencePorosity = targetSubRegion.getField< fields::mpm::particleReferencePorosity >();

                // Scale constitutive model crackspeed so damage does not evolve while healing
                if( constitutiveModel.hasWrapper( "crackSpeed" ) )
                {
                  real64 & crackSpeed = constitutiveModel.getReference< real64 >( "crackSpeed" );
                  crackSpeed *= 1e-100;
                }

                forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
                {
                  localIndex const p = activeParticleIndices[pp];

                  real64 temp[3] = {};
                  // real64 particleDamageGradientNormalized[3];
                  // LvArray::tensorOps::copy< 3 >( particleDamageGradientNormalized, particleDamageGradient[p] );
                  // norm = LvArray::tensorOps::l2Norm< 3 >( particleDamageGradientNormalized );
                  // // CC: zero stress component in ZZ for plane strain?

                  // // Does the damage gradient need to be normalized since we only care about the sign of the stress magnitude along direction?
                  // real64 normalStress = 0.0
                  // if(norm > 1e-20 ){
                  //   LvArray::tensorOps::normalize< 3 >( particleDamageGradientNormalized ); // CC: debug ( should probably be normalized compute of damage field for particles )
                  //   LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, particleStress[p], particleDamageGradientNormalized );
                  //   normalStress = LvArray::tensorOps::AiBi< 3 >( particleDamageGradientNormalized, temp );
                  // }

                  LvArray::tensorOps::Ri_eq_symAijBj< 3 >( temp, particleStress[p], particleDamageGradient[p] );
                  real64 normalStress = LvArray::tensorOps::AiBi< 3 >( particleDamageGradient[p], temp );

                  real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
                  if( ( healType == 1 || healType == 3 || ( healType == 0 && ( normalStress < 0.0 || detF < 1.0 ) ) ) && particleDamage[p] > 0.0 )
                  {
                    particleCrystalHealFlag[p] = 1;

                    if( detF > 1.0 && healType == 3 )
                    {

                    // If there is a porosity model, healing can modify the material to be
                    // a porous material with a reference volume equal to the current volume,
                    // and a porosity that will allow compaction to the original material density
                    // Since we don't yet have porosity, this shouldn't be used.

                    // particleReferencePorosity[p] = 1.0 - 1.0 / detF;

                    // This def grad scaling needs to be modified for plane strain.  and it won't work
                    // well for anisotropic materials, so instead it will only be active for healType 3
                    // so we can still use 0,1 for the other cases.
                      real64 power = planeStrain ? 0.5 : 1.0 / 3.0;
                      real64 scaling = LvArray::math::pow( detF, power );


                      LvArray::tensorOps::scale< 3, 3 >( particleDeformationGradient[p], 1 / scaling );

                      LvArray::tensorOps::scale< 3 >( particleRVectors[p][0], scaling );
                      LvArray::tensorOps::scale< 3 >( particleRVectors[p][1], scaling );

                      if( !planeStrain )
                      {
                        LvArray::tensorOps::scale< 3 >( particleRVectors[p][2], scaling );
                      }

                      particleReferenceVolume[p] *= detF;
                    }

                  }
                  else
                  {
                    particleCrystalHealFlag[p] = 0;
                  }
                } );
                crystalHeal.setMarkedParticlesToHeal( 1 );
              }
              else
              {
                // Heal particles by gradually reducing their damage over a user-specified interval.
                // Explicitly force exactly zero damage at the end.
                if( time_n - dt/2 <= eventTime + eventInterval && eventTime + eventInterval < time_n + dt / 2 )
                {
                  forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
                  {
                    localIndex const p = activeParticleIndices[pp];
                    if( particleCrystalHealFlag[p] == 1 )
                    {
                      particleDamage[p] = 0.0;
                    }
                  } );
                  event.setIsComplete( 1 );
                }
                else
                {
                  forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
                  {
                    localIndex const p = activeParticleIndices[pp];
                    if( particleCrystalHealFlag[p] == 1 )
                    {
                      particleDamage[p] *= fmax( 0.0, 1.0 - dt * 20.0 * ( time_n - eventTime ) / ( eventInterval * eventInterval ) ); // Recursive form
                    }
                  } );
                }
              }

              if( event.isComplete() )
              {
                if( constitutiveModel.hasWrapper( "crackSpeed" ) )
                {
                  real64 & crackSpeed = constitutiveModel.getReference< real64 >( "crackSpeed" );
                  crackSpeed *= 1e100;
                }
              }
            }
          }
        } );
      }

      if( event.getName() == "MachineSample" )
      {
          MachineSampleMPMEvent & machineSample = dynamicCast< MachineSampleMPMEvent & >( event );

          real64 hEl[3] = {};
          LvArray::tensorOps::copy< 3 >(hEl, m_hEl);
          real64 xGlobalMin[3] = {};
          LvArray::tensorOps::copy< 3 >(xGlobalMin, m_xGlobalMin);
          real64 xGlobalMax[3] = {};
          LvArray::tensorOps::copy< 3 >(xGlobalMax, m_xGlobalMax);
          real64 domainExtent[3] = {};
          LvArray::tensorOps::copy< 3 >(domainExtent, m_domainExtent);

          particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
          {
            arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
            arrayView1d< int > const particleDeleteFlag = subRegion.getField< fields::mpm::particleDeleteFlag >();

            SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

            string sampleType = machineSample.getSampleType(); // cannot pass string into parallelDevicePolicy forAll use enum instead

            real64 gaugeRadius = machineSample.getGaugeRadius();
            real64 gaugeLength = machineSample.getGaugeLength();
            real64 filletRadius = machineSample.getFilletRadius();
            real64 machiningLength = 2*filletRadius + gaugeLength;

            real64 diskRadius = machineSample.getDiskRadius();

            real64 domainCenter[3] = {};
            LvArray::tensorOps::copy< 3 >( domainCenter, xGlobalMin);
            LvArray::tensorOps::add< 3 >( domainCenter, xGlobalMax);
            LvArray::tensorOps::scale< 3 >( domainCenter, 0.5 );

            forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
            {
              localIndex const p = activeParticleIndices[pp];

              // Dogbone is milled along y direction
              if( sampleType == "dogbone")
              {
                // Particle is inside gauge section
                if( particlePosition[p][1] > xGlobalMin[1] + (domainExtent[1] - machiningLength)/2 &&
                    particlePosition[p][1] < xGlobalMax[1] - (domainExtent[1] - machiningLength)/2 )
                {
                  real64 distSqr = LvArray::math::pow( particlePosition[p][0] - domainCenter[0], 2 ) + LvArray::math::pow( particlePosition[p][2] - domainCenter[2], 2 );
                  //Is particle inside the gauge section
                  if( particlePosition[p][1] > xGlobalMin[1] + (domainExtent[1] - gaugeLength)/2 &&
                      particlePosition[p][1] < xGlobalMax[1] - (domainExtent[1] - gaugeLength)/2 )
                  {
                    // Check distance in XZ plane from radius of gauge
                    if( distSqr  > gaugeRadius * gaugeRadius )
                    {
                      particleDeleteFlag[p] = 1;
                    }
                  }
                  //Check if particle is outside of fillet radius
                  else
                  {
                    real64 y = LvArray::math::abs( particlePosition[p][1] - domainCenter[1] ) - gaugeLength/2;
                    real64 rr = filletRadius - LvArray::math::sqrt( ( filletRadius * filletRadius - y * y ) ) + gaugeRadius;
                    if( distSqr > rr * rr)
                    {
                      particleDeleteFlag[p] = 1;
                    }
                  }
                }
              }

              if( sampleType == "brazilDisk")
              {
                real64 distSqr = LvArray::math::pow( particlePosition[p][0] - domainCenter[0], 2 ) + LvArray::math::pow( particlePosition[p][1] - domainCenter[1], 2 )  + LvArray::math::pow( particlePosition[p][2] - domainCenter[2], 2 );

                if( distSqr > diskRadius * diskRadius )
                {
                  particleDeleteFlag[p] = 1;
                }
              }

              if( sampleType == "cylinder" )
              {
                real64 distSqr = LvArray::math::pow( particlePosition[p][0] - domainCenter[0], 2 ) + LvArray::math::pow( particlePosition[p][2] - domainCenter[2], 2 );

                if( distSqr  > gaugeRadius * gaugeRadius )
                {
                  particleDeleteFlag[p] = 1;
                }

              }

            });
          });
          event.setIsComplete( 1 );
      }

      if( event.getName() == "FrictionCoefficientSwap" )
      {
        GEOS_LOG_RANK_0("Performing friction coefficient swap");
        FrictionCoefficientSwapMPMEvent & frictionCoefficientSwap = dynamicCast< FrictionCoefficientSwapMPMEvent & >( event );

        m_frictionCoefficient = frictionCoefficientSwap.getFrictionCoefficient();
        m_frictionCoefficientTable = frictionCoefficientSwap.getFrictionCoefficientTable();

        initializeFrictionCoefficients();

        event.setIsComplete( 1 );
      }

      if( event.getName() == "BodyForceUpdate" )
      {
        BodyForceUpdateMPMEvent & bodyForceUpdate = dynamicCast< BodyForceUpdateMPMEvent & >( event );

        LvArray::tensorOps::copy< 3 >( m_bodyForce, bodyForceUpdate.getBodyForce() );

        event.setIsComplete( 1 );
      }

      if( event.getName() == "DeformationUpdate" )
      {
        DeformationUpdateMPMEvent & deformationUpdate = dynamicCast< DeformationUpdateMPMEvent & >( event );

        int oldPrescribedFTable = m_prescribedFTable;
        int oldPrescribedBoundaryFTable = m_prescribedBoundaryFTable;
        int oldStressControl[3];
        GEOS_UNUSED_VAR( oldPrescribedFTable );
        GEOS_UNUSED_VAR( oldPrescribedBoundaryFTable );
        GEOS_UNUSED_VAR( oldStressControl );
        LvArray::tensorOps::copy< 3 >( oldStressControl, m_stressControl );

        // CC: TODO need to add special handling for turning off and on stress control because the F values will not be known
        m_prescribedFTable = deformationUpdate.getPrescribedFTable();
        m_prescribedBoundaryFTable = deformationUpdate.getPrescribedBoundaryFTable();
        m_stressControl = deformationUpdate.getStressControl();
      }

      if( event.getName() == "CohesiveZoneReference" )
      {
        m_referenceCohesiveZone = 1;
        m_enableCohesiveLaws = 1;

        event.setIsComplete( 1 );
      }
    }
  } );

}

void SolidMechanicsMPM::performMaterialSwap( ParticleManager & particleManager,
                                             string sourceRegionName,
                                             string destinationRegionName )
{
  // Material swap is performed by copying particle data from source region and constitutive model to destination region and constitutive model
  // The destination constiutive model is initialized as an empty particle subRegion

  // Should this really be a set? Only unique fields are allowable
  std::set< std::string > overwrittableFields = { fields::mpm::particleMass::key(),
                                                  fields::mpm::particleDensity::key() };

  ParticleRegion & sourceParticleRegion = particleManager.getRegion< ParticleRegion >( sourceRegionName );
  ParticleRegion & destinationParticleRegion = particleManager.getRegion< ParticleRegion > ( destinationRegionName );
  auto & sourceSubRegions = sourceParticleRegion.getSubRegions();
  auto & destinationSubRegions = destinationParticleRegion.getSubRegions();

  for( int r=0; r < sourceSubRegions.size(); ++r)
  {
    // Assumes indices of each subregions correspond to the same particle types (e.g. single point, CPTI, CPDI, etc.)
    ParticleSubRegion & sourceSubRegion = dynamicCast< ParticleSubRegion & >( *sourceSubRegions[r] );
    ParticleSubRegion & destinationSubRegion = dynamicCast< ParticleSubRegion & >( *destinationSubRegions[r] );

    // Resize destination subregion
    destinationSubRegion.resize( sourceSubRegion.size() );

    arrayView1d< real64 const > const sourceParticleVolume = sourceSubRegion.getParticleVolume();
    arrayView1d< real64 const > const sourceParticleReferenceVolume = sourceSubRegion.getField< fields::mpm::particleReferenceVolume >();

    // Get constitutive handles for density and state variables
    string const & sourceSolidMaterialName = sourceSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & sourceConstitutiveModel = getConstitutiveModel< ContinuumBase >( sourceSubRegion, sourceSolidMaterialName );
    arrayView3d< real64 const > const sourceOldStress = sourceConstitutiveModel.getReference< array3d< real64 > >( constitutive::ContinuumBase::viewKeyStruct::oldStressString() );

    string const & destinationSolidMaterialName = destinationSubRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & destinationConstitutiveModel = getConstitutiveModel< ContinuumBase >( destinationSubRegion, destinationSolidMaterialName );
    real64 const destinationDefaultConstitutiveDensity = destinationConstitutiveModel.getReference< real64 >( constitutive::ContinuumBase::viewKeyStruct::defaultDensityString() );
    arrayView3d< real64 > const destinationOldStress = destinationConstitutiveModel.getReference< array3d< real64 > >( constitutive::ContinuumBase::viewKeyStruct::oldStressString() );

    // Copy old stress from constitutive model separately for hypoelastic materials (is this stress really representative after a material swap?)
    forAll< serialPolicy >( sourceSubRegion.size(), [&]( localIndex const p )
    {
      LvArray::tensorOps::copy< 6 >( destinationOldStress[p][0], sourceOldStress[p][0] );
    } );

    sourceSubRegion.forWrappers( [&]( WrapperBase & sourceWrapper )
    {
      string const fieldName = sourceWrapper.getName();

      //Filter out only particle fields for copy by prefix
      if( fieldName.substr(0, 8) == "particle" )
      {
        GEOS_LOG_RANK( fieldName );

        WrapperBase & destinationWrapper = destinationSubRegion.getWrapperBase( fieldName );

        bool overwriteField = overwrittableFields.count( fieldName ) > 0; // Cannot use contains (C++20)


        types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
        {
          using ArrayType = camp::first< decltype( tupleOfTypes ) >;
          using T = typename ArrayType::ValueType;

          auto const sourceArray = Wrapper< ArrayType >::cast( sourceWrapper ).reference().toViewConst();
          auto destinationArray = Wrapper< ArrayType >::cast( destinationWrapper ).reference().toView();

          GEOS_ERROR_IF( sourceArray.size() != destinationArray.size(), "During material swap " << fieldName << "  fields did not have the same size!");

          forAll< serialPolicy >( sourceArray.size( 0 ), [&] GEOS_HOST ( localIndex const p )
          {
            // For scalar particle fields need to do assignment manually
            if constexpr( ArrayType::NDIM == 1 )
            {
              // If wrapper name is among those that should be overriden by new sub region skip (e.g. mass, density, etc.)
              // We currently only overwrite scalar quantities, but may need to adjust if we overwrite nonscalar quantities
              if ( overwriteField )
              {
                // Is there a better way to do this?
                if( fieldName == fields::mpm::particleMass::key() )
                {
                  destinationArray[p] = destinationDefaultConstitutiveDensity * sourceParticleReferenceVolume[p];
                }

                if( fieldName == fields::mpm::particleDensity::key() )
                {
                  destinationArray[p] = destinationDefaultConstitutiveDensity * sourceParticleReferenceVolume[p] / sourceParticleVolume[p];
                }
              }
              else
              {
                destinationArray[p] = sourceArray[p];
              }
            }
            else
            {
              auto sourceSlice = sourceArray[p];
              auto destinationSlice = destinationArray[p];
              LvArray::forValuesInSliceWithIndices( destinationSlice, [slice=sourceSlice] ( T & val, auto const ... indices )
              {
                val = slice( indices ... );
              } );
            }
          } );

        }, sourceWrapper );
      }
    } );

    // Remove particles from subregion since they now reside in destination subregion
    // Need to make a set to pass to erase ( maybe subregion needs a clear that deletes all particles or the like)
    // If we are resizing to zero do we need to erase indices?
    sourceSubRegion.resize( 0 );

    sourceSubRegion.setActiveParticleIndices();
    destinationSubRegion.setActiveParticleIndices();
  }
}

void SolidMechanicsMPM::syncGridFields( stdVector< std::string > const & fieldNames,
                                        DomainPartition & domain,
                                        NodeManager & nodeManager,
                                        MeshLevel & mesh,
                                        MPI_Op op )
{
  GEOS_MARK_FUNCTION;

  MPI_iCommData iComm;
  iComm.resize( domain.getNeighbors().size() );

  // (0) Bring grid fields to host
  for( auto const & name : fieldNames )
  {
    WrapperBase & wrapper = nodeManager.getWrapperBase( name );
    wrapper.move( LvArray::MemorySpace::host, true );
  }

  // (1) Initialize
  FieldIdentifiers fieldsToBeSynced;
  fieldsToBeSynced.addFields( FieldLocation::Node, fieldNames );
  stdVector< NeighborCommunicator > & neighbors = domain.getNeighbors();

  // (2) Swap send and receive indices so we can sum from ghost to master
  for( size_t n=0; n<neighbors.size(); ++n )
  {
    int const neighborRank = neighbors[n].neighborRank();
    array1d< localIndex > & nodeGhostsToReceive = nodeManager.getNeighborData( neighborRank ).ghostsToReceive();
    array1d< localIndex > & nodeGhostsToSend = nodeManager.getNeighborData( neighborRank ).ghostsToSend();
    array1d< localIndex > temp = nodeGhostsToSend;
    nodeGhostsToSend = nodeGhostsToReceive;
    nodeGhostsToReceive = temp;
  }

  // (3) Additive sync
  CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldsToBeSynced, mesh, neighbors, iComm, true );
  parallelDeviceEvents packEvents;
  CommunicationTools::getInstance().asyncPack( fieldsToBeSynced, mesh, neighbors, iComm, true, packEvents );
  waitAllDeviceEvents( packEvents );
  CommunicationTools::getInstance().asyncSendRecv( neighbors, iComm, true, packEvents );
  parallelDeviceEvents unpackEvents;
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, iComm, true, unpackEvents, op ); // needs an extra argument to
                                                                                                        // indicate unpack operation
  // (4) Swap send and receive indices back so we can sync from master to ghost
  for( size_t n=0; n<neighbors.size(); ++n )
  {
    int const neighborRank = neighbors[n].neighborRank();
    array1d< localIndex > & nodeGhostsToReceive = nodeManager.getNeighborData( neighborRank ).ghostsToReceive();
    array1d< localIndex > & nodeGhostsToSend = nodeManager.getNeighborData( neighborRank ).ghostsToSend();
    array1d< localIndex > temp = nodeGhostsToSend;
    nodeGhostsToSend = nodeGhostsToReceive;
    nodeGhostsToReceive = temp;
  }

  // (5) Perform sync
  CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldsToBeSynced, mesh, neighbors, iComm, true );
  parallelDeviceEvents packEvents2;
  CommunicationTools::getInstance().asyncPack( fieldsToBeSynced, mesh, neighbors, iComm, true, packEvents2 );
  waitAllDeviceEvents( packEvents2 );
  CommunicationTools::getInstance().asyncSendRecv( neighbors, iComm, true, packEvents2 );
  parallelDeviceEvents unpackEvents2;
  CommunicationTools::getInstance().finalizeUnpack( mesh, neighbors, iComm, true, unpackEvents2 );
}

void SolidMechanicsMPM::singleFaceVectorFieldSymmetryBC( const int face,
                                                         arrayView3d< real64 > const & vectorMultiField,
                                                         arrayView3d< real64 > const & dVectorMultiField, //Change in vectorMultiField
                                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                         Group & nodeSets )
{
  // This is a helper function for enforcing symmetry BCs on a single face and is meant to be called by other functions, not directly by the
  // solver:
  //   * enforceGridVectorFieldSymmetryBC calls this on all grid faces
  //   * applyEssentialBCs calls this on faces that aren't moving (moving faces due to F-table need special treatment)
  int nEl[3] = {};
  LvArray::tensorOps::copy< 3 >( nEl, m_nEl );
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  array1d< SortedArray< localIndex > > & boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );

  for( localIndex fieldIndex=0; fieldIndex < m_numVelocityFields; ++fieldIndex )
  {
    // Face-associated quantities
    int dir0 = face / 2;           // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)
    int dir1 = ( dir0 + 1 ) % 3;   // 1, 1, 2, 2, 0, 0
    int dir2 = ( dir0 + 2 ) % 3;   // 2, 2, 0, 0, 1, 1
    int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1

    // Enforce BCs on boundary nodes
    SortedArrayView< localIndex const > const boundaryNodesView = boundaryNodes[face].toView();
    int const numBoundaryNodes = boundaryNodesView.size();
    // Probably not a big enough loop to warrant parallelization
    forAll< serialPolicy >( numBoundaryNodes, [=] GEOS_HOST ( localIndex const gg )
    {
      int const g = boundaryNodesView[gg];
      dVectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[g][fieldIndex][dir0];
      vectorMultiField[g][fieldIndex][dir0] = 0.0;
    } );

    // Perform field reflection on buffer nodes
    SortedArrayView< localIndex const > const bufferNodesView = bufferNodes[face].toView();
    int const numBufferNodes = bufferNodesView.size();
    // Probably not a big enough loop to warrant parallelization
    forAll< serialPolicy >( numBufferNodes, [=] GEOS_HOST ( localIndex const gg )
    {
      int const g = bufferNodesView[gg];
      int ijk[3];
      ijk[dir0] = positiveNormal * ( nEl[dir0] - 2 ) + ( 1 - positiveNormal ) * ( 2 );
      ijk[dir1] = LvArray::math::round( ( gridPosition[g][dir1] - xLocalMin[dir1] ) / hEl[dir1] );
      ijk[dir2] = LvArray::math::round( ( gridPosition[g][dir2] - xLocalMin[dir2] ) / hEl[dir2] );

      localIndex gFrom = ijkMap[ijk[0]][ijk[1]][ijk[2]];

      real64 previousVectorMultiField[3] = {};
      previousVectorMultiField[dir0] = vectorMultiField[g][fieldIndex][dir0];
      previousVectorMultiField[dir1] = vectorMultiField[g][fieldIndex][dir1];
      previousVectorMultiField[dir2] = vectorMultiField[g][fieldIndex][dir2];

      vectorMultiField[g][fieldIndex][dir0] = -vectorMultiField[gFrom][fieldIndex][dir0]; // Negate component aligned with surface normal
      vectorMultiField[g][fieldIndex][dir1] =  vectorMultiField[gFrom][fieldIndex][dir1];
      vectorMultiField[g][fieldIndex][dir2] =  vectorMultiField[gFrom][fieldIndex][dir2];

      dVectorMultiField[g][fieldIndex][dir0] += vectorMultiField[g][fieldIndex][dir0] - previousVectorMultiField[dir0];
      dVectorMultiField[g][fieldIndex][dir1] += vectorMultiField[g][fieldIndex][dir1] - previousVectorMultiField[dir1];
      dVectorMultiField[g][fieldIndex][dir2] += vectorMultiField[g][fieldIndex][dir2] - previousVectorMultiField[dir2];
    } );
  }
}

void SolidMechanicsMPM::enforceGridVectorFieldSymmetryBC( arrayView3d< real64 > const & vectorMultiField,
                                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                          Group & nodeSets )
{
  for( int face=0; face < 6; ++face )
  {
    if( m_boundaryConditionTypes[face] == 1 || m_boundaryConditionTypes[face] == 2 ) // || m_boundaryConditionTypes[face] == 3 )// What about for BC = 3 (contact)?
    {
      array3d< real64 > dVectorMultiField( vectorMultiField.size(0), vectorMultiField.size(1), vectorMultiField.size(2) ); // dumby variable, CC TODO there must be a better solution
      singleFaceVectorFieldSymmetryBC( face, vectorMultiField, dVectorMultiField, gridPosition, nodeSets );
    }
  }
}

void SolidMechanicsMPM::applyEssentialBCs( const real64 dt,
                                           const real64 time_n,
                                           NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // Get grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView1d< int const > const gridGhostRank = nodeManager.ghostRank();
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 > const gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVelocityString() );
  arrayView3d< real64 > const gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );

  arrayView2d< real64 const > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
  arrayView3d< real64 const > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );
  arrayView3d< real64 const > const gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  // Get node sets
  Group & nodeSets = nodeManager.sets();
  array1d< SortedArray< localIndex > > & boundaryNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::boundaryNodesString() );
  array1d< SortedArray< localIndex > > & bufferNodes = nodeSets.getReference< array1d< SortedArray< localIndex > > >( viewKeyStruct::bufferNodesString() );

  // Impose BCs on each face while gathering reaction forces
  real64 localFaceReactions[6];
  for( int face = 0; face < 6; ++face )
  {
    // TODO Eventually perform cast to BC enum type!
    switch( m_boundaryConditionTypes[face] )
    {
      case 0: // Outflow
        break; // Do nothing
      case 1: // Symmetry
        {
          // CC: TODO add gridDVelocity update to theses for XPIC
          singleFaceVectorFieldSymmetryBC( face, gridVelocity, gridDVelocity, gridPosition, nodeSets );

          //Dumby paramter, we do not need the change in grid acceleration
          array3d< real64 > gridDAcceleration( gridAcceleration.size(0), gridAcceleration.size(1), gridAcceleration.size(2)); //CC: TODO Probably want to avoid allocating a lot of memory just for the dummy variable
          singleFaceVectorFieldSymmetryBC( face, gridAcceleration, gridDAcceleration, gridPosition, nodeSets );
        }
        break;
      case 2: // Moving
      case 3: // Contact
        if ( m_prescribedBoundaryFTable == 1 || m_stressControl[0] == 1 || m_stressControl[1] == 1 || m_stressControl[2] == 1 || m_boundaryConditionTypes[face] == 3 ) // Double check stress control (think we only need to check if stress control for the direction of faces are on not all them)
        {
          for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex )
          {
            // Face-associated quantities
            int dir0 = face / 2;           // 0, 0, 1, 1, 2, 2 (x-, x+, y-, y+, z-, z+)
            int dir1 = (dir0 + 1) % 3;     // 1, 1, 2, 2, 0, 0
            int dir2 = (dir0 + 2) % 3;     // 2, 2, 0, 0, 1, 1
            int positiveNormal = face % 2; // even => (-) => 0, odd => (+) => 1

            // Enforce BCs on boundary nodes using F-table
            SortedArrayView< localIndex const > const boundaryNodesView = boundaryNodes[face].toView();
            int const numBoundaryNodes = boundaryNodesView.size();
            // Probably not a big enough loop to warrant parallelization
            forAll< serialPolicy >( numBoundaryNodes, [&, gridPosition, gridVelocity, gridDVelocity, gridMass] GEOS_HOST ( localIndex const gg )
            {
              int const g = boundaryNodesView[gg];

              // If boundary condition type is 3 (e.g. contact, default is sticky), performs a check if normal grid component of grid velocity is moving into plane, if so flips component and adds domain velocity
              real64 prescribedVelocity = gridVelocity[g][fieldIndex][dir0];
              if( m_boundaryConditionTypes[face] == 3 )
              {
                // Currently fails if contact is only enforced when surface position is at 0 for inContact condition on grid boundary too, likely needs some soft scaling
                real64 surfacePosition = gridCenterOfVolume[g][fieldIndex][dir0];
                if( gridSurfaceFieldMass[g][fieldIndex] > m_smallMass )
                {
                  surfacePosition = gridSurfacePosition[g][fieldIndex][dir0];
                }
                bool inContact = ( gridVelocity[g][fieldIndex][dir0] * ( -1.0 + 2.0 * positiveNormal ) > 0.0 ) && ( surfacePosition * ( -1.0 + 2.0 * positiveNormal ) > 0.0 );

                if( inContact )
                {
                  prescribedVelocity = inContact * ( m_domainL[dir0] * gridPosition[g][dir0] );
                // prescribedVelocity = inContact * ( -m_boundaryFaceCoefficientsOfRestitution[face]*gridVelocity[g][fieldIndex][dir0] + m_domainL[dir0] * gridPosition[g][dir0] );
                }

                // Enforce friction in transverse directions
                real64 mu = m_boundaryFaceFrictionCoefficients[face];
                real64 frictionalForce = mu * fmax(-gridAcceleration[g][face][dir0] * ( -1.0 + 2.0 * positiveNormal ), 0.0);

                real64 inPlaneSpeed = LvArray::math::sqrt( gridVelocity[g][fieldIndex][dir1] * gridVelocity[g][fieldIndex][dir1] + gridVelocity[g][fieldIndex][dir2] * gridVelocity[g][fieldIndex][dir2] );
                real64 r1 = gridVelocity[g][fieldIndex][dir1] / inPlaneSpeed;
                real64 r2 = gridVelocity[g][fieldIndex][dir2] / inPlaneSpeed;

                real64 da1 = r1 * frictionalForce;
                real64 da2 = r2 * frictionalForce;

                gridDVelocity[g][fieldIndex][dir1] = da1 * dt;
                gridDVelocity[g][fieldIndex][dir2] = da2 * dt;

                gridVelocity[g][fieldIndex][dir1] += gridDVelocity[g][fieldIndex][dir1];
                gridVelocity[g][fieldIndex][dir1] += gridDVelocity[g][fieldIndex][dir2];

                gridAcceleration[g][fieldIndex][dir1] -= da1;
                gridAcceleration[g][fieldIndex][dir2] -= da2;
              }
              else
              {
                prescribedVelocity = m_domainL[dir0] * gridPosition[g][dir0];

                if(m_enablePrescribedBoundaryTransverseVelocities[face] == 1)
                {
                  real64 prescribedTransverseVelocity1 = m_prescribedBoundaryTransverseVelocities[face][0];
                  gridDVelocity[g][fieldIndex][dir1] = prescribedTransverseVelocity1 - gridVelocity[g][fieldIndex][dir1];
                  real64 accelerationForTransverseBC1 = gridDVelocity[g][fieldIndex][dir1] / dt; // acceleration needed to satisfy BC along transverse directions
                  gridVelocity[g][fieldIndex][dir1] = prescribedTransverseVelocity1;
                  gridAcceleration[g][fieldIndex][dir1] += accelerationForTransverseBC1;

                  real64 prescribedTransverseVelocity2 = m_prescribedBoundaryTransverseVelocities[face][1];
                  gridDVelocity[g][fieldIndex][dir2] = prescribedTransverseVelocity2 - gridVelocity[g][fieldIndex][dir2];
                  real64 accelerationForTransverseBC2 = gridDVelocity[g][fieldIndex][dir2] / dt; // acceleration needed to satisfy BC along transverse directions
                  gridVelocity[g][fieldIndex][dir2] = prescribedTransverseVelocity2;
                  gridAcceleration[g][fieldIndex][dir2] += accelerationForTransverseBC2;
                }
              }

              gridDVelocity[g][fieldIndex][dir0] = prescribedVelocity - gridVelocity[g][fieldIndex][dir0]; // CC: TODO double check this, because it overrides the change in velocity that might have been written during enforceContact
              real64 accelerationForBC = gridDVelocity[g][fieldIndex][dir0] / dt; // acceleration needed to satisfy BC
              gridVelocity[g][fieldIndex][dir0] = prescribedVelocity;
              gridAcceleration[g][fieldIndex][dir0] += accelerationForBC;

              if( gridGhostRank[g] <= -1 ) // so we don't double count reactions at partition boundaries
              {
                localFaceReactions[face] += accelerationForBC * gridMass[g][fieldIndex];
              }
            } );

            // Perform field reflection on buffer nodes - accounts for moving boundary effects
            SortedArrayView< localIndex const > const bufferNodesView = bufferNodes[face].toView();
            int const numBufferNodes = bufferNodesView.size();
            // Possibly not a big enough loop to warrant parallelization
            forAll< serialPolicy >( numBufferNodes, [&, gridPosition, gridVelocity, gridDVelocity, gridAcceleration] GEOS_HOST ( localIndex const gg )
              {
                int const g = bufferNodesView[gg];

                // Initialize grid ijk indices
                int ijk[3];
                ijk[dir1] = LvArray::math::round((gridPosition[g][dir1] - m_xLocalMin[dir1]) / m_hEl[dir1] );
                ijk[dir2] = LvArray::math::round((gridPosition[g][dir2] - m_xLocalMin[dir2]) / m_hEl[dir2] );

                // Grab the node index that we're copying from
                ijk[dir0] = positiveNormal * (m_nEl[dir0] - 2) + (1 - positiveNormal) * (2);
                localIndex gFrom = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

                // Grab the associated boundary node index for moving boundary correction
                ijk[dir0] = positiveNormal * (m_nEl[dir0] - 1) + (1 - positiveNormal) * (1);
                localIndex gBoundary = m_ijkMap[ijk[0]][ijk[1]][ijk[2]];

                //Store previous velocity to compute change in velocity
                real64 gridPreviousVelocity[3] = {};
                gridPreviousVelocity[dir0] = gridVelocity[g][fieldIndex][dir0];
                gridPreviousVelocity[dir1] = gridVelocity[g][fieldIndex][dir1];
                gridPreviousVelocity[dir2] = gridVelocity[g][fieldIndex][dir2];

                // Calculate velocity, Negate component aligned with surface normal and correct for moving boundary
                gridVelocity[g][fieldIndex][dir0] = -gridVelocity[gFrom][fieldIndex][dir0] + 2.0 * gridVelocity[gBoundary][fieldIndex][dir0];
                gridVelocity[g][fieldIndex][dir1] = gridVelocity[gFrom][fieldIndex][dir1];
                gridVelocity[g][fieldIndex][dir2] = gridVelocity[gFrom][fieldIndex][dir2];

                // Compute change in velocity for XPIC calculations
                gridDVelocity[g][fieldIndex][dir0] += gridVelocity[g][fieldIndex][dir0] - gridPreviousVelocity[dir0];
                gridDVelocity[g][fieldIndex][dir1] += gridVelocity[g][fieldIndex][dir1] - gridPreviousVelocity[dir1];
                gridDVelocity[g][fieldIndex][dir2] += gridVelocity[g][fieldIndex][dir2] - gridPreviousVelocity[dir2];

                // Calculate acceleration, Negate component aligned with surface normal and correct for moving boundary
                gridAcceleration[g][fieldIndex][dir0] = -gridAcceleration[gFrom][fieldIndex][dir0] + 2.0 * gridAcceleration[gBoundary][fieldIndex][dir0];
                gridAcceleration[g][fieldIndex][dir1] = gridAcceleration[gFrom][fieldIndex][dir1];
                gridAcceleration[g][fieldIndex][dir2] = gridAcceleration[gFrom][fieldIndex][dir2];
              } );
          }
        }
        break;
      default:
        GEOS_ERROR("Unrecognized boundary condition type in MPM Solver!");
        break;
    }
  }

  // Compute change in grid velocities for XPIC calculations
  // TODO for multifield contact corrections

  // Reduce reaction forces from all partitions
  real64 globalFaceReactions[6];
  for( int face = 0; face < 6; face++ )
  {
    MPI_Allreduce( &localFaceReactions[face],
                   &globalFaceReactions[face],
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   MPI_COMM_GEOS );
  }
  LvArray::tensorOps::copy< 6 >( m_globalFaceReactions, globalFaceReactions );

  // Get end-of-step domain dimensions - note that m_domainExtent is updated later
  real64 length, width, height;
  length = m_domainExtent[0] * (1.0 + m_domainL[0] * dt);
  width  = m_domainExtent[1] * (1.0 + m_domainL[1] * dt);
  height = m_domainExtent[2] * (1.0 + m_domainL[2] * dt);

  // Write global reactions to file
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 && m_reactionHistory == 1 )
  {
    if(time_n + dt >= m_nextReactionWriteTime )
    {
      std::ofstream file;
      // can't enable exception now because of gcc bug that raises ios_base::failure with useless message
      // file.exceptions(file.exceptions() | std::ios::failbit);
      file.open( "reactionHistory.csv", std::ios::out | std::ios::app );
      if( file.fail() )
      {
        throw std::ios_base::failure( std::strerror( errno ) );
      }
      // make sure write fails with exception if something is wrong
      file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
      file << std::setprecision( std::numeric_limits< long double >::digits10 )
          << time_n + dt << ","
          << m_domainF[0] << ","
          << m_domainF[1] << ","
          << m_domainF[2] << ","
          << length << ","
          << width << ","
          << height << ","
          << globalFaceReactions[0] << ","
          << globalFaceReactions[1] << ","
          << globalFaceReactions[2] << ","
          << globalFaceReactions[3] << ","
          << globalFaceReactions[4] << ","
          << globalFaceReactions[5] << ","
          << m_domainL[0] << ","
          << m_domainL[1] << ","
          << m_domainL[2] << std::endl;
      file.close();
      m_nextReactionWriteTime += m_reactionWriteInterval;
    }
  }
}

// Should this be removed from contact and included with checks inside particleToGrid?
void SolidMechanicsMPM::computeGridSurfaceNormals( ParticleManager & particleManager,
                                                   NodeManager & nodeManager )
{

  //CC: fields for on-the-fly mapNodes calculations
  int const numDims = m_numDims;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numContactGroups = m_numContactGroups;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  // arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  // arrayView2d< real64 > const gridSurfaceNormalWeights = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightsString() );
  // arrayView2d< real64 > const gridSurfaceNormalWeightNormalization = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightNormalizationString() );
  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    // Get views to mapping arrays
    int const explicitSurfaceNormalInfluence = m_explicitSurfaceNormalInfluence;
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    ParticleType const particleType = subRegion.getParticleType();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp ) // Can parallelize with atomics
    {
      localIndex const p = activeParticleIndices[pp];

      // On-the-fly shape function computation
      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues);

      // Map to grid
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        localIndex const mappedNode = mappedNodes[g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        for( int i = 0; i < numDims; ++i )
        {
          real64 particleContributionToGrid = 0.0;
          particleContributionToGrid += shapeFunctionGradientValues[g][i] * particleVolume[p];
          if( particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ) // Update this with enum type for type safety to specifically only implement 2 and 3 for now (those with explicit surface normals)
          {
            // Also maps explicit particle surface normals which will dominate if m_explicitSurfaceNormalInfluence is large
            // If particle surface normal was disabled due to damage or CPDI domain scaling (e.g. zeroed) then the following does not add anything to gridSurfaceNormal
            particleContributionToGrid += explicitSurfaceNormalInfluence * particleSurfaceNormal[p][i] * shapeFunctionValues[g] * particleVolume[p];
          }
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceNormal[mappedNode][fieldIndex][i], particleContributionToGrid );
        }
      }
    } ); // particle loop
  } ); // subregion loop

  gridSurfaceNormal.move( LvArray::MemorySpace::host );
}

void SolidMechanicsMPM::computeGridSurfacePositions( ParticleManager & particleManager,
                                                     NodeManager & nodeManager )
{
  int const numDims = m_numDims;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numContactGroups = m_numContactGroups;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin);
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  // arrayView3d< real64 const > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView2d< real64 > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
  arrayView3d< real64 > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    // arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    ParticleType const particleType = subRegion.getParticleType();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      real64 particleContributionToGrid = 0.0;

      // This may need to change if we add more surface flag options assumes all surface flags with value greater than one have explicit surface positions and normals defined
      // if( particleSurfaceFlag[p] > 1 )
      if( LvArray::tensorOps::l2Norm< 3 >( particleSurfaceNormal[p] ) > 1e-12 )
      {
        // On-the-fly shape function computation
        localIndex mappedNodes[64];
        real64 shapeFunctionValues[64];
        real64 shapeFunctionGradientValues[64][3];
        mapNodesAndComputeShapeFunctions( ijkMap,
                                          xLocalMin,
                                          hEl,
                                          particleType,
                                          particlePosition[p],
                                          particleRVectors[p],
                                          gridPosition,
                                          mappedNodes,
                                          shapeFunctionValues,
                                          shapeFunctionGradientValues );

        // Map to grid
        for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
        {
          localIndex const mappedNode = mappedNodes[g];
          real64 const shapeFunctionValue = shapeFunctionValues[g];

          localIndex const fieldIndex = partitionField( numContactGroups,
                                                 damageFieldPartitioning,
                                                 particleGroup[p],
                                                 particleDamageGradient[p],
                                                 particleSurfaceNormal[p],
                                                 gridDamageGradient[mappedNode] );

          particleContributionToGrid = shapeFunctionValue * particleMass[p];
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceFieldMass[mappedNode][fieldIndex], particleContributionToGrid );

          real64 surfacePositionRelativeToNode[3];
          LvArray::tensorOps::copy< 3 >( surfacePositionRelativeToNode, particleSurfacePosition[p] );
          LvArray::tensorOps::add< 3 >( surfacePositionRelativeToNode, particlePosition[p] );
          LvArray::tensorOps::subtract< 3 >( surfacePositionRelativeToNode, gridPosition[mappedNode] );

          // TODO: Check which of these produces better results
          // We need to use the surface position relative to the grid node
          // Two possible options map relative surface position either using particle surface normal or grid surface normal
          real64 surfacePositionAlongNormal = LvArray::tensorOps::AiBi< 3 >( surfacePositionRelativeToNode, particleSurfaceNormal[p] );

          for( int i=0; i < numDims; ++i )
          {
            particleContributionToGrid = shapeFunctionValue * particleMass[p] * surfacePositionAlongNormal * particleSurfaceNormal[p][i];
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfacePosition[mappedNode][fieldIndex][i], particleContributionToGrid );
          }
        }
      }

    } ); // particle loop
  } ); // subregion loop

  gridSurfacePosition.move( LvArray::MemorySpace::host );
}

void SolidMechanicsMPM::normalizeGridSurfaceNormalsAndPositions( NodeManager & nodeManager )
{
  localIndex const numVelocityFields = m_numVelocityFields;
  int const planeStrain = m_planeStrain;
  real64 const smallMass = m_smallMass;

  arrayView2d< real64 const > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView2d< real64 const > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
  arrayView3d< real64 > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );

  forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      arraySlice1d< real64 > const surfaceNormal = gridSurfaceNormal[g][fieldIndex];
      if( gridMass[g][fieldIndex] > smallMass )
      {
        real64 norm = planeStrain == 1 ? LvArray::math::sqrt( surfaceNormal[0] * surfaceNormal[0] + surfaceNormal[1] * surfaceNormal[1] ) : LvArray::tensorOps::l2Norm< 3 >( surfaceNormal );
        if( norm > 1e-12 ) // TODO: Pick a good global threhsold, probably should be a user settable input
        {
          LvArray::tensorOps::scale< 3 >( surfaceNormal, 1.0 / norm );
        }

        if( gridSurfaceFieldMass[g][fieldIndex] > 1e-12)
        {
          LvArray::tensorOps::scale< 3 >( gridSurfacePosition[g][fieldIndex], 1.0 / gridSurfaceFieldMass[g][fieldIndex] );
        }
        else
        {
          LvArray::tensorOps::fill< 3 >( gridSurfacePosition[g][fieldIndex], 0.0 );
        }
        continue;
      }

      LvArray::tensorOps::fill< 3 >( surfaceNormal, 0.0);
    }
  } );

  // gridSurfaceNormal.move( LvArray::MemorySpace::host );
  // gridSurfacePosition.move( LvArray::MemorySpace::host );
}

// void SolidMechanicsMPM::normalizeGridSurfacePositions( NodeManager & nodeManager )
// {
//   real64 const smallMass = m_smallMass;
//   int const numVelocityFields = m_numVelocityFields;

//   arrayView2d< real64 const > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
//   arrayView3d< real64 > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );

//   forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
//   {
//     for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
//     {
//       if( gridSurfaceFieldMass[g][fieldIndex] > smallMass )
//       {
//         LvArray::tensorOps::scale< 3 >( gridSurfacePosition[g][fieldIndex], 1.0 / gridSurfaceFieldMass[g][fieldIndex] );
//         continue;
//       }

//       LvArray::tensorOps::fill< 3 >( gridSurfacePosition[g][fieldIndex], 0.0 );
//     }
//   } );

//   gridSurfacePosition.move( LvArray::MemorySpace::host );
// }

void SolidMechanicsMPM::computeNodalAreas( NodeManager & nodeManager )
{
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  localIndex const numVelocityFields = m_numVelocityFields;
  localIndex const numSurfaceIntegrationPoints = m_numSurfaceIntegrationPoints;

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 const > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView3d< real64 const > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );
  arrayView2d< real64 > const gridSurfaceArea = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceAreaString() );

  real64 const L = LvArray::tensorOps::l2Norm< 3 >( hEl );
  real64 const dA = pow( 2 * L / numSurfaceIntegrationPoints, 2 ); 

  GEOS_LOG_RANK("hEl: {" << hEl[0] << ", " << hEl[1] << ", " << hEl[2] <<  "}, L: " << L << ", dA: " << dA);

 real64 total_area = 0.0;
 for( localIndex g = 0; g < nodeManager.size(); ++g )
  {
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      gridSurfaceArea[g][fieldIndex] = 0.0; //Initialize this to zero

      // if( LvArray::tensorOps::l2Norm< 3 >( gridSurfaceNormal[g][fieldIndex] ) < 1e-16)
      // {
      //   continue;
      // }

      // Construct basis for integration of surface plane
      real64 n[3];
      if( m_useNodePosForArea == 1)
      {
        LvArray::tensorOps::copy< 3 >( n, gridPosition[g] );
      }
      else
      {
        LvArray::tensorOps::copy< 3 >( n, gridSurfaceNormal[g][fieldIndex] );
      }

      real64 r = LvArray::tensorOps::l2Norm< 3 >( n );
      if( abs(r) > 1e-16 )
      {
        LvArray::tensorOps::scale< 3 >( n, 1/r );
      }
      else
      {
        continue;
      }

      real64 distanceToSurface = 0.0;
      if(m_useNodePosForArea == 1)
      {
        distanceToSurface = 0.4 - r;
      }
      else
      {
        distanceToSurface = LvArray::tensorOps::AiBi< 3 >( gridSurfacePosition[g][fieldIndex], n );
      }

      real64 s1[3];
      real64 s2[3];
      computeOrthonormalBasis( n,
                               s1,
                               s2 );

      // Numerically integrate surface plane using trilinear grid shape functions
      for( localIndex i=0; i < numSurfaceIntegrationPoints; ++i )
      {
        real64 eta = 2.0 * i / ( numSurfaceIntegrationPoints - 1 ) - 1.0; // Natural coordinate of surface
        for( localIndex j=0; j < numSurfaceIntegrationPoints; ++j )
        {
          real64 xi = 2.0 * j / ( numSurfaceIntegrationPoints - 1 ) - 1.0; // Natural coordinate of surface

          real64 surfacePoint[3];
          surfacePoint[0] = distanceToSurface * n[0] + L * eta * s1[0] + L * xi * s2[0];
          surfacePoint[1] = distanceToSurface * n[1] + L * eta * s1[1] + L * xi * s2[1];
          surfacePoint[2] = distanceToSurface * n[2] + L * eta * s1[2] + L * xi * s2[2];

          if( ( -hEl[0] < surfacePoint[0] ) &&
              ( surfacePoint[0] < hEl[0] ) &&
              ( -hEl[1] < surfacePoint[1] ) &&
              ( surfacePoint[1] < hEl[1] ) &&
              ( -hEl[2] < surfacePoint[2] ) &&
              ( surfacePoint[2] < hEl[2] ) )
          {
            gridSurfaceArea[g][fieldIndex] += ( 1 - fabs( surfacePoint[0] / hEl[0] ) ) * ( 1 - fabs( surfacePoint[1] / hEl[1]) ) * ( 1 - fabs( surfacePoint[2] / hEl[2]) );
          }
        }
      }
      // nodalArea[g][fieldIndex] *= dA * nodalTotalVolume[g] / ( hEl[0] * hEl[1] * hEl[2] );
      GEOS_UNUSED_VAR( gridMaterialVolume );
      gridSurfaceArea[g][fieldIndex] *= dA;
      total_area += fieldIndex == 0 ? gridSurfaceArea[g][fieldIndex] : 0.0;
      GEOS_LOG_RANK_IF( fieldIndex == 0, g << ", " <<
                     fieldIndex << ", " <<
                     gridPosition[g][0] << ", " << gridPosition[g][1] << ", " << gridPosition[g][2] << ", " <<
                     gridSurfaceNormal[g][fieldIndex][0] << ", " << gridSurfaceNormal[g][fieldIndex][1] << ", " << gridSurfaceNormal[g][fieldIndex][2] << ", " << 
                     gridSurfacePosition[g][fieldIndex][0] << ", " << gridSurfacePosition[g][fieldIndex][1] << ", " << gridSurfacePosition[g][fieldIndex][2] << ", " <<
                     n[0] << ", " << n[1] << ", " << n[2] << ", " << 
                     distanceToSurface << ", " << 
                     gridSurfaceArea[g][fieldIndex]);
    }
  }
  GEOS_LOG_RANK("Total surface area for field 0 = " << total_area);
}

void SolidMechanicsMPM::computeGridSurfaceNormalWeights( ParticleManager & particleManager,
                                                         NodeManager & nodeManager )
{
  // Local copies for on-the-fly mapNodes calculations
  // int const numDims = m_numDims;
  real64 const smallMass = m_smallMass;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numVelocityFields = m_numVelocityFields;
  int const numContactGroups = m_numContactGroups;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  // arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 > const gridSurfaceNormalWeights = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightsString() );
  arrayView2d< real64 > const gridSurfaceNormalWeightNormalization = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightNormalizationString() );
  arrayView3d< real64 const > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  // arrayView1d< real64 const > const gridSurfaceMass = nodeManager.getReference< array1d< real64 > >( viewKeyStruct::gridSurfaceMassString() );
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // arrayView1d< int const > particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    // arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp ) // Can parallelize with atomics
    {
      localIndex const p = activeParticleIndices[pp];

      real64 particleContributionToGrid;

      // Map to grid
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        localIndex const mappedNode = mappedNodes[pp][g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        real64 surfaceNormal[3] = {};
        LvArray::tensorOps::copy< 3 >( surfaceNormal, particleSurfaceNormal[p] );

        particleContributionToGrid = LvArray::tensorOps::AiBi< 3 >( gridSurfaceNormal[mappedNode][fieldIndex], surfaceNormal ) * shapeFunctionValues[pp][g] * particleMass[p];
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceNormalWeights[mappedNode][fieldIndex], particleContributionToGrid );

        if( LvArray::tensorOps::l2NormSquared< 3 >( particleSurfaceNormal[p] ) > 1e-16 )
        {
          particleContributionToGrid = shapeFunctionValues[pp][g] * particleMass[p];
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceNormalWeightNormalization[mappedNode][fieldIndex], particleContributionToGrid );
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop

  // Loop over grid nodes and their fields to normalize the grid surface normal weights
  forAll< serialPolicy >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      if( gridSurfaceNormalWeightNormalization[g][fieldIndex] > smallMass ) // small mass threshold
      {
        gridSurfaceNormalWeights[g][fieldIndex] /= gridSurfaceNormalWeightNormalization[g][fieldIndex];
        continue;
      }

      gridSurfaceNormalWeights[g][fieldIndex] = 0.0;
    }
  } );
}

void SolidMechanicsMPM::initializeFrictionCoefficients()
{
  if( m_frictionCoefficientTable.size(0) != 0 )
  {
    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) != m_frictionCoefficientTable.size(1), "frictionCoefficientTable must be square.");
    GEOS_ERROR_IF( m_frictionCoefficientTable.size(0) != m_numContactGroups, "frictionCoefficientTable must have the same number of rows and columns as the number of contact groups.");

    for(int i = 0; i < m_numContactGroups; ++i)
    {
      for(int j = i+1; i < m_numContactGroups; ++j)
      {
        GEOS_ERROR_IF( LvArray::math::abs(m_frictionCoefficientTable[i][j] - m_frictionCoefficientTable[j][i]) > DBL_EPSILON, "Off-diagonal friction coefficients must match" );
      }
    }
    return;
  }
  else
  {
    if( isEqual( m_frictionCoefficient, -1.0) )
    // if( static_cast<int>( m_frictionCoefficient ) == -1 )
    {
      m_frictionCoefficient = 0.0;
    }

    GEOS_ERROR_IF( m_frictionCoefficient < 0.0, "Friction coefficient must be positive.");

    m_frictionCoefficientTable.resize( m_numContactGroups, m_numContactGroups );
    for(int i = 0; i < m_numContactGroups; ++i)
    {
      for(int j = 0; j < m_numContactGroups; ++j)
      {
        m_frictionCoefficientTable[i][j] = m_frictionCoefficient;
      }
    }
  }
}

void SolidMechanicsMPM::computeContactForces( real64 const dt,
                                              NodeManager & nodeManager )
{
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );

  ContactNormalTypeOption const contactNormalType = m_contactNormalType;
  ContactGapCorrectionOption const contactGapCorrection = m_contactGapCorrection;
  OverlapCorrectionOption const overlapCorrection = m_overlapCorrection;

  int const planeStrain = m_planeStrain;
  int const preventCZInterpentration = m_preventCZInterpentration;
  int const treatFullyDamagedAsSingleField = m_treatFullyDamagedAsSingleField;
  int const useSurfacePositionForContact = m_useSurfacePositionForContact;

  int const numContactGroups = m_numContactGroups;
  int const numVelocityFields = m_numVelocityFields;
 
  real64 const smallMass = m_smallMass;
  real64 const neighborRadius = m_neighborRadius;
  real64 const separabilityMinDamage = m_separabilityMinDamage;
  real64 const thinFeatureDFGThreshold = m_thinFeatureDFGThreshold;

  array2d< real64 > frictionCoefficientTableCopy( numContactGroups, numContactGroups );
  for( int i = 0; i < numContactGroups; ++i)
  {
    for( int j = 0; j < numContactGroups; ++j)
    {
      frictionCoefficientTableCopy[i][j] = m_frictionCoefficientTable[i][j];
    }
  }
  arrayView2d< real64 const > const frictionCoefficientTable = frictionCoefficientTableCopy;

  // Grid fields
  arrayView2d< real64 const > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 const > const gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView2d< real64 const > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 const > const gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView3d< real64 const > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 const > const gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 const > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView2d< real64 const > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
  arrayView3d< real64 const > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );
  arrayView3d< real64 const > const gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 const > const gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );
  arrayView2d< int const > const gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );
  arrayView2d< real64 const > const gridSurfaceNormalWeights = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightsString() );
  arrayView3d< real64 > const gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridContactForceString() );
  arrayView2d< real64 > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    // Initialize gridContactForce[g] to zero. TODO: This shouldn't be necessary?
    // This looks to be zeroed every timestep in initializeGridFields, need to test if removing this breaks anything
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      LvArray::tensorOps::fill< 3 >( gridContactForce[g][fieldIndex], 0.0);
    }

    for( localIndex A = 0; A < numVelocityFields - 1; ++A )
    {
      for( localIndex B = A + 1; B < numVelocityFields; ++B )
      {
        // Make sure both fields in the pair are active
        bool active = ( gridMass[g][A] > smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][A] ) > 1.0e-16 )
                      and
                      ( gridMass[g][B] > smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridSurfaceNormal[g][B] ) > 1.0e-16 ); // CC: Should grid surface normal min magnitude be DBL_MIN instead?

        real64 frictionCoefficient = frictionCoefficientTable[A % numContactGroups][B % numContactGroups];

        if( active )
        {
          // Evaluate the separability criterion for the contact pair.
          int separable = 0;
          int useCohesiveTangentialForces = 0;
          // int zeroTangentialForces 0; // When using contact compressive forces to prevent interpentration, ignore tangential forces // Don't think this is needed if friction coefficient is set to zero
          if( gridCohesiveFieldFlag[g][A] && gridCohesiveFieldFlag[g][B] ) // Eventually this will need to check that both fields correspond to pairs for which a cohesive law is defined
          {
            if( preventCZInterpentration != 1 )
            {
              return;
            }

            separable = 1;
            frictionCoefficient = 0.0;
            useCohesiveTangentialForces = 1;
          }
          else
          {
            separable = evaluateSeparabilityCriterion( planeStrain,
                                                       numContactGroups,
                                                       treatFullyDamagedAsSingleField,
                                                       separabilityMinDamage,
                                                       thinFeatureDFGThreshold,
                                                       neighborRadius,
                                                       A,
                                                       B,
                                                       gridDamage[g][A],
                                                       gridDamage[g][B],
                                                       gridMaxDamage[g][A],
                                                       gridMaxDamage[g][B],
                                                       gridDamageGradient[g],
                                                       gridCenterOfMass[g][A],
                                                       gridCenterOfMass[g][B] );
          }

          computePairwiseNodalContactForce( contactNormalType,
                                            contactGapCorrection,
                                            overlapCorrection,
                                            hEl,
                                            planeStrain,
                                            smallMass,
                                            useSurfacePositionForContact,
                                            useCohesiveTangentialForces,
                                            separable,
                                            dt,
                                            frictionCoefficient,
                                            gridMass[g][A],
                                            gridMass[g][B],
                                            gridMaterialVolume[g][A],
                                            gridMaterialVolume[g][B],
                                            gridVelocity[g][A],
                                            gridVelocity[g][B],
                                            gridMomentum[g][A],
                                            gridMomentum[g][B],
                                            gridSurfaceNormal[g][A],
                                            gridSurfaceNormal[g][B],
                                            gridSurfaceFieldMass[g][A],
                                            gridSurfaceFieldMass[g][B],
                                            gridSurfacePosition[g][A],
                                            gridSurfacePosition[g][B],
                                            gridCenterOfMass[g][A],
                                            gridCenterOfMass[g][B],
                                            gridCenterOfVolume[g][A],
                                            gridCenterOfVolume[g][B],
                                            gridSurfaceNormalWeights[g][A],
                                            gridSurfaceNormalWeights[g][B],
                                            gridContactForce[g][A],
                                            gridContactForce[g][B] );
        }
      }
    }
  } );
}

void SolidMechanicsMPM::computePairwiseNodalContactForce( ContactNormalTypeOption const & contactNormalType,
                                                          ContactGapCorrectionOption const & contactGapCorrection,
                                                          OverlapCorrectionOption const & overlapCorrection,
                                                          real64 const (& hEl) [3],
                                                          int const & planeStrain,
                                                          real64 const & smallMass,
                                                          int const & useSurfacePositionForContact,
                                                          int const & useCohesiveTangentialForces,
                                                          int & separable,
                                                          real64 const & dt,
                                                          real64 const & frictionCoefficient,
                                                          real64 const & mA, // Mass of field A
                                                          real64 const & mB, // Mass of field B
                                                          real64 const & VA,
                                                          real64 const & VB,
                                                          arraySlice1d< real64 const > const vA,
                                                          arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( vB ),
                                                          arraySlice1d< real64 const > const qA,
                                                          arraySlice1d< real64 const > const qB,
                                                          arraySlice1d< real64 const > const nA, // Surface normal of field A
                                                          arraySlice1d< real64 const > const nB, // Surface normal of field B
                                                          real64 const spmA, // Surface field mass of field A
                                                          real64 const spmB, // Surface field mass of field B
                                                          arraySlice1d< real64 const > const spA, // Surface position of field A
                                                          arraySlice1d< real64 const > const spB, // Surface position of field B
                                                          arraySlice1d< real64 const > const xA, // Center of mass of field A
                                                          arraySlice1d< real64 const > const xB, // Center of mass of field B
                                                          arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( centerOfVolumeA ), // Center of volume of field A
                                                          arraySlice1d< real64 const > const GEOS_UNUSED_PARAM( centerOfVolumeB ), // Center of volume of field B
                                                          real64 const & wA, // Surface normal weights of field A
                                                          real64 const & wB, // Surface normal weights of field A
                                                          arraySlice1d< real64 > const fA,
                                                          arraySlice1d< real64 > const fB )
{
  // Total mass for the contact pair.
  real64 mAB = mA + mB;

  // Outward normal of field A with respect to field B.
  real64 nAB[3] = {};

  // contact normal averaging: 0: simple, 1: mass weighted, 2: follow larger mass (useful for platens)
  switch( contactNormalType )
  {
    case ContactNormalTypeOption::Difference:
      LvArray::tensorOps::copy< 3 >( nAB, nA );
      LvArray::tensorOps::subtract< 3 >( nAB, nB );
      break;
    case ContactNormalTypeOption::MassWeighted:
      LvArray::tensorOps::scaledCopy< 3 >( nAB, nA, mA );
      LvArray::tensorOps::scaledAdd< 3 >( nAB, nB, -mB );
      break;
    case ContactNormalTypeOption::LargerMass:
      if( mA > mB )
      {
        LvArray::tensorOps::copy< 3 >( nAB, nA );
      }
      else
      {
        LvArray::tensorOps::scaledCopy< 3 >( nAB, nB, -1 );
      }
      break;
    case ContactNormalTypeOption::Mixed:
      {
        // If density is similar use mass weighted average
        real64 rhoA = mA / VA;
        real64 rhoB = mB / VB;
        if ( LvArray::math::abs( rhoA - rhoB ) < 0.1 * rhoA )
        {
          for( int i=0; i<3; ++i )
          {
            LvArray::tensorOps::scaledCopy< 3 >( nAB, nA, mA );
            LvArray::tensorOps::scaledAdd< 3 >( nAB, nB, -mB );
          }
        }
        else
        {
          // if one field is significantly more dense,
          // Use the surface normal for whichever field has higher density.
          // This should be good for corners against flat surfaces.
          if( rhoA > rhoB )
          {
            LvArray::tensorOps::copy< 3 >( nAB, nA );
          }
          else
          {
            LvArray::tensorOps::scaledCopy< 3 >( nAB, nB, -1 );
          }
        }
        break;
      }
    case ContactNormalTypeOption::Aligned:
      {
        real64 tempA[3];
        real64 tempB[3];
        // LvArray::tensorOps::scaledCopy< 3 >( tempA, nA, pow(wA, m_contactNormalExponent) );
        // LvArray::tensorOps::scaledCopy< 3 >( tempB, nB, pow(wB, m_contactNormalExponent) );
        real64 threshold = 0.9;
        real64 xxA = LvArray::math::min(LvArray::math::max((wA-threshold)/threshold, 0.0), 1.0);
        real64 xxB = LvArray::math::min(LvArray::math::max((wB-threshold)/threshold, 0.0), 1.0);
        LvArray::tensorOps::scaledCopy< 3 >(tempA, nA, 3*LvArray::math::pow(xxA,2)-2*LvArray::math::pow(xxA,3));
        LvArray::tensorOps::scaledCopy< 3 >(tempB, nB, 3*LvArray::math::pow(xxB,2)-2*LvArray::math::pow(xxB,3));
        LvArray::tensorOps::copy< 3 >( nAB, tempA );
        LvArray::tensorOps::subtract< 3 >( nAB, tempB );
      }
      break;
    default:
      GEOS_ERROR( "Unrecognized contact normal type!" );
      break;
  }

  // Vector pointing from CoM of A field to CoM of B field, stretched by grid spacing
  // for( int i=0; i<3; ++i )
  // {
  //   nAB[i] = (xB[i] - xA[i]) / m_hEl[i];
  // }

  // Normalize the effective surface normal
  if( planeStrain == 1 )
  {
    nAB[2] = 0.0;
  }
  real64 norm = LvArray::math::sqrt( nAB[0] * nAB[0] + nAB[1] * nAB[1] + nAB[2] * nAB[2] );

  // // CC: Is this the best solution?
  // // Need to prevent crashing from randomly defined normals that might arise when running with damage field gradient partitioning
  // // If the normal is small, then assume no separation force and return
  // if ( norm < 1e-20 )
  // {
  //     fA[0] = 0.0;
  //     fA[1] = 0.0;
  //     fA[2] = 0.0;
  //     fB[0] = 0.0;
  //     fB[1] = 0.0;
  //     fB[2] = 0.0;
  //     return;
  // }

  // nAB[0] /= norm;
  // nAB[1] /= norm;
  // nAB[2] /= norm;

  if( norm > 1e-20 )
  {
    LvArray::tensorOps::scale< 3 >( nAB, 1 / norm );
  }
  else
  {
    // If normals are randomly defined as in the case of a fully damaged region, just default to normal A
    // since the two fields should not be separable
    LvArray::tensorOps::copy< 3 >( nAB, nA );
    nAB[2] = 0.0;
  }

  // Total momentum for the contact pair.
  real64 qAB[3];
  qAB[0] = qA[0] + qB[0];
  qAB[1] = qA[1] + qB[1];
  qAB[2] = qA[2] + qB[2];

  // Center-of-mass velocity for the contact pair.
  real64 vAB[3];
  vAB[0] = qAB[0] / mAB;
  vAB[1] = qAB[1] / mAB;
  vAB[2] = qAB[2] / mAB;

  // Compute s1AB and s2AB, to form an orthonormal basis. This uses the method by E. Herbold
  // to ensure consistency between surfaces.
  real64 s1AB[3], s2AB[3]; // Tangential vectors for the contact pair
  computeOrthonormalBasis( nAB, s1AB, s2AB );

  // Compute force decomposition, declare increment in fA from "this" contact pair
  real64 fnor =  ( mA / dt ) * ( (vAB[0] - vA[0]) * nAB[0]  + (vAB[1] - vA[1]) * nAB[1]  + (vAB[2] - vA[2]) * nAB[2] ),
         ftan1 = ( mA / dt ) * ( (vAB[0] - vA[0]) * s1AB[0] + (vAB[1] - vA[1]) * s1AB[1] + (vAB[2] - vA[2]) * s1AB[2] ),
         ftan2 = ( mA / dt ) * ( (vAB[0] - vA[0]) * s2AB[0] + (vAB[1] - vA[1]) * s2AB[1] + (vAB[2] - vA[2]) * s2AB[2] );
  real64 dfA[3];

  // // CC: to fix field partitioning with explicit surface normals at sharp tips such as nanoindenter tip
  // // TODO: confirm this works
  // real64 dC[3] = {};
  // LvArray::tensorOps::copy< 3 >( dC, centerOfVolumeA );
  // LvArray::tensorOps::subtract< 3 >(dC, centerOfVolumeB );
  // separable |= LvArray::tensorOps::AiBi< 3 >(dC, nAB) > 0;

  // Check for separability, and enforce either slip, or no-slip contact, accordingly
  if( separable == 0 )
  {
    // Surfaces are bonded, treat as single velocity field by applying normal force to prevent
    // interpenetration, and tangential force to prevent slip.
    dfA[0] = fnor * nAB[0] + ftan1 * s1AB[0] + ftan2 * s2AB[0];
    dfA[1] = fnor * nAB[1] + ftan1 * s1AB[1] + ftan2 * s2AB[1];
    dfA[2] = fnor * nAB[2] + ftan1 * s1AB[2] + ftan2 * s2AB[2];
    fA[0] += dfA[0];
    fA[1] += dfA[1];
    fA[2] += dfA[2];
    fB[0] -= dfA[0];
    fB[1] -= dfA[1];
    fB[2] -= dfA[2];
    return;
  }

  // Surfaces are separable. For frictional contact, apply a normal force to
  // prevent interpenetration, and tangential force to prevent slip unless f_tan > mu*f_nor

  // Calculate the contact gap between the fields
  real64 gap0;
  if( planeStrain == 1 )
  {
    // gap0 = ( hEl[0]*hEl[1] ) / LvArray::math::sqrt( hEl[1]*hEl[1]*nAB[0]*nAB[0] + hEl[0]*hEl[0]*nAB[1]*nAB[1] );
    gap0 = 1/LvArray::math::sqrt( LvArray::math::pow(nAB[0]/hEl[0],2) + LvArray::math::pow(nAB[1]/hEl[1],2) );
  }
  else
  {
    // CC: debug
    // currently failes with polymer model (suspect compliant materials) when z direction boundary type is 2
    //This needs to be corrected
    // gap0 = (hEl[0]*hEl[1]*hEl[2]) /
    //        LvArray::math::sqrt( hEl[2]*hEl[2]*nAB[2]*nAB[2]*( hEl[1]*hEl[1]*nAB[0]*nAB[0] + hEl[0]*hEl[0]*nAB[1]*nAB[1] ) + hEl[0]*hEl[0]*hEl[1]*hEl[1]*( nAB[0]*nAB[0] + nAB[1]*nAB[1] ) );

    // Elliptical solution
    // Gives element sizes in each direction and reasonable approximations in others
    gap0 = 1/LvArray::math::sqrt(LvArray::math::pow(nAB[0]/hEl[0],2) + LvArray::math::pow(nAB[1]/hEl[1],2) + LvArray::math::pow(nAB[2]/hEl[2],2));

  }

  // TODO: A fudge factor of 0.67 on gap0 makes diagonal surfaces (wrt grid) close better I think, but this is more general
  // real64 gap = (xB[0] - xA[0]) * nAB[0] + (xB[1] - xA[1]) * nAB[1] + (xB[2] - xA[2]) * nAB[2] - gap0;

  // In case of asymmetric interfaces (explicit normals and positions only on one side) use the center of mass instead of surface position
  real64 gapScale = 0.0;
  real64 surfacePosA[3] = {};
  if ( spmA > smallMass && useSurfacePositionForContact )
  {
    LvArray::tensorOps::copy< 3 >( surfacePosA, spA);
  }
  else
  {
    LvArray::tensorOps::copy< 3 >( surfacePosA, xA );
    gapScale += 0.5;
  }

  real64 surfacePosB[3] = {};
  if( spmB > smallMass && useSurfacePositionForContact )
  {
    LvArray::tensorOps::copy< 3 >( surfacePosB, spB);
  }
  else
  {
    LvArray::tensorOps::copy< 3 >( surfacePosB, xB );
    gapScale += 0.5;
  }
    // gap = LvArray::tensorOps::AiBi< 3 >(nAB, surfacePosB) - LvArray::tensorOps::AiBi< 3 >(nAB, surfacePosA) - gapScale*gap0;
  // }

  real64 gap = (surfacePosB[0] - surfacePosA[0]) * nAB[0] + (surfacePosB[1] - surfacePosA[1]) * nAB[1] + (surfacePosB[2] - surfacePosA[2]) * nAB[2] - gapScale*gap0;

  real64 contact = 0.0;
  real64 temp[3] = {};
  LvArray::tensorOps::copy< 3 >( temp, vA );
  LvArray::tensorOps::subtract< 3 >( temp, vAB );
  real64 test = LvArray::tensorOps::AiBi< 3 >( temp, nAB);

  // real64 test = (vA[0] - vAB[0]) * nAB[0] + (vA[1] - vAB[1]) * nAB[1] + (vA[2] - vAB[2]) * nAB[2];

  switch( contactGapCorrection )
  {
    // Simplely check if component of field velocities will result in interpenetration if uncorrected
    case ContactGapCorrectionOption::Simple:
      contact = test > 0.0 ? 1.0 : 0.0;
      break;
    // Check materials are interpenetrating and velocities will results in further interpenetration
    case ContactGapCorrectionOption::Implicit:
      contact = test > 0.0 && gap < 0.0 ? 1.0 : 0.0;
      break;
    //
    case ContactGapCorrectionOption::Softened:
      if ( test > 0 )
      {
        // realT gap0 = normalSpacing*cellSpacing;
        if (gap <= 0.0)
        {
          contact = 1.0;
        }
        else if (gap < gap0)
        {
          contact = 1.0 - gap/gap0;
        }
      }
      break;
    default:
      GEOS_ERROR( "Unknown contact gap correction type specified" );
      break;
  }

  // Modify normal contact force
  fnor *= contact;

  // Additional compressive normal force added to correct overlap, this is added only to the normal force and
  // doesn't affect the shear calculation (which in FEM cuts down on spurious stress oscillations).
  real64 fgap = 0.0;
  if( overlapCorrection == OverlapCorrectionOption::NormalForce )
  {
    real64 cellSpacing[3];
    LvArray::tensorOps::copy< 3 >( cellSpacing, hEl);

    real64 normalSpacing[3];
    normalSpacing[0] = LvArray::math::abs( nAB[0] );
    normalSpacing[1] = LvArray::math::abs( nAB[1] );
    normalSpacing[2] = LvArray::math::abs( nAB[2] );

    real64 cellVolume = hEl[0] * hEl[1] * hEl[2];
    real64 cellLength = LvArray::tensorOps::AiBi< 3 >( normalSpacing, cellSpacing );
    real64 cellArea = cellVolume / cellLength;
    real64 overlap = planeStrain ? ( VA + VB - 0.5*cellVolume ) / cellArea : ( VA + VB - cellVolume ) / cellArea;
    if( overlap > 0.0 )
    {
      fgap = -1.0 * overlap * fmin( mA, mB ) / ( dt * dt ); // This could be -0.5*overlap*mAB/dt**2, but that might be less stable if one mass tiny
    }
  }

  // Determine force for tangential sticking
  real64 ftanMag = sqrt( ftan1 * ftan1 + ftan2 * ftan2 );

  // Get direction of tangential contact force
  real64 sAB[3];
  if( ftanMag > 0.0 )
  {
    sAB[0] = (s1AB[0] * ftan1 + s2AB[0] * ftan2) / ftanMag;
    sAB[1] = (s1AB[1] * ftan1 + s2AB[1] * ftan2) / ftanMag;
    sAB[2] = (s1AB[2] * ftan1 + s2AB[2] * ftan2) / ftanMag;
  }
  else
  {
    sAB[0] = 0.0;
    sAB[1] = 0.0;
    sAB[2] = 0.0;
  }

  if( useCohesiveTangentialForces == 1)
  {
    ftan1 = 0.0;
    ftan2 = 0.0;
  }

  real64 ftan = LvArray::math::min( frictionCoefficient * LvArray::math::abs( fnor ), ftanMag ); // This goes to zero when contact=0 due to the std::min
  dfA[0] = ( fnor + fgap ) * nAB[0] + ftan * sAB[0];
  dfA[1] = ( fnor + fgap ) * nAB[1] + ftan * sAB[1];
  dfA[2] = ( fnor + fgap ) * nAB[2] + ftan * sAB[2];
  fA[0] += dfA[0];
  fA[1] += dfA[1];
  fA[2] += dfA[2];
  fB[0] -= dfA[0];
  fB[1] -= dfA[1];
  fB[2] -= dfA[2];
}

GEOS_HOST_DEVICE
void SolidMechanicsMPM::computeOrthonormalBasis( const real64 * e1, // input "normal" unit vector.
                                                 real64 * e2,       // output "tangential" unit vector.
                                                 real64 * e3 )      // output "tangential" unit vector.
{
  // This routine takes in a normalized vector and gives two orthogonal vectors.
  // It is rather arbitrary, in general, how this is done; however, this routine
  // is written in a consistent way. This is important, for example, if you
  // have a face on one processor and a ghost face on another processor and
  // you try to define a basis that is EXACTLY the same for debugging purposes.
  // This came up when trying to debug shear traction components and this
  // routine helped a lot.

  real64 e1x = e1[0],
         e1y = e1[1],
         e1z = e1[2],
         e2x,
         e2y,
         e2z;

  // find the maximum pair of values from the first vector.
  real64 maxp1 = LvArray::math::abs( e1x ) + LvArray::math::abs( e1y ),
         maxp2 = LvArray::math::abs( e1y ) + LvArray::math::abs( e1z ),
         maxp3 = LvArray::math::abs( e1x ) + LvArray::math::abs( e1z );

  // This part amounts to setting the smallest component of the input vector to 0.
  // Then, the orthogonal vector is a simple operation in the 2D plane.
  if( maxp1 >= maxp2 && maxp1 >= maxp3 )
  {
    e2x = e1y;
    e2y = -e1x;
    e2z = 0.0;
  }
  else if( maxp2 >= maxp1 && maxp2 >= maxp3 )
  {
    e2y = e1z;
    e2z = -e1y;
    e2x = 0.0;
  }
  else
  {
    e2z = e1x;
    e2x = -e1z;
    e2y = 0.0;
  }

  // Normalize the first base vector
  real64 e2norm = LvArray::math::sqrt( e2x * e2x + e2y * e2y + e2z * e2z );
  e2[0] = e2x / e2norm;
  e2[1] = e2y / e2norm;
  e2[2] = e2z / e2norm;

  // The second base vector is simply a cross product of the first two.  For
  // safety sake, this may need to be normalized to combat issues of numerical
  // precision.
  e3[0] =  e1y * e2z - e1z * e2y;
  e3[1] = -e1x * e2z + e1z * e2x;
  e3[2] =  e1x * e2y - e1y * e2x;

  real64 e3norm = LvArray::math::sqrt( e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2] );
  e3[0] /= e3norm;
  e3[1] /= e3norm;
  e3[2] /= e3norm;
}

void SolidMechanicsMPM::setGridFieldLabels( NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // TODO: Find a way to automatically loop over these fields.
  // I feel like a dimension check would work.
  // Generate labels
  stdVector< std::string > fieldLabels( m_numVelocityFields );
  std::generate( fieldLabels.begin(), fieldLabels.end(), [i=0]() mutable { return "velocityField" + std::to_string( ++i ); } );
  string const axesLabels[] = { "X", "Y", "Z" };
  string const tensorLabels[] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

  // Apply labels to scalar multi-fields
  stdVector< std::string > keys2d = {
                                        viewKeyStruct::gridMassString(),
                                        viewKeyStruct::gridDamageString(),
                                        viewKeyStruct::gridMaxDamageString(),
                                        viewKeyStruct::gridMaterialVolumeString(),
                                        viewKeyStruct::gridSurfaceFieldMassString(),
                                        // viewKeyStruct::gridSurfaceNormalWeightsString(),
                                        // viewKeyStruct::gridSurfaceNormalWeightNormalizationString(),
                                        // viewKeyStruct::gridMassWeightedDamageString(),
                                      };
  for( auto const & key: keys2d )
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array2d< real64 > >( key );
    wrapper.setDimLabels( 1, fieldLabels );
  }

  nodeManager.getWrapper< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() ).setDimLabels( 1, fieldLabels );
  nodeManager.getWrapper< array2d< real64 > >( viewKeyStruct::gridPrincipalExplicitSurfaceNormalString() ).setDimLabels( 1, axesLabels );
  nodeManager.getWrapper< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() ).setDimLabels( 1, axesLabels );

  // Apply labels to vector multi-fields
  stdVector< std::string > keys3d = { viewKeyStruct::gridVelocityString(), // 3
                                      viewKeyStruct::gridDVelocityString(), // 3
                                      viewKeyStruct::gridDisplacementString(), // 3
                                      viewKeyStruct::gridMomentumString(), // 3
                                      viewKeyStruct::gridAccelerationString(), // 3
                                      viewKeyStruct::gridInternalForceString(), // 3
                                      viewKeyStruct::gridExternalForceString(), // 3
                                      viewKeyStruct::gridContactForceString(), // 3
                                      viewKeyStruct::gridSurfaceNormalString(), // 3
                                      viewKeyStruct::gridCenterOfMassString(), // 3
                                      viewKeyStruct::gridCenterOfVolumeString(), // 3
                                      viewKeyStruct::gridSurfaceNormalString(),
                                      viewKeyStruct::gridSurfacePositionString(),
                                      // viewKeyStruct::gridExplicitSurfaceNormalString(),
                                      viewKeyStruct::gridReferenceAreaVectorString(),
                                      viewKeyStruct::gridReferenceSurfacePositionString(),
                                      // viewKeyStruct::gridReferenceMaterialVolumeString(),
                                      viewKeyStruct::gridCohesiveAreaString(),
                                      viewKeyStruct::gridCohesiveForceString()
                                      };

  for( auto const & key: keys3d )
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array3d< real64 > >( key );
    wrapper.setDimLabels( 1, fieldLabels );
    wrapper.setDimLabels( 2, axesLabels );
  }

  stdVector< std::string > keysVoigt = {
                                          viewKeyStruct::gridNormalStressString()
                                        };

  for( auto const & key: keysVoigt )
  {
    WrapperBase & wrapper = nodeManager.getWrapper< array3d< real64 > >( key );
    wrapper.setDimLabels( 1, fieldLabels );
    wrapper.setDimLabels( 2, axesLabels );
  }
}

void SolidMechanicsMPM::resizeGrid( SpatialPartition & partition,
                                    NodeManager & nodeManager,
                                    real64 const dt )
{
  GEOS_MARK_FUNCTION;

  // Modify SpatialPartition class members
  partition.updateSizes( m_domainL, dt );

  // Modify SolidMechanicsMPM class members
  real64 ratio[3];
  for( int i = 0; i < 3; ++i )
  {
    // Incremental stretch
    ratio[i] = 1.0 + m_domainL[i] * dt;

    // Modify SolidMechanicsMPM class members
    m_hEl[i] *= ratio[i];
    m_xLocalMin[i] *= ratio[i];
    m_xLocalMax[i] *= ratio[i];
    m_xLocalMinNoGhost[i] *= ratio[i];
    m_xLocalMaxNoGhost[i] *= ratio[i];
    m_xGlobalMin[i] *= ratio[i];
    m_xGlobalMax[i] *= ratio[i];
    m_partitionExtent[i] *= ratio[i];
    m_domainExtent[i] *= ratio[i];
  }

  // Update nodal positions
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  forAll< serialPolicy >( gridPosition.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    gridPosition[g][0] *= ratio[0];
    gridPosition[g][1] *= ratio[1];
    gridPosition[g][2] *= ratio[2];
  } );
}

void SolidMechanicsMPM::solverProfiling( std::string label )
{
  if( m_solverProfiling == 1 )
  {
    MPI_Barrier( MPI_COMM_GEOS );
    m_profilingTimes.push_back( MPI_Wtime() );
    m_profilingLabels.push_back( label );
  }
}

void SolidMechanicsMPM::setParticlesConstitutiveNames( ParticleSubRegionBase & subRegion ) const
{
  subRegion.registerWrapper< string >( viewKeyStruct::solidMaterialNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  // Check if this should be changed to constitutive base since we have introduced the gas material which does not inherit from solidBase
  string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
  solidMaterialName = getConstitutiveName< SolidBase >( subRegion );
  GEOS_ERROR_IF( solidMaterialName.empty(), GEOS_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}

void SolidMechanicsMPM::setConstitutiveNamesCallSuper( ParticleSubRegionBase & subRegion ) const
{
  PhysicsSolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( viewKeyStruct::solidMaterialNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
  solidMaterialName = PhysicsSolverBase::getConstitutiveName< ContinuumBase >( subRegion );
  GEOS_ERROR_IF( solidMaterialName.empty(), GEOS_FMT( "ContinuumBase model not found on subregion {}", subRegion.getName() ) );
}

// Eventually the bin sorting subroutine from computeNeighborList should be broken out into it's own function
// void SolidMechanicsMPM::particleBinSort()
// {
//   // Expand bin limits by neighbor radius to account for the buffer zone of ghost particles outside the patch limits
//   real64 const neighborRadius = m_neighborRadius;
//   real64 const neighborRadiusSquared = neighborRadius * neighborRadius;
//   real64 xmin = m_xLocalMinNoGhost[0] - neighborRadius,
//          xmax = m_xLocalMaxNoGhost[0] + neighborRadius,
//          ymin = m_xLocalMinNoGhost[1] - neighborRadius,
//          ymax = m_xLocalMaxNoGhost[1] + neighborRadius,
//          zmin = m_xLocalMinNoGhost[2] - neighborRadius,
//          zmax = m_xLocalMaxNoGhost[2] + neighborRadius;
  
//   // Initialize bin sort
//   real64 binWidth = m_binSizeMultiplier * neighborRadius;
//   int nxbins = LvArray::math::ceil( ( xmax - xmin ) / binWidth ),
//       nybins = LvArray::math::ceil( ( ymax - ymin ) / binWidth ),
//       nzbins = m_planeStrain ? 1 : LvArray::math::ceil( ( zmax - zmin ) / binWidth );
//   int nbins = nxbins * nybins * nzbins;
//   real64 dx = ( xmax - xmin ) / nxbins,
//          dy = ( ymax - ymin ) / nybins,
//          dz = ( zmax - zmin ) / nzbins;

//   localIndex totalNumberOfParticles = particleManager.getNumberOfParticles();
//   localIndex totalNumberOfBins = nbins * m_numberOfSubRegions;

//   // Initializes array with totalNumberOfBins inner arrays and default sizes 0 (first dimension is the bins, second are particles indices per bin)
//   ArrayOfArrays< localIndex > bins(totalNumberOfBins, 0); 
  
//   // Precompute number of particles in each bin and resize the bins arrayOfArrays

//   // ParticleSubRegions are not access in gpu kernels so we must make copies of all relevant fields
//   // such as the particle centers for neighbor checking
//   array2d< real64 > allParticleCenters(totalNumberOfParticles, 3);

//   // TODO: We create a view here, but we may want to have views with different access restrictions (e.g. const) separately in each step
//   arrayView2d< real64 > const allParticleCentersView = allParticleCenters.toView();
  
//   array1d< localIndex > subRegionSizes(m_numberOfSubRegions);
//   array1d< localIndex > regionIndicesOfSubRegions(m_numberOfSubRegions);
//   array1d< localIndex > subRegionIndicesInRegions(m_numberOfSubRegions);

//   arrayView1d< localIndex > const subRegionSizesView = subRegionSizes.toView();
//   arrayView1d< localIndex > const regionIndicesOfSubRegionsView = regionIndicesOfSubRegions.toView();
//   arrayView1d< localIndex > const subRegionIndicesInRegionsView = subRegionIndicesInRegions.toView();
//   localIndex particleIndexOffset = 0;
//   localIndex subRegionIndex = 0;
//   particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
//   {
//     // arrayView1d< int const > const particleRank = subRegion.getParticleRank();

//     RAJA::MultiReduceSum< serialMultiReduce, localIndex > binSizeReduction( nbins, 0 );
   
//     arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();

//     // GEOS_LOG_RANK("binSizeReductionView.size(): " << binSizeReductionView.size() << ", nbins: " << nbins << ", xmin: {"<< xmin << ", " << ymin << ", " << zmin << "}"  << ", xmax: {"<< xmax << ", " << ymax << ", " << zmax << "}");
//     // GEOS_LOG_RANK("nbins: " << nbins << ", xmin: {"<< xmin << ", " << ymin << ", " << zmin << "}"  << ", xmax: {"<< xmax << ", " << ymax << ", " << zmax << "}");
    
//     // SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
//     // SortedArrayView< localIndex const > const inactiveParticleIndices = subRegion.inactiveParticleIndices();

//     // When compiling on ruby running with sequential reduction policy, race condition in bin size counting causes memory issue
//     // Fix add export OMP_NUM_THREADS=1 before srun, need to check compilation variables to ensure GEOS_USE_OPENMP is disabled on ruby
//     // forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
//     forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST ( localIndex const p )
//     {
//       // Particle bin ijk indices
//       localIndex i = LvArray::math::floor( ( particlePosition[p][0] - xmin ) / dx );
//       localIndex j = LvArray::math::floor( ( particlePosition[p][1] - ymin ) / dy );
//       localIndex k = LvArray::math::floor( ( particlePosition[p][2] - zmin ) / dz );
//       localIndex binIndex = i + j * nxbins + k * nxbins * nybins;
      
//       binSizeReduction[binIndex] += 1;
//       // if( binIndex > nbins ) 
//       // {
//       //   GEOS_LOG_LEVEL_BY_RANK(2, "binIndex: " << binIndex << ", i,j,k: " << i << ", " << j << ", " <<  k <<  ", nbins: " << nbins << ", p: " << p << ", pos: {" << particlePosition[p][0] << ", "  << particlePosition[p][1] << ", " << particlePosition[p][2] << "}" << ", pRank: " << particleRank[p] << ", active: " << activeParticleIndices.contains(p) << ", inactive:" << inactiveParticleIndices.contains(p) );
//       // }

//       // Copy particle data for neighbor search later
//       LvArray::tensorOps::copy< 3 >( allParticleCentersView[particleIndexOffset + p], particlePosition[p] );
//     } );
//     particleIndexOffset += subRegion.size();

//     // Not sure this would be worth doing in parallel (also resizing is not done on device)
//     for( int b = 0; b < nbins; ++b )
//     {
//       bins.resizeArray( subRegionIndex * nbins + b,
//                         binSizeReduction[b].get() );
//     }
    
//     ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
 
//     regionIndicesOfSubRegionsView[subRegionIndex] = region.getIndexInParent();
//     subRegionIndicesInRegionsView[subRegionIndex] = subRegion.getIndexInParent();
//     subRegionSizesView[subRegionIndex] = subRegion.size();
//     ++subRegionIndex;

//     allParticleCentersView.move( LvArray::MemorySpace::host ); // Must move the particle center data back from device explicitly
//   } );

//   // Populate bins with particle data

//   // Stored the current count of particles per bin during population
//   array1d< localIndex > binCount(totalNumberOfBins);

//   // Initialize bin count to zero
//   for( localIndex b = 0; b < totalNumberOfBins; ++b )
//   {
//     binCount[b] = 0;
//   }

//   ArrayOfArraysView< localIndex > binsView = bins.toView();

//   solverProfiling( "Populate bins..." );
//   GEOS_LOG_LEVEL_BY_RANK(2, "Populate bins..." );

//   subRegionIndex = 0;
//   particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
//   {
//     arrayView1d< localIndex > const binCountView = binCount.toView();

//     arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
//     // Running this in parallel threading over particles would result in race condition during particle index assignment to bins
//     // We could introduce atomics, but it's unclear whether there would be any substantial performance improvement, should still test
//     // Alternatively we could thread over bins, but it may be faster to just compute in serial since each bin would have to check all particles in subregion
//     forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
//     {
//       // Particle bin ijk indices
//       localIndex i = LvArray::math::floor( ( particlePosition[p][0] - xmin ) / dx );
//       localIndex j = LvArray::math::floor( ( particlePosition[p][1] - ymin ) / dy );
//       localIndex k = LvArray::math::floor( ( particlePosition[p][2] - zmin ) / dz );
//       localIndex binIndex = nbins * subRegionIndex + i + j * nxbins + k * nxbins * nybins;
      
//       binsView[binIndex][binCountView[binIndex]] = p;
//       ++binCountView[binIndex];
//     } );

//     binsView.move( LvArray::MemorySpace::host );
//     binCountView.move( LvArray::MemorySpace::host );

//     ++subRegionIndex;
//   } );
// }

real64 SolidMechanicsMPM::computeNeighborList( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // Time this function
  real64 tStart = MPI_Wtime();

  // Expand bin limits by neighbor radius to account for the buffer zone of ghost particles outside the patch limits
  real64 neighborRadius = m_neighborRadius;
  real64 neighborRadiusSquared = neighborRadius * neighborRadius;
  real64 xmin = m_xLocalMinNoGhost[0] - neighborRadius,
         xmax = m_xLocalMaxNoGhost[0] + neighborRadius,
         ymin = m_xLocalMinNoGhost[1] - neighborRadius,
         ymax = m_xLocalMaxNoGhost[1] + neighborRadius,
         zmin = m_xLocalMinNoGhost[2] - neighborRadius,
         zmax = m_xLocalMaxNoGhost[2] + neighborRadius;
  
  // Initialize bin sort
  real64 binWidth = m_binSizeMultiplier * neighborRadius;
  int nxbins = LvArray::math::ceil( ( xmax - xmin ) / binWidth ),
      nybins = LvArray::math::ceil( ( ymax - ymin ) / binWidth ),
      nzbins = m_planeStrain ? 1 : LvArray::math::ceil( ( zmax - zmin ) / binWidth );
  int nbins = nxbins * nybins * nzbins;
  real64 dx = ( xmax - xmin ) / nxbins,
         dy = ( ymax - ymin ) / nybins,
         dz = ( zmax - zmin ) / nzbins;

  localIndex totalNumberOfParticles = particleManager.getNumberOfParticles();
  localIndex totalNumberOfBins = nbins * m_numberOfSubRegions;

  // Initializes array with totalNumberOfBins inner arrays and default sizes 0 (first dimension is the bins, second are particles indices per bin)
  ArrayOfArrays< localIndex > bins(totalNumberOfBins, 0); 
  
  // Precompute number of particles in each bin and resize the bins arrayOfArrays

  // ParticleSubRegions are not access in gpu kernels so we must make copies of all relevant fields
  // such as the particle centers for neighbor checking
  array2d< real64 > allParticleCenters(totalNumberOfParticles, 3);

  // TODO: We create a view here, but we may want to have views with different access restrictions (e.g. const) separately in each step
  arrayView2d< real64 > const allParticleCentersView = allParticleCenters.toView();
  
  array1d< localIndex > subRegionSizes(m_numberOfSubRegions);
  array1d< localIndex > regionIndicesOfSubRegions(m_numberOfSubRegions);
  array1d< localIndex > subRegionIndicesInRegions(m_numberOfSubRegions);

  // solverProfiling( "Count bin sizes..." );
  // GEOS_LOG_LEVEL_BY_RANK(2, "Count bin sizes..." );

  arrayView1d< localIndex > const subRegionSizesView = subRegionSizes.toView();
  arrayView1d< localIndex > const regionIndicesOfSubRegionsView = regionIndicesOfSubRegions.toView();
  arrayView1d< localIndex > const subRegionIndicesInRegionsView = subRegionIndicesInRegions.toView();
  localIndex particleIndexOffset = 0;
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // arrayView1d< int const > const particleRank = subRegion.getParticleRank();

    RAJA::MultiReduceSum< serialMultiReduce, localIndex > binSizeReduction( nbins, 0 );
   
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();

    // GEOS_LOG_RANK("binSizeReductionView.size(): " << binSizeReductionView.size() << ", nbins: " << nbins << ", xmin: {"<< xmin << ", " << ymin << ", " << zmin << "}"  << ", xmax: {"<< xmax << ", " << ymax << ", " << zmax << "}");
    // GEOS_LOG_RANK("nbins: " << nbins << ", xmin: {"<< xmin << ", " << ymin << ", " << zmin << "}"  << ", xmax: {"<< xmax << ", " << ymax << ", " << zmax << "}");
    
    // SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    // SortedArrayView< localIndex const > const inactiveParticleIndices = subRegion.inactiveParticleIndices();

    // When compiling on ruby running with sequential reduction policy, race condition in bin size counting causes memory issue
    // Fix add export OMP_NUM_THREADS=1 before srun, need to check compilation variables to ensure GEOS_USE_OPENMP is disabled on ruby
    // forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
    forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST ( localIndex const p )
    {
      // Particle bin ijk indices
      localIndex i = LvArray::math::floor( ( particlePosition[p][0] - xmin ) / dx );
      localIndex j = LvArray::math::floor( ( particlePosition[p][1] - ymin ) / dy );
      localIndex k = LvArray::math::floor( ( particlePosition[p][2] - zmin ) / dz );
      localIndex binIndex = i + j * nxbins + k * nxbins * nybins;
      
      binSizeReduction[binIndex] += 1;
      // if( binIndex > nbins ) 
      // {
      //   GEOS_LOG_LEVEL_BY_RANK(2, "binIndex: " << binIndex << ", i,j,k: " << i << ", " << j << ", " <<  k <<  ", nbins: " << nbins << ", p: " << p << ", pos: {" << particlePosition[p][0] << ", "  << particlePosition[p][1] << ", " << particlePosition[p][2] << "}" << ", pRank: " << particleRank[p] << ", active: " << activeParticleIndices.contains(p) << ", inactive:" << inactiveParticleIndices.contains(p) );
      // }

      // Copy particle data for neighbor search later
      LvArray::tensorOps::copy< 3 >( allParticleCentersView[particleIndexOffset + p], particlePosition[p] );
    } );
    particleIndexOffset += subRegion.size();

    // Not sure this would be worth doing in parallel (also resizing is not done on device)
    for( int b = 0; b < nbins; ++b )
    {
      bins.resizeArray( subRegionIndex * nbins + b,
                        binSizeReduction[b].get() );
    }
    
    ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
 
    regionIndicesOfSubRegionsView[subRegionIndex] = region.getIndexInParent();
    subRegionIndicesInRegionsView[subRegionIndex] = subRegion.getIndexInParent();
    subRegionSizesView[subRegionIndex] = subRegion.size();
    ++subRegionIndex;

    allParticleCentersView.move( LvArray::MemorySpace::host ); // Must move the particle center data back from device explicitly
  } );

  // Populate bins with particle data

  // Stored the current count of particles per bin during population
  array1d< localIndex > binCount(totalNumberOfBins);

  // Initialize bin count to zero
  for( localIndex b = 0; b < totalNumberOfBins; ++b )
  {
    binCount[b] = 0;
  }

  ArrayOfArraysView< localIndex > binsView = bins.toView();

  subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< localIndex > const binCountView = binCount.toView();

    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    // Running this in parallel threading over particles would result in race condition during particle index assignment to bins
    // We could introduce atomics, but it's unclear whether there would be any substantial performance improvement, should still test
    // Alternatively we could thread over bins, but it may be faster to just compute in serial since each bin would have to check all particles in subregion
    forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
    {
      // Particle bin ijk indices
      localIndex i = LvArray::math::floor( ( particlePosition[p][0] - xmin ) / dx );
      localIndex j = LvArray::math::floor( ( particlePosition[p][1] - ymin ) / dy );
      localIndex k = LvArray::math::floor( ( particlePosition[p][2] - zmin ) / dz );
      localIndex binIndex = nbins * subRegionIndex + i + j * nxbins + k * nxbins * nybins;
      
      binsView[binIndex][binCountView[binIndex]] = p;
      ++binCountView[binIndex];
    } );

    binsView.move( LvArray::MemorySpace::host );
    binCountView.move( LvArray::MemorySpace::host );

    ++subRegionIndex;
  } );

  // Precompute number of neighbors for each particles
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegionA )
  {
    // Get 'this' particle's location
    arrayView2d< real64 const > const xA = subRegionA.getParticleCenter();

    // Find neighbors of all particles and count them
    SortedArrayView< localIndex const > const subRegionAActiveParticleIndices = subRegionA.activeParticleIndices();
 
    // Create an array to store neighbor counts for all particles
    array1d< localIndex > neighborCounts(subRegionAActiveParticleIndices.size());
    arrayView1d< localIndex > const neighborCountsView = neighborCounts.toView();

    forAll< parallelDevicePolicy<> >( subRegionAActiveParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      // Locating neighbors of particle with index a
      localIndex a = subRegionAActiveParticleIndices[pp];

      // Bin ijk indices bounding a sphere of radius m_neighborRadius centered at 'this' particle
      localIndex imin = LvArray::math::floor( ( xA[a][0] - neighborRadius - xmin ) / dx );
      localIndex jmin = LvArray::math::floor( ( xA[a][1] - neighborRadius - ymin ) / dy );
      localIndex kmin = LvArray::math::floor( ( xA[a][2] - neighborRadius - zmin ) / dz );
      localIndex imax = LvArray::math::floor( ( xA[a][0] + neighborRadius - xmin ) / dx );
      localIndex jmax = LvArray::math::floor( ( xA[a][1] + neighborRadius - ymin ) / dy );
      localIndex kmax = LvArray::math::floor( ( xA[a][2] + neighborRadius - zmin ) / dz );

      // Adjust bin ijk indices if necessary
      imin = LvArray::math::max( imin, 0 );
      imax = LvArray::math::min( imax, nxbins-1 );
      jmin = LvArray::math::max( jmin, 0 );
      jmax = LvArray::math::min( jmax, nybins-1 );
      kmin = LvArray::math::max( kmin, 0 );
      kmax = LvArray::math::min( kmax, nzbins-1 );

      neighborCountsView[pp] = 0; //Think this might have been a bug leading to seg fault when a index was greated than size of neighborCountsView (e.g. active particle indices size)
      localIndex particleIndexOffset2 = 0;
      // Loop all over subRegions
      for( localIndex subRegionIndex2 = 0; subRegionIndex2 < m_numberOfSubRegions; ++subRegionIndex2)
      {
        // Loop over bins
        for( localIndex iBin=imin; iBin<=imax; ++iBin )
        {
          for( localIndex jBin=jmin; jBin<=jmax; ++jBin )
          {
            for( localIndex kBin=kmin; kBin<=kmax; ++kBin )
            {
              localIndex binIndex = subRegionIndex2 * nbins + iBin + jBin * nxbins + kBin * nxbins * nybins;

              for( localIndex bb = 0; bb < binsView.sizeOfArray( binIndex ); ++bb )
              {
                localIndex b = binsView[binIndex][bb];
                real64 xBA[3];
                LvArray::tensorOps::copy< 3 >(xBA, allParticleCentersView[particleIndexOffset2 + b]);
                LvArray::tensorOps::subtract< 3 >(xBA, xA[a]);
                if( LvArray::tensorOps::l2NormSquared< 3 >(xBA) <= neighborRadiusSquared )
                {
                  // Update neighbor count for this particles
                  ++neighborCountsView[pp];
                }
              }
            }
          }
        }
        particleIndexOffset2 += subRegionSizesView[subRegionIndex2];
      }
    } );

    // Must explicitly move array back to host
    subRegionSizesView.move(LvArray::MemorySpace::host);
    neighborCountsView.move(LvArray::MemorySpace::host);

    OrderedVariableToManyParticleRelation & neighborList = subRegionA.neighborList();
    neighborList.freeOnDevice(); // just being careful
    neighborList.resizeFromCapacities( neighborCounts ); // Resize inner arrays from neighborCounts
  } );

  // Perform neighbor search over appropriate bins (populate neighborlist with particle data)
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegionA )
  {
    // Get neighbor list views
    OrderedVariableToManyParticleRelation & neighborList = subRegionA.neighborList();
    // arrayView1d< localIndex const > const & numParticles = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex > neighborRegions = neighborList.m_toParticleRegion.toView();
    ArrayOfArraysView< localIndex > neighborSubRegions = neighborList.m_toParticleSubRegion.toView();
    ArrayOfArraysView< localIndex > neighborIndices = neighborList.m_toParticleIndex.toView();

    // Get 'this' particle's location
    arrayView2d< real64 const > const xA = subRegionA.getParticleCenter();

    // Find neighbors of 'this' particle
    SortedArrayView< localIndex const > const subRegionAActiveParticleIndices = subRegionA.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( subRegionAActiveParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex a = subRegionAActiveParticleIndices[pp];

      // Bin ijk indices bounding a sphere of radius neighborRadius centered at 'this' particle
      localIndex imin = LvArray::math::floor( ( xA[a][0] - neighborRadius - xmin ) / dx );
      localIndex jmin = LvArray::math::floor( ( xA[a][1] - neighborRadius - ymin ) / dy );
      localIndex kmin = LvArray::math::floor( ( xA[a][2] - neighborRadius - zmin ) / dz );
      localIndex imax = LvArray::math::floor( ( xA[a][0] + neighborRadius - xmin ) / dx );
      localIndex jmax = LvArray::math::floor( ( xA[a][1] + neighborRadius - ymin ) / dy );
      localIndex kmax = LvArray::math::floor( ( xA[a][2] + neighborRadius - zmin ) / dz );

      // Adjust bin ijk indices if necessary
      imin = LvArray::math::max( imin, 0 );
      imax = LvArray::math::min( imax, nxbins-1 );
      jmin = LvArray::math::max( jmin, 0 );
      jmax = LvArray::math::min( jmax, nybins-1 );
      kmin = LvArray::math::max( kmin, 0 );
      kmax = LvArray::math::min( kmax, nzbins-1 );

      // Inner subregion loop
      localIndex neighborCount = 0;
      localIndex particleIndexOffset3 = 0;
      for( localIndex subRegionIndex2 = 0; subRegionIndex2 < m_numberOfSubRegions; ++subRegionIndex2 )
      {
        // Loop over bins
        for( localIndex iBin=imin; iBin<=imax; ++iBin )
        {
          for( localIndex jBin=jmin; jBin<=jmax; ++jBin )
          {
            for( localIndex kBin=kmin; kBin<=kmax; ++kBin )
            {
              localIndex binIndex = subRegionIndex2 * nbins + iBin + jBin * nxbins + kBin * nxbins * nybins;
              
              for( localIndex bb = 0; bb < binsView.sizeOfArray( binIndex ); ++bb )
              {
                localIndex b = binsView[binIndex][bb];
                real64 xBA[3];
                LvArray::tensorOps::copy< 3 >( xBA, allParticleCentersView[particleIndexOffset3 + b] );
                LvArray::tensorOps::subtract< 3 >( xBA, xA[a] );
                if( LvArray::tensorOps::l2NormSquared< 3 >(xBA) <= neighborRadiusSquared )
                {
                  neighborRegions[pp][neighborCount] = regionIndicesOfSubRegionsView[subRegionIndex2];
                  neighborSubRegions[pp][neighborCount] = subRegionIndicesInRegionsView[subRegionIndex2];
                  neighborIndices[pp][neighborCount] = b;
                  // neighborRegions[a][neighborCount] = regionIndicesOfSubRegionsView[subRegionIndex2];
                  // neighborSubRegions[a][neighborCount] = subRegionIndicesInRegionsView[subRegionIndex2];
                  // neighborIndices[a][neighborCount] = b;
                  ++neighborCount;
                }
              }
            }
          }
        }
        particleIndexOffset3 += subRegionSizesView[subRegionIndex2];
      }
    } );
    
    // TODO: Do these need to be here if I do not need to print these to console for debugging?
    neighborRegions.move( LvArray::MemorySpace::host );
    neighborSubRegions.move( LvArray::MemorySpace::host );
    neighborIndices.move( LvArray::MemorySpace::host );
  } );

  return( MPI_Wtime() - tStart );
}

void SolidMechanicsMPM::optimizeBinSort( ParticleManager & particleManager )
{
  // Each partition determines its optimal multiplier which results in the minimum time for neighbor list construction
  // The global multiplier is set by a weighted average of each partition's multiplier, with the number of particles
  // on the partition being the weight factor.

  // Start
  int optimalMultiplier = 1;

  // Identify the largest possible multiplier - we can limit this if it's prohibitive to check larger multipliers in large 3D sims
  real64 xL = m_xLocalMaxNoGhost[0] - m_xLocalMinNoGhost[0] + 2 * m_neighborRadius;
  real64 yL = m_xLocalMaxNoGhost[1] - m_xLocalMinNoGhost[1] + 2 * m_neighborRadius;
  real64 zL = m_xLocalMaxNoGhost[2] - m_xLocalMinNoGhost[2] + 2 * m_neighborRadius;
  int maxMultiplier = LvArray::math::max( LvArray::math::ceil( xL / m_neighborRadius ), LvArray::math::max( LvArray::math::ceil( yL / m_neighborRadius ), LvArray::math::ceil( zL / m_neighborRadius ) ));
  maxMultiplier = LvArray::math::max( maxMultiplier, 1 );

  // Identify this partition's optimal multiplier
  real64 minTime = DBL_MAX;
  for( int multiplier = 1; multiplier <= maxMultiplier; ++multiplier )
  {
    m_binSizeMultiplier = multiplier;
    real64 sortingTime = computeNeighborList( particleManager ); // put 0 here for debugging
    if( sortingTime < minTime )
    {
      minTime = sortingTime;
      optimalMultiplier = multiplier;
    }
  }

  // MPI comms
  int globalNumberOfParticles;
  real64 globalWeightedMultiplier;
  int localNumberOfParticles = particleManager.getNumberOfParticles();
  real64 localWeightedMultiplier = optimalMultiplier * localNumberOfParticles;
  MPI_Allreduce( &localWeightedMultiplier,
                 &globalWeightedMultiplier,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 MPI_COMM_GEOS );
  MPI_Allreduce( &localNumberOfParticles,
                 &globalNumberOfParticles,
                 1,
                 MPI_INT,
                 MPI_SUM,
                 MPI_COMM_GEOS );

  // Set bin size multiplier
  m_binSizeMultiplier = LvArray::math::max( (int) LvArray::math::round( globalWeightedMultiplier / globalNumberOfParticles ), 1 );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 SolidMechanicsMPM::kernel( real64 const & r ) // distance from particle to query point.
{
  return kernelDevice( r,
                       m_neighborRadius,
                       m_planeStrain );
}

real64 SolidMechanicsMPM::kernelDevice( real64 const & r, // distance from particle to query point.
                                        real64 const & neighborRadius,
                                        int const & planeStrain ) 
{
  // Compute the value of a particle kernel function at some point.

  const real64 R = neighborRadius; // kernel Radius
  const real64 norm = planeStrain == 1 ? 1.0610329539459689051 / ( R * R ) : 1.1936620731892150183 / ( R * R * R ); // kernel should integrate to unity

  if( r < R )
  {
    return norm * ( 1 - 3.0 * r * r / ( R * R ) + 2.0 * r * r * r / ( R * R * R ) );
  }
  else
  {
    return ( 0.0 );
  }
}

void SolidMechanicsMPM::kernelGradient( arraySlice1d< real64 const > const & x,   // query point
                                        arraySlice1d< real64 const > const & xp,  // particle location
                                        real64 const & r,                         // distance from particle to query point.
                                        real64 (& result)[3] )
{
  real64 temp[3];
  LvArray::tensorOps::copy< 3 >( temp, xp );
  kernelGradientDevice( x,
                        temp,
                        r,
                        m_neighborRadius,
                        m_planeStrain,
                        result );
}

void SolidMechanicsMPM::kernelGradientDevice( arraySlice1d< real64 const > const & x,   // query point
                                              real64 const (& xp)[3], //arraySlice1d< real64 const > const & xp,  // particle location
                                              real64 const & r,                         // distance from particle to query point.
                                              real64 const & neighborRadius,
                                              int const & planeStrain,
                                              real64 (& result) [3] )
{
  // Compute the value of a particle kernel function gradient at some point.

  const real64 R = neighborRadius; // kernel Radius
  const real64 norm = planeStrain ? 1.0610329539459689051 / ( R * R ) : 1.1936620731892150183 / ( R * R * R ); // kernel should integrate to unity

  if( r < R )
  {
    real64 const s = norm * 6.0 * ( r - R ) / ( R * R * R );
    real64 relativePosition[3];
    LvArray::tensorOps::copy< 3 >( relativePosition, x );
    LvArray::tensorOps::subtract< 3 >( relativePosition, xp );
    LvArray::tensorOps::scaledCopy< 3 >( result, relativePosition, s );
  }
  else
  {
    LvArray::tensorOps::fill< 3 >( result, 0.0 );
  }
}

real64 SolidMechanicsMPM::computeKernelField( arraySlice1d< real64 const > const x,  // query point
                                              localIndex const numNeighbors,
                                              arrayView2d< real64 const > const & xp,  // List of neighbor particle locations.
                                              arrayView1d< real64 const > const & Vp,  // List of neighbor particle volumes.
                                              arrayView1d< real64 const > const & fp ) // scalar field values (e.g. damage) at neighbor particles
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.

  // Initialize
  real64 relativePosition[3];
  real64 kernelVal,
         f = 0.0,
         k = 0.0,
         r;

  // Sum
  for( localIndex p = 0; p < numNeighbors; ++p )
  {
    relativePosition[0] = x[0] - xp[p][0];
    relativePosition[1] = x[1] - xp[p][1];
    relativePosition[2] = x[2] - xp[p][2];
    r = sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );
    kernelVal = kernel( r );
    k += Vp[p] * kernelVal;
    f += Vp[p] * fp[p] * kernelVal;
  }

  // Return the normalized kernel field (which eliminates edge effects)
  if( k > 0.0 )
  {
    return ( f / k );
  }
  else
  {
    return ( 0.0 );
  }
}

real64 SolidMechanicsMPM::computeKernelFieldDevice( arraySlice1d< real64 const > const x,  // query point
                                                    real64 const & neighborRadius,
                                                    int const & planeStrain,
                                                    localIndex const & numNeighbors,
                                                    arraySlice1d< localIndex const > const regionIndices,
                                                    arraySlice1d< localIndex const > const subRegionIndices,
                                                    arraySlice1d< localIndex const > const particleIndices,
                                                    ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                                    ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                                    ParticleManager::ParticleView< arrayView1d< real64 const > > particleDamageView,
                                                    ParticleManager::ParticleView< arrayView1d< int const > > particleSurfaceFlagView,
                                                    ParticleManager::ParticleView< arrayView1d< int const > > particleCohesiveZoneFlag )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.

  // Initialize
  real64 relativePosition[3];
  real64 f = 0.0,
         k = 0.0;

  // Sum
  for( localIndex p = 0; p < numNeighbors; ++p )
  {
    localIndex regionIndex = regionIndices[p];
    localIndex subRegionIndex = subRegionIndices[p];
    localIndex particleIndex = particleIndices[p];

    real64 Vp = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
    
    LvArray::tensorOps::copy< 3 >( relativePosition , x );
    LvArray::tensorOps::subtract< 3 >( relativePosition , particlePositionView[regionIndex][subRegionIndex][particleIndex]);
    real64 r = LvArray::tensorOps::l2Norm< 3 >( relativePosition );
    
    real64 fp = particleDamageView[regionIndex][subRegionIndex][particleIndex]; // Damage for dfg
    if( particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex] == 1 ||
        particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex] == 2 ||
        particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex] == 3 ||
        particleCohesiveZoneFlag[regionIndex][subRegionIndex][particleIndex] == 1 )
    {
      fp = 1.0;
    }

    real64 kernelVal = kernelDevice( r, neighborRadius, planeStrain );
    k += Vp * kernelVal;
    f += Vp * fp * kernelVal;
  }

  // Return the normalized kernel field (which eliminates edge effects)
  if( k > 0.0 )
  {
    return ( f / k );
  }
  else
  {
    return ( 0.0 );
  }
}

void SolidMechanicsMPM::computeKernelFieldGradient( arraySlice1d< real64 const > const x, // query point
                                                    localIndex const numNeighbors,
                                                    arrayView2d< real64 const > const & xp,           // List of neighbor particle locations.
                                                    arrayView1d< real64 const > const & Vp,           // List of neighbor particle volumes.
                                                    arrayView1d< real64 const > const & fp,           // scalar field values (e.g. damage) at neighbor particles
                                                    arraySlice1d< real64 > const result )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 kernelVal,
         f = 0.0,
         k = 0.0,
         r;

  // Gradient of the scalar field
  real64 relativePosition[3],
         fGrad[3] = {},
         kGrad[3] = {},
         kernelGradVal[3];

  for( int p = 0; p < numNeighbors; ++p )
  {
    relativePosition[0] = x[0] - xp[p][0];
    relativePosition[1] = x[1] - xp[p][1];
    relativePosition[2] = x[2] - xp[p][2];
    r = sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );

    kernelVal = kernel( r );
    k += Vp[p] * kernelVal;
    f += Vp[p] * fp[p] * kernelVal;

    kernelGradient( x, xp[p], r, kernelGradVal );
    kGrad[0] += kernelGradVal[0] * Vp[p];
    kGrad[1] += kernelGradVal[1] * Vp[p];
    kGrad[2] += kernelGradVal[2] * Vp[p];
    fGrad[0] += kernelGradVal[0] * Vp[p] * fp[p];
    fGrad[1] += kernelGradVal[1] * Vp[p] * fp[p];
    fGrad[2] += kernelGradVal[2] * Vp[p] * fp[p];
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelField = f/k;
    result[0] = fGrad[0] / k - f * kGrad[0] / (k * k);
    result[1] = fGrad[1] / k - f * kGrad[1] / (k * k);
    result[2] = fGrad[2] / k - f * kGrad[2] / (k * k);
  }
  else
  {
    //kernelField = 0.0;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;
  }
}

void SolidMechanicsMPM::computeKernelFieldGradientDevice( arraySlice1d< real64 const > const x, // query point
                                                          real64 const & neighborRadius,
                                                          int const & planeStrain,
                                                          localIndex const & numNeighbors,
                                                          arraySlice1d< localIndex const > const regionIndices,
                                                          arraySlice1d< localIndex const > const subRegionIndices,
                                                          arraySlice1d< localIndex const > const particleIndices,
                                                          ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                                          ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                                          ParticleManager::ParticleView< arrayView1d< real64 const > > particleDamageView,
                                                          ParticleManager::ParticleView< arrayView1d< int const > > particleSurfaceFlagView,
                                                          ParticleManager::ParticleView< arrayView1d< int const > > particleCohesiveZoneFlag,
                                                          arraySlice1d< real64 > const result )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 f = 0.0,
         k = 0.0;

  // Gradient of the scalar field
  real64 fGrad[3] = {},
         kGrad[3] = {};

  for( int n = 0; n < numNeighbors; ++n )
  {
    localIndex const regionIndex = regionIndices[n];
    localIndex const subRegionIndex = subRegionIndices[n];
    localIndex const particleIndex = particleIndices[n];

    real64 const Vp = particleVolumeView[regionIndex][subRegionIndex][particleIndex];

    real64 xp[3];
    LvArray::tensorOps::copy< 3 >( xp, particlePositionView[regionIndex][subRegionIndex][particleIndex] );

    real64 fp = particleDamageView[regionIndex][subRegionIndex][particleIndex]; // Damage for dfg
    if( particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex] == 1 ||
        particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex] == 2 ||
        particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex] == 3 ||
        particleCohesiveZoneFlag[regionIndex][subRegionIndex][particleIndex] == 1 )
    {
      fp = 1.0;
    }

    real64 relativePosition[3];
    LvArray::tensorOps::copy< 3 >( relativePosition, x );
    LvArray::tensorOps::subtract< 3 >( relativePosition, xp );
    real64 const r = LvArray::tensorOps::l2Norm< 3 >( relativePosition ); // sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );

    real64 const kernelVal = kernelDevice( r, neighborRadius, planeStrain );
    k += Vp * kernelVal;
    f += Vp * fp * kernelVal;

    real64 kernelGradVal[3];
    kernelGradientDevice( x, xp, r, neighborRadius, planeStrain, kernelGradVal );
    LvArray::tensorOps::scaledAdd< 3 >( kGrad, kernelGradVal, Vp );
    LvArray::tensorOps::scaledAdd< 3 >( fGrad, kernelGradVal, fp * Vp );
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelField = f/k;
    result[0] = fGrad[0] / k - f * kGrad[0] / (k * k);
    result[1] = fGrad[1] / k - f * kGrad[1] / (k * k);
    result[2] = fGrad[2] / k - f * kGrad[2] / (k * k);
  }
  else
  {
    //kernelField = 0.0;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;
  }
}

void SolidMechanicsMPM::computeKernelFieldDivergenceDevice( arraySlice1d< real64 const > const x, // query point
                                                            real64 const & neighborRadius,
                                                            int const & planeStrain,
                                                            localIndex const & numNeighbors,
                                                            arraySlice1d< localIndex const > const regionIndices,
                                                            arraySlice1d< localIndex const > const subRegionIndices,
                                                            arraySlice1d< localIndex const > const particleIndices,
                                                            ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                                            ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                                            ParticleManager::ParticleView< arrayView2d< real64 const > > particleVectorView,
                                                            real64 & result )
{
 // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 f[3] = {},
         k = 0.0;

  // Gradient of the scalar field
  real64 fGrad[3] = {},
         kGrad[3] = {};

  for( int p = 0; p < numNeighbors; ++p )
  {
    localIndex const regionIndex = regionIndices[p];
    localIndex const subRegionIndex = subRegionIndices[p];
    localIndex const particleIndex = particleIndices[p];

    real64 const Vp = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
    
    real64 xp[3];
    LvArray::tensorOps::copy< 3 >( xp, particlePositionView[regionIndex][subRegionIndex][particleIndex] );

    arraySlice1d< real64 const > const fp = particleVectorView[regionIndex][subRegionIndex][particleIndex];

    real64 relativePosition[3];
    LvArray::tensorOps::copy< 3 >( relativePosition, x );
    LvArray::tensorOps::subtract< 3 >( relativePosition, xp );
    real64 const r = LvArray::tensorOps::l2Norm< 3 >( relativePosition );

    real64 const kernelVal = kernelDevice( r, neighborRadius, planeStrain );
    k += Vp * kernelVal;
    LvArray::tensorOps::scaledAdd< 3 >( f, fp, kernelVal * Vp);

    real64 kernelGradVal[3];
    kernelGradientDevice( x, xp, r, neighborRadius, planeStrain, kernelGradVal );
    LvArray::tensorOps::scaledAdd< 3 >( kGrad, kernelGradVal, Vp);

    fGrad[0] += kernelGradVal[0] * Vp * fp[0];
    fGrad[1] += kernelGradVal[1] * Vp * fp[1];
    fGrad[2] += kernelGradVal[2] * Vp * fp[2];
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelDivergence = f/k;
    result = fGrad[0] / k - f[0] * kGrad[0] / (k * k) + fGrad[1] / k - f[1] * kGrad[1] / (k * k) + fGrad[2] / k - f[2] * kGrad[2] / (k * k);
  }
  else
  {
    //kernelDivergence = 0.0;
    result = 0.0;
  }
}

void SolidMechanicsMPM::computeKernelVectorGradient( arraySlice1d< real64 const > const x,  // query point
                                                     localIndex const numNeighbors,
                                                     arrayView2d< real64 const > const & xp,            // List of neighbor particle locations.
                                                     arrayView1d< real64 const > const & Vp,            // List of neighbor particle volumes.
                                                     arrayView2d< real64 const > const & fp,            // vector field values (e.g. velocity) at neighbor particles
                                                     arraySlice2d< real64 > const result )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 kernelVal,
         f[3],
         k = 0.0,
         r;

  // Gradient of the scalar field
  real64 relativePosition[3],
         fGrad[3][3],
         kGrad[3],
         kernelGradVal[3];

  for( int p = 0; p < numNeighbors; ++p )
  {
    for( int i = 0; i < 3; ++i )
    {
      relativePosition[i] = x[i] - xp[p][i];
    }
    r = sqrt( relativePosition[0] * relativePosition[0] + relativePosition[1] * relativePosition[1] + relativePosition[2] * relativePosition[2] );

    kernelVal = kernel( r );
    k += Vp[p] * kernelVal;
    for( int i = 0; i < 3; ++i )
    {
      f[i] += Vp[p] * fp[p][i] * kernelVal;
    }

    kernelGradient( x, xp[p], r, kernelGradVal );
    for( int i = 0; i < 3; ++i )
    {
      kGrad[i] += kernelGradVal[i] * Vp[p];
      for( int j = 0; j < 3; ++j )
      {
        fGrad[i][j] += fp[p][i] * kernelGradVal[j] * Vp[p];
      }
    }
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelField = f/k;
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        result[i][j] = fGrad[i][j] / k - f[i] * kGrad[j] / (k * k);
      }
    }
  }
  else
  {
    //kernelField = 0.0;
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        result[i][j] = 0.0;
      }
    }
  }
}

void SolidMechanicsMPM::computeKernelVectorGradientDevice( arraySlice1d< real64 const > const x,       // query point
                                                           real64 const & neighborRadius,
                                                           int const & planeStrain,
                                                           localIndex const numNeighbors,
                                                           arraySlice1d< localIndex const > const regionIndices,
                                                           arraySlice1d< localIndex const > const subRegionIndices,
                                                           arraySlice1d< localIndex const > const particleIndices,
                                                           ParticleManager::ParticleView< arrayView1d< real64 const > > particleVolumeView,
                                                           ParticleManager::ParticleView< arrayView2d< real64 const > > particlePositionView,
                                                           ParticleManager::ParticleView< arrayView2d< real64 const > > particleDisplacementView,
                                                           arraySlice2d< real64 > const result )
{
  // Compute the kernel scalar field at a point, for a given list of neighbor particles.
  // The lists xp, fp, and the length np could refer to all the particles in the patch,
  // but generally this function will be evaluated with x equal to some particle center,
  // and xp, fp, will be lists for just the neighbors of the particle.
  // TODO: Modify to also "return" the kernel field value

  // Scalar kernel field values
  real64 f[3] = {},
         k = 0.0;;

  // Gradient of the scalar field
  real64 fGrad[3][3] = {},
         kGrad[3] = {};

  for( int n = 0; n < numNeighbors; ++n )
  {
    localIndex const regionIndex = regionIndices[n];
    localIndex const subRegionIndex = subRegionIndices[n];
    localIndex const particleIndex = particleIndices[n];

    real64 const Vp = particleVolumeView[regionIndex][subRegionIndex][particleIndex];

    real64 xp[3];
    LvArray::tensorOps::copy< 3 >( xp, particlePositionView[regionIndex][subRegionIndex][particleIndex] );
    
    real64 fp[3];
    LvArray::tensorOps::copy< 3 >( fp, particleDisplacementView[regionIndex][subRegionIndex][particleIndex] );

    real64 relativePosition[3];
    LvArray::tensorOps::copy< 3 >( relativePosition, x );
    LvArray::tensorOps::subtract< 3 >( relativePosition, xp );
    real64 r = LvArray::tensorOps::l2Norm< 3 >( relativePosition );

    real64 const kernelVal = kernelDevice( r, neighborRadius, planeStrain );
    k += Vp * kernelVal;
    LvArray::tensorOps::scaledAdd< 3 >( f, fp, Vp * kernelVal );

    real64 kernelGradVal[3];
    kernelGradientDevice( x, xp, r, neighborRadius, planeStrain, kernelGradVal );
    LvArray::tensorOps::scaledAdd< 3 >( kGrad, kernelGradVal, Vp );
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        fGrad[i][j] += fp[i] * kernelGradVal[j] * Vp;
      }
    }
  }

  // Return the normalized kernel field gradient (which eliminates edge effects)
  if( k > 0.0 )
  {
    //kernelField = f/k;
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        result[i][j] = fGrad[i][j] / k - f[i] * kGrad[j] / (k * k);
      }
    }
  }
  else
  {
    //kernelField = 0.0;
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        result[i][j] = 0.0;
      }
    }
  }
}

void SolidMechanicsMPM::computeDamageFieldGradient( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const planeStrain = m_planeStrain;
  real64 const neighborRadius = m_neighborRadius;

  // Get accessors for volume, position, damage, surface flag, and cohesive flag
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > const particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > const particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > const particleDamageAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleDamage" );
  ParticleManager::ParticleViewAccessor< arrayView1d< int const > > const particleSurfaceFlagAccessor = particleManager.constructArrayViewAccessor< int, 1 >( "particleSurfaceFlag" );
  ParticleManager::ParticleViewAccessor< arrayView1d< int const > > const particleCohesiveZoneFlagAccessor = particleManager.constructArrayViewAccessor< int, 1 >( "particleCohesiveZoneFlag" );

  // Get views of accessors
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > const particleVolumeView = particleVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particlePositionView = particlePositionAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > const particleDamageView = particleDamageAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< int const > > const particleSurfaceFlagView = particleSurfaceFlagAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< int const > > const particleCohesiveZoneFlagView = particleCohesiveZoneFlagAccessor.toNestedViewConst();

  // Perform neighbor operations
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();

    // Get const views for array of arrays
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and damage field gradient
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Call kernel field gradient function
      computeKernelFieldGradientDevice( particlePosition[p],
                                        neighborRadius,
                                        planeStrain,
                                        numNeighborsAll[pp],
                                        neighborRegions[pp],
                                        neighborSubRegions[pp],
                                        neighborIndices[pp],
                                        particleVolumeView,
                                        particlePositionView,
                                        particleDamageView,
                                        particleSurfaceFlagView,
                                        particleCohesiveZoneFlagView,
                                        particleDamageGradient[p] );
    } );

    // Required to avoid invalid pointer error of nested arrayViews
    particleVolumeView.move( LvArray::MemorySpace::host );
    particlePositionView.move( LvArray::MemorySpace::host );
    particleDamageView.move( LvArray::MemorySpace::host );
    particleSurfaceFlagView.move( LvArray::MemorySpace::host );
    particleCohesiveZoneFlagView.move( LvArray::MemorySpace::host );

    numNeighborsAll.move( LvArray::MemorySpace::host );
    neighborRegions.move( LvArray::MemorySpace::host );
    neighborSubRegions.move( LvArray::MemorySpace::host );
    neighborIndices.move( LvArray::MemorySpace::host );

    // for(int pp=0; pp < activeParticleIndices.size(); ++pp)
    // {
    //   localIndex const p = activeParticleIndices[pp];

    //   if(cycleNumber > 9190 && cycleNumber < 9290)
    //   {
    //   printf("subRegion: %d, p: %d, neighbors: %d\n", subRegionIndex, p, numNeighborsAll[p]);
    //   }
    //   // for( int n = 0; n < numNeighborsAll[p]; ++n)
    //   // {
    //     // int const regionIndex = neighborRegions[p][n];
    //     // int const subRegionIndex = neighborSubRegions[p][n];
    //     // int const particleIndex = neighborIndices[p][n];
    //   forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( localIndex const ppp )
    //   {
    //     if(cycleNumber > 9190 && cycleNumber < 9290)
    //     {
    //       for (int n = 0; n < numNeighborsAll[p]; ++n)
    //       {
    //         int const regionIndex = neighborRegions[p][n];
    //         int const subRegionIndex = neighborSubRegions[p][n];
    //         int const particleIndex = neighborIndices[p][n];
    //         printf("\t(regionIndex=%d, subRegionIndex=%d, particleIndex=%d\n",
    //                regionIndex,
    //                subRegionIndex,
    //                particleIndex);
    //       }
    //     }
    //   } );

    //   forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( localIndex const ppp )
    //   {
    //     if(cycleNumber > 9190 && cycleNumber < 9290)
    //     {
    //       for (int n = 0; n < numNeighborsAll[p]; ++n)
    //       {
    //         int const regionIndex = neighborRegions[p][n];
    //         int const subRegionIndex = neighborSubRegions[p][n];
    //         int const particleIndex = neighborIndices[p][n];
    //         printf("\tpVolume: %f, pPos: {%f, %f, %f}, pDamage: %f, pSurfFlag: %d, pCZFlag: %d)\n",
    //                 particleVolumeView[regionIndex][subRegionIndex][particleIndex],
    //                 particlePositionView[regionIndex][subRegionIndex][particleIndex][0],
    //                 particlePositionView[regionIndex][subRegionIndex][particleIndex][1],
    //                 particlePositionView[regionIndex][subRegionIndex][particleIndex][2],
    //                 particleDamageView[regionIndex][subRegionIndex][particleIndex],
    //                 particleSurfaceFlagView[regionIndex][subRegionIndex][particleIndex],
    //                 particleCohesiveZoneFlagView[regionIndex][subRegionIndex][particleIndex]);
    //       }
    //     }
    //   } );

    //   forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( localIndex const ppp )
    //   {
    //     // Call kernel field gradient function
    //     computeKernelFieldGradientDevice( particlePosition[p],
    //                                       neighborRadius,
    //                                       planeStrain,
    //                                       numNeighborsAll[p],
    //                                       neighborRegions[p],
    //                                       neighborSubRegions[p],
    //                                       neighborIndices[p],
    //                                       particleVolumeView,
    //                                       particlePositionView,
    //                                       particleDamageView,
    //                                       particleSurfaceFlagView,
    //                                       particleCohesiveZoneFlagView,
    //                                       particleDamageGradient[p] );
    //   } );

    //   // Required to avoid invalid pointer error of nested arrayViews
    //   particleVolumeView.move( LvArray::MemorySpace::host );
    //   particlePositionView.move( LvArray::MemorySpace::host );
    //   particleDamageView.move( LvArray::MemorySpace::host );
    //   particleSurfaceFlagView.move( LvArray::MemorySpace::host );
    //   particleCohesiveZoneFlagView.move( LvArray::MemorySpace::host );

    //   numNeighborsAll.move( LvArray::MemorySpace::host );
    //   neighborRegions.move( LvArray::MemorySpace::host );
    //   neighborSubRegions.move( LvArray::MemorySpace::host );
    //   neighborIndices.move( LvArray::MemorySpace::host );
    // // }
    // if(cycleNumber > 9190 && cycleNumber < 9290 )
    // {
    // printf("Exited - p: %d\n", p);
    // }
  // }
  ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::estimateParticleSurfacePosition( ParticleManager & particleManager )
{
  real64 const neighborRadius = m_neighborRadius;
  int const planeStrain = m_planeStrain;
  
  // Get accessors for volume, position, damage, surface flag, and cohesive flag
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > const particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > const particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > const particleDamageAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleDamage" );
  ParticleManager::ParticleViewAccessor< arrayView1d< int const > > const particleSurfaceFlagAccessor = particleManager.constructArrayViewAccessor< int, 1 >( "particleSurfaceFlag" );
  ParticleManager::ParticleViewAccessor< arrayView1d< int const > > const particleCohesiveZoneFlagAccessor = particleManager.constructArrayViewAccessor< int, 1 >( "particleCohesiveZoneFlag" );

  // Get views of accessors
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > const particleVolumeView = particleVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particlePositionView = particlePositionAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > const particleDamageView = particleDamageAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< int const > > const particleSurfaceFlagView = particleSurfaceFlagAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< int const > > const particleCohesiveZoneFlagView = particleCohesiveZoneFlagAccessor.toNestedViewConst();

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();

    // Get const views for array of arrays
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 > const particleEstimatedSurfacePosition = subRegion.getField< fields::mpm::particleEstimatedSurfacePosition >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      real64 particleNonlocalDamage = computeKernelFieldDevice( particlePosition[p], // query point
                                                                neighborRadius,
                                                                planeStrain,
                                                                numNeighborsAll[pp],
                                                                neighborRegions[pp],
                                                                neighborSubRegions[pp],
                                                                neighborIndices[pp],
                                                                particleVolumeView,
                                                                particlePositionView,
                                                                particleDamageView,
                                                                particleSurfaceFlagView,
                                                                particleCohesiveZoneFlagView );
                    
      real64 particleDistanceToSurface = inverseKernel( particleNonlocalDamage );

      real64 temp[3];
      LvArray::tensorOps::copy< 3 >( temp, particleDamageGradient[p] );
      LvArray::tensorOps::normalize< 3 >( temp );
      LvArray::tensorOps::scaledCopy< 3 >( particleEstimatedSurfacePosition[p], temp, particleDistanceToSurface );
    } );
  } );
}

// Compute the value r from the w by inverting kernel w(r)
real64 SolidMechanicsMPM::inverseKernel( real64 const & d ) // value of kernel function
{
	// This is a fast newton solver to invert the cubic kernel function, and return the
	// distance from the value of the kernel field.
    real64 x = 0.5; // initial guess
    real64 x_prev = x - 1; // previous value of x
    real64 eps = 1e-3;

    while (fabs(x - x_prev) > eps) {
        x_prev = x;
        x -= (1 + 2 * x * x * x - 3 * x * x - d) / (6 * x * x - 6 * x);
    }

    // returns distance from crack tip:
    return x*m_neighborRadius;
}

void SolidMechanicsMPM::updateSurfaceFlagOverload( ParticleManager & particleManager )
{
  // Surface flags are overloaded so we can visualize surfaces and damage features simultaneously
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      if( particleSurfaceFlag[p] < 2 )
      {
        if( particleDamage[p] > 0.0 ) // Activate damage field if any particles in domain have damage.
        {
          particleSurfaceFlag[p] = 1;
        }
        else
        {
          particleSurfaceFlag[p] = 0;
        }
      }
    } );
  } );
}

void SolidMechanicsMPM::projectDamageFieldGradientToGrid( DomainPartition & domain,
                                                          ParticleManager & particleManager,
                                                          NodeManager & nodeManager,
                                                          MeshLevel & mesh )
{
  GEOS_MARK_FUNCTION;

  int const numDims = m_numDims;

  // Grid nodes gain the damage field gradient of the particle mapping to them with the largest damage field gradient

  // Get grid fields
  arrayView2d< real64 > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get particle fields
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    // arrayView1d< int const > const particleSurfaceFlag = subRegion.getField< fields::mpm::particleSurfaceFlag >();
    arrayView2d< real64 const> const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get nodes this particle maps to
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp ) // Parallelize with atomics/reduction
    {
      localIndex const p = activeParticleIndices[pp];

      // Currently assume surface flags 2 and 3 and above have surface normals defined from PFW
      // This seems more efficient, e.g. avoid computing l2Norm for every particle, but may need to be changed if additional surface flag types are added
      real64 damageGradient[3];
      if( LvArray::tensorOps::l2NormSquared< 3 >(particleSurfaceNormal[p]) > 1e-12 ) // particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ) // I think this needs to be check
      {
        LvArray::tensorOps::copy< 3 >( damageGradient, particleSurfaceNormal[p] );
        LvArray::tensorOps::scale< 3 >( damageGradient, ( 2 + particleID[p] ) / m_neighborRadius );
      }
      else
      {
        LvArray::tensorOps::copy< 3 >( damageGradient, particleDamageGradient[p] );
      }

      // Map to grid
      for( localIndex const & g: mappedNodes[pp] )
      {
        if( LvArray::tensorOps::l2NormSquared< 3 >( damageGradient ) > LvArray::tensorOps::l2NormSquared< 3 >( gridDamageGradient[g] ) )
        {
          for( localIndex i=0; i < numDims; ++i )
          {
            gridDamageGradient[g][i] = damageGradient[i];
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop

  // Sync damage gradient field
  syncGridFields( { viewKeyStruct::gridDamageGradientString() }, domain, nodeManager, mesh, MPI_MAX );
}

void SolidMechanicsMPM::computeParticleFieldMappings( ParticleManager & particleManager,
                                                      NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
 
  // Grid fields
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< integer > const mappedFields = m_mappedFields[subRegionIndex];
    
    // Particle fields
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( integer g = 0; g < numberOfVerticesPerParticle * 8; ++g )
      {
        localIndex const mappedNode = mappedNodes[pp][g];

        mappedFields[pp][g] = partitionField( numContactGroups,
                                              damageFieldPartitioning,
                                              particleGroup[p],
                                              particleDamageGradient[p],
                                              particleSurfaceNormal[p],
                                              gridDamageGradient[mappedNode] );
      }
    } );

    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::updateDeformationGradient( real64 dt,
                                                   ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const FSubcycles = m_FSubcycles;
  int const planeStrain = m_planeStrain;
  int const exactJIntegration = m_exactJIntegration;

  real64 dtSub = dt / m_FSubcycles;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();

    // Update F
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // UPDATE DEFORMATION GRADIENT L = Fdot * Finv -> F_new = F_old + Fdot*dt = F_old + L*F_old*dt
      real64 oldDeformationGradient[3][3] = {};
      LvArray::tensorOps::copy< 3, 3 >( oldDeformationGradient, particleDeformationGradient[p] );

      // Integrate F with L
      real64 previousDeformationGradient[3][3] = {};
      for( int iter = 0 ; iter < FSubcycles ; ++iter )
      {
        LvArray::tensorOps::copy< 3, 3 >( previousDeformationGradient, particleDeformationGradient[p] );

        real64 temp[3][3] = {};
        LvArray::tensorOps::copy< 3, 3 >( temp, previousDeformationGradient );

        real64 tempFDot[3][3] = {};
        LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( tempFDot, particleVelocityGradient[p], temp );
        LvArray::tensorOps::scale< 3, 3 >( tempFDot, dtSub);
        LvArray::tensorOps::add< 3, 3 >( particleDeformationGradient[p], tempFDot );
      }

      if ( exactJIntegration == 1 )
      { // This will improve accuracy with nearly incompressible materials
        // Exact integration of J with tr(L)
        real64 Jold = LvArray::tensorOps::determinant< 3 >( oldDeformationGradient );
        real64 Jnew = Jold * exp( LvArray::tensorOps::trace< 3 >( particleVelocityGradient[p] ) * dt );

        // Modify F to reduce volumetric integration error
        real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
        real64 scale = 1.0;
        if( detF > 0.0 && Jnew >= 0.0 )
        {
          real64 power = planeStrain ? 0.5 : 1.0/3.0;
          scale = pow( Jnew / detF, power );
        }
        LvArray::tensorOps::scale< 3, 3 >( particleDeformationGradient[p], scale );
      }

      real64 temp[3][3] = {};
      LvArray::tensorOps::copy< 3, 3 >( temp, particleDeformationGradient[p] );

      real64 minusOldDeformationGradient[3][3] = {};
      LvArray::tensorOps::scaledCopy< 3, 3 >( minusOldDeformationGradient, oldDeformationGradient, -1.0);
      LvArray::tensorOps::add< 3, 3 >( temp, minusOldDeformationGradient );
      LvArray::tensorOps::scale< 3, 3 >( temp, 1 / dt );
      LvArray::tensorOps::copy< 3, 3 >( particleFDot[p], temp );
    } );

    particleDeformationGradient.move( LvArray::MemorySpace::host );
    particleFDot.move( LvArray::MemorySpace::host );
  } );
}

// Check whether it is faster to do one forAll with if checks of constitutive members inside rather than a bunch of individual forAlls for each field
void SolidMechanicsMPM::updateConstitutiveModelDependencies( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get needed particle fields
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // arrayView2d< real64 const > const particleMaterialDirection = subRegion.getParticleMaterialDirection();
    // arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    // arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    // arrayView1d< real64 const > const particleTemperature = subRegion.getParticleTemperature();
    // arrayView1d< real64 const > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();

    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );

    // arrayView1d< real64 > constitutiveLengthScale;
    // arrayView2d< real64 > constitutiveMaterialDirection;
    // arrayView3d< real64 > constitutiveDeformationGradient;
    // arrayView3d< real64 > constitutiveVelocityGradient;
    // arrayView1d< real64 > constitutiveTemperature;
    // arrayView1d< real64 > constitutiveVolume;
    // arrayView2d< real64 > constitutiveDensity;
    // arrayView2d< real64 > constitutiveJacobian;

    // int hasLengthScale = constitutiveModel.hasWrapper( "lengthScale" );
    // int hasMaterialDirection = constitutiveModel.hasWrapper( "materialDirection" );
    // int hasDeformationGradient = constitutiveModel.hasWrapper( "deformationGradient" );
    // int hasVelocityGradient = constitutiveModel.hasWrapper( "velocityGradient" );
    // int hasTemperature = constitutiveModel.hasWrapper( "temperature" );
    // int hasVolume = constitutiveModel.hasWrapper( "volume" );
    // int hasDensity = constitutiveModel.hasWrapper( "density" );
    // int hasJacobian = constitutiveModel.hasWrapper( "jacobian" );

    // if( hasLengthScale )
    // {
    //   constitutiveLengthScale = constitutiveModel.getReference< array1d< real64 > >( "lengthScale" );
    // }

    // if( hasMaterialDirection )
    // {
    //   constitutiveMaterialDirection = constitutiveModel.getReference< array2d< real64 > >( "materialDirection" );
    // }

    // if( hasDeformationGradient )
    // {
    //   constitutiveDeformationGradient = constitutiveModel.getReference< array3d< real64 > >( "deformationGradient" );
    // }

    // if( hasVelocityGradient )
    // {
    //   constitutiveVelocityGradient = constitutiveModel.getReference< array3d< real64 > >( "velocityGradient" );
    // }

    // if( hasVolume )
    // {
    //   constitutiveVolume = constitutiveModel.getReference< array1d< real64 > >( "volume" );
    // }

    // if( hasDensity )
    // {
    //   constitutiveDensity = constitutiveModel.getReference< array2d< real64 > >( "density" );
    // }

    // if( hasJacobian )
    // {
    //   constitutiveJacobian = constitutiveModel.getReference< array2d< real64 > >( "jacobian" );
    // }

    // forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    // {
    //   localIndex const p = activeParticleIndices[pp];

    //   if( hasLengthScale )
    //   {
    //     constitutiveLengthScale[p] = pow( particleVolume[p], 1.0 / 3.0 );
    //   }

    //   if( hasMaterialDirection )
    //   {
    //     LvArray::tensorOps::copy< 3 >( constitutiveMaterialDirection[p], particleMaterialDirection[p] );
    //   }

    //   if( hasDeformationGradient )
    //   {
    //     LvArray::tensorOps::copy< 3, 3 >( constitutiveDeformationGradient[p], particleDeformationGradient[p] );
    //   }

    //   if( hasVelocityGradient )
    //   {
    //     LvArray::tensorOps::copy< 3, 3 >( constitutiveVelocityGradient[p], particleVelocityGradient[p] );
    //   }

    //   if( hasTemperature )
    //   {
    //     constitutiveTemperature[p] = particleTemperature[p];
    //   }

    //   if( hasVolume )
    //   {
    //     constitutiveVolume[p] = particleVolume[p];
    //   }

    //   if( hasDensity )
    //   {
    //     constitutiveDensity[p][0] = particleDensity[p];
    //   }

    //   if( hasJacobian )
    //   {
    //     constitutiveJacobian[p][0] = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
    //   }
    // } );

    // Pass whatever data the constitutive models may need
    if( constitutiveModel.hasWrapper( "lengthScale" ) ) // Fragile code because someone could change this key without our knowledge. TODO: Make an
                                                        // integrated test that checks this
    {
      arrayView1d< real64 > const lengthScale = constitutiveModel.getReference< array1d< real64 > >( "lengthScale" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        lengthScale[p] = pow( particleVolume[p], 1.0 / 3.0 );
      } );
    }

    if(  constitutiveModel.hasWrapper( "materialDirection" ) )
    {
      // CC: Todo add check for fiber vs plane update to material direction
      arrayView2d< real64 const > const particleMaterialDirection = subRegion.getParticleMaterialDirection();
      arrayView2d< real64 > const constitutiveMaterialDirection = constitutiveModel.getReference< array2d< real64 > >( "materialDirection" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 3 >( constitutiveMaterialDirection[p], particleMaterialDirection[p] );
      } );
    }

    if(  constitutiveModel.hasWrapper( "deformationGradient" ) )
    {
      arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      arrayView3d< real64 > const constitutiveDeformationGradient = constitutiveModel.getReference< array3d< real64 > >( "deformationGradient" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 3, 3 >(constitutiveDeformationGradient[p], particleDeformationGradient[p]);
      } );
    }

    if(  constitutiveModel.hasWrapper( "velocityGradient" ) )
    {
      arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
      arrayView3d< real64 > const constitutiveVelocityGradient = constitutiveModel.getReference< array3d< real64 > >( "velocityGradient" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 3, 3 >(constitutiveVelocityGradient[p], particleVelocityGradient[p]);
      } );
    }

    if(  constitutiveModel.hasWrapper( "temperature" ) )
    {
      arrayView1d< real64 const > const particleTemperature = subRegion.getParticleTemperature();
      arrayView1d< real64 > const constitutiveTemperature = constitutiveModel.getReference< array1d< real64 > >( "temperature" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveTemperature[p] = particleTemperature[p];
      } );
    }

    if(  constitutiveModel.hasWrapper( "volume" ) )
    {
      arrayView1d< real64 > const constitutiveVolume = constitutiveModel.getReference< array1d< real64 > >( "volume" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveVolume[p] = particleVolume[p];
      } );
    }

    if(  constitutiveModel.hasWrapper( "density" ) )
    {
      arrayView1d< real64 const > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
      arrayView2d< real64 > const constitutiveDensity = constitutiveModel.getReference< array2d< real64 > >( "density" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveDensity[p][0] = particleDensity[p];
      } );
    }

    if(  constitutiveModel.hasWrapper( "jacobian" ) )
    {
      arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      arrayView2d< real64 > const constitutiveJacobian = constitutiveModel.getReference< array2d< real64 > >( "jacobian" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        constitutiveJacobian[p][0] = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );
      } );
    }
  } );
}

void SolidMechanicsMPM::updateStress( real64 dt,
                                      ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & solid = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );

    // Get particle kinematic fields that are fed into constitutive model
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 const > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    int hyperelasticUpdate = 0;
    if ( solid.getCatalogName() == "HyperelasticMMS" || solid.getCatalogName() == "Hyperelastic" || solid.getCatalogName() == "Chiumenti" )
    {
      hyperelasticUpdate = 1;
    }

    // Call constitutive model
    ConstitutivePassThruMPM< ContinuumBase >::execute( solid, [&] ( auto & castedSolid )
    {
      using SolidType = TYPEOFREF( castedSolid );
      typename SolidType::KernelWrapper constitutiveModelWrapper = castedSolid.createKernelUpdates();
      solidMechanicsMPMKernels::StateUpdateKernel::launch< parallelDevicePolicy<> >( subRegion.activeParticleIndices(),
                                                                                     constitutiveModelWrapper,
                                                                                     dt,
                                                                                     hyperelasticUpdate,
                                                                                     particleDeformationGradient,
                                                                                     particleFDot,
                                                                                     particleVelocityGradient,
                                                                                     particleStress );
    } );
  } );
}

void SolidMechanicsMPM::particleKinematicUpdate( const real64 dt,
                                                 ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const numDims = m_numDims;
  real64 const minParticleJacobian = m_minParticleJacobian;
  real64 const maxParticleJacobian = m_maxParticleJacobian;
  real64 const maxParticleVelocitySquared = m_maxParticleVelocitySquared;
  real64 const numeric_max = std::numeric_limits< real64 >::max();
  int const useReferenceVectorsForParticleUpdate = m_useReferenceVectorsForParticleUpdate;

  RAJA::ReduceSum< parallelDeviceReduce, int > numParticlesIllConditionedJacobian( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, int > numParticlesVelocityOverflowed( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, int > numParticlesOverMaxVelocity( 0 );

  // Update particle volume and density
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView1d< real64 > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
    arrayView2d< real64 > const particleMaterialDirection = subRegion.getParticleMaterialDirection();
    arrayView2d< real64 > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView2d< real64 > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

    arrayView1d< int > const particleDeleteFlag = subRegion.getField< fields::mpm::particleDeleteFlag >();
    arrayView1d< real64 const > const particleReferenceVolume = subRegion.getField< fields::mpm::particleReferenceVolume >();
    arrayView1d< real64 > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 const > const particleReferenceRVectors = subRegion.getField< fields::mpm::particleReferenceRVectors >();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 const > const particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView2d< real64 const > const particleReferenceMaterialDirection = subRegion.getField< fields::mpm::particleReferenceMaterialDirection >();
    arrayView2d< real64 const > const particleReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
    arrayView2d< real64 const > const particleReferenceSurfacePosition = subRegion.getField< fields::mpm::particleReferenceSurfacePosition >();
    arrayView2d< real64 const > const particleReferenceSurfaceTraction = subRegion.getField< fields::mpm::particleReferenceSurfaceTraction >();

    // Get constitutive model reference to check if material represents a fiber
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );
    bool isFiber = constitutiveModel.hasWrapper( "isFiber" ); // dumby variable whose value doesn't matter only that it is defined ( possibly a better way to do this )

    // Update volume and r-vectors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      real64 detF = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p] );

      bool flaggedForDeletion = false;
      if( detF <= minParticleJacobian || detF >= maxParticleJacobian )
      {
        numParticlesIllConditionedJacobian += 1;
        flaggedForDeletion = true;
      }

      // With body forces and surface tractions, particles that become detached can accelerate sufficiently to overflow the velocity squared
      // Here we detect if particle velocities will overflow when squared and flag them for deletion to avoid erroring out
      real64 particleSpeedSquared = 0;
      for( int d = 0; d < numDims; ++d  )
      {
        real64 addSqr = particleVelocity[p][d] * particleVelocity[p][d];
        if( particleSpeedSquared > numeric_max - addSqr )
        {
          numParticlesVelocityOverflowed += 1;
          flaggedForDeletion = true;
          break;
        }
        particleSpeedSquared += addSqr;
      }

      if( !flaggedForDeletion && particleSpeedSquared > maxParticleVelocitySquared )
      {
        numParticlesOverMaxVelocity += 1;
        flaggedForDeletion = true;
      }

      if( !flaggedForDeletion )
      {
        particleVolume[p] = particleReferenceVolume[p] * detF;
        particleDensity[p] = particleMass[p] / particleVolume[p];

        if ( useReferenceVectorsForParticleUpdate == 0 )
        {
          real64 deformationGradient[3][3];
          LvArray::tensorOps::copy< 3, 3 >( deformationGradient, particleDeformationGradient[p] );

          // Update particle surface normals
          real64 deformationGradientCofactor[3][3];
          LvArray::tensorOps::cofactor< 3 >( deformationGradientCofactor, deformationGradient );

          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleSurfaceNormal[p], deformationGradientCofactor, particleReferenceSurfaceNormal[p] );
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleSurfaceTraction[p], deformationGradientCofactor, particleReferenceSurfaceTraction[p] );

          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleSurfacePosition[p], deformationGradient, particleReferenceSurfacePosition[p] );

          // Update the material direction
          real64 materialDirection[3];
          if( isFiber )
          {
            LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( materialDirection, deformationGradient, particleReferenceMaterialDirection[p] );
          }
          else
          {
            LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( materialDirection, deformationGradientCofactor, particleReferenceMaterialDirection[p] );
          }
          LvArray::tensorOps::copy< 3 >(particleMaterialDirection[p], materialDirection);
        }
        else 
        {
          // Compute derivative of deformation gradient cofactor
          //
          real64 invF[3][3];
          LvArray::tensorOps::invert< 3 >( invF, particleDeformationGradient[p] );

          real64 invFT[3][3];
          LvArray::tensorOps::transpose< 3, 3 >( invFT, invF );

          real64 G[3][3];
          LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( G, invF, particleFDot[p] );

          real64 trG = LvArray::tensorOps::trace< 3 >( G );

          real64 GT[3][3];
          LvArray::tensorOps::transpose< 3, 3 >( GT, G );
          // LvArray::tensorOps::scale< 3, 3 >( GT, -1.0 );

          real64 dcofFdt[3][3];
          LvArray::tensorOps::addIdentity< 3 >( dcofFdt , trG );
          LvArray::tensorOps::scale< 3, 3 >( GT, -1.0 );
          LvArray::tensorOps::add< 3, 3 >( dcofFdt, GT );

          real64 out[3][3];
          LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( out, invFT, dcofFdt );
          LvArray::tensorOps::scale< 3, 3 >( out, detF );

          real64 dcofF[3][3];
          LvArray::tensorOps::scaledCopy< 3, 3 >( dcofF, out, dt );

          real64 dF[3][3];
          LvArray::tensorOps::scaledCopy< 3, 3 >( dF, particleDeformationGradient[p], dt );

          // Update quantities
          real64 temp[3]= { 0.0 };

          // Particle surface normal         
          LvArray::tensorOps::copy< 3 >( temp, particleSurfaceNormal[p] );
          LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( particleSurfaceNormal[p], dcofF, temp );

          // Particle surface traction
          LvArray::tensorOps::copy< 3 >( temp, particleSurfaceTraction[p] );
          LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( particleSurfaceTraction[p], dcofF, temp );

          // Particle surface position
          LvArray::tensorOps::copy< 3 >( temp, particleSurfacePosition[p] );
          LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( particleSurfacePosition[p], dF, temp );

          // Particle material direction
          real64 materialDirection[3];
          if( isFiber )
          {
            LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( materialDirection, dF, particleMaterialDirection[p] );
          }
          else
          {
            LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( materialDirection, dcofF, particleMaterialDirection[p] );
          }
          LvArray::tensorOps::copy< 3 >(particleMaterialDirection[p], materialDirection);
        }

        // Renormalize particle surface normals
        real64 norm = LvArray::tensorOps::l2Norm< 3 >( particleSurfaceNormal[p] );
        if( norm > 1e-16)
        {
          LvArray::tensorOps::scale< 3 >( particleSurfaceNormal[p], 1/norm );
        }
      }
      else
      {
        particleDeleteFlag[p] = 1;
        particleVolume[p] = particleReferenceVolume[p];
        particleDensity[p] = particleMass[p] / particleReferenceVolume[p];
        LvArray::tensorOps::copy< 3, 3 >( particleRVectors[p], particleReferenceRVectors[p] );
      }
    } );

    // Remove particles to be deleted from active indicies by reconstructing them
    subRegion.setActiveParticleIndices();
  } );

  int numParticlesIllConditionedJacobianGlobal = MpiWrapper::sum( numParticlesIllConditionedJacobian.get() );
  int numParticlesVelocityOverflowedGlobal = MpiWrapper::sum( numParticlesVelocityOverflowed.get() );
  int numParticlesOverMaxVelocityGlobal = MpiWrapper::sum( numParticlesOverMaxVelocity.get() );

  // MpiWrapper::allReduce< int >( &numParticlesIllConditionedJacobian,
  //                               &numParticlesIllConditionedJacobianGlobal,
  //                               1,
  //                               MPI_SUM,
  //                               MPI_COMM_GEOS );

  // MpiWrapper::allReduce< int >( &numParticlesVelocityOverflowed,
  //                               &numParticlesVelocityOverflowedGlobal,
  //                               1,
  //                               MPI_SUM,
  //                               MPI_COMM_GEOS );

  // MpiWrapper::allReduce< int >( &numParticlesOverMaxVelocity,
  //                               &numParticlesOverMaxVelocityGlobal,
  //                               1,
  //                               MPI_SUM,
  //                               MPI_COMM_GEOS );

  GEOS_LOG_RANK_0_IF( numParticlesIllConditionedJacobianGlobal > 0, "Flagged " << numParticlesIllConditionedJacobianGlobal  << " particles with unreasonable Jacobian (J<" << m_minParticleJacobian << " or J>" << m_maxParticleJacobian << ") for deletion!" );
  GEOS_LOG_RANK_0_IF( numParticlesVelocityOverflowedGlobal > 0, "Flagged " << numParticlesVelocityOverflowedGlobal << " particles velocity squared overflow for deletion!" );
  GEOS_LOG_RANK_0_IF( numParticlesOverMaxVelocityGlobal > 0, "Flagged " << numParticlesOverMaxVelocityGlobal << " particles with unreasonable velocity (v " << m_maxParticleVelocity << ") for deletion!" );

  // Compute particles R vectors
  computeRVectors( particleManager );
}

void SolidMechanicsMPM::computeAndWriteBoxAverage( const real64 dt,
                                                   const real64 time_n,
                                                   ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const computeInternalEnergyAndTemperature = m_computeInternalEnergyAndTemperature;

  // real64 boxPlasticStrain[6];
  // real64 boxStress[6]; // we sum stress * volume in particles, additive sync, then divide by box volume.

  // real64 boxMass = 0.0; // we sum particle mass, additive sync, then divide by box volume
  // real64 boxParticleReferenceVolume = 0.0;
  // real64 boxDamage = 0.0; // we sum damage * reference volume, additive sync, then divide by total reference volume in box
  // real64 boxInternalEnergy = 0.0;
  // real64 boxKineticEnergy = 0.0;
  // real64 boxMatVolume = 0.0; // Sum volume of all particles

  RAJA::MultiReduceSum< RAJA::seq_multi_reduce, real64 > boxPlasticStrain(6, 0.0);
  RAJA::MultiReduceSum< RAJA::seq_multi_reduce, real64 > boxStress(6, 0.0);
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxMass( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxParticleReferenceVolume( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxDamage( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxInternalEnergy( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxKineticEnergy( 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxMatVolume( 0.0 );

  // real64 boxInternalForce[3];
  // real64 boxCohesiveForce[3];
  // real64 boxContactForce[3];

  real64 boxAverageMin[3];
  LvArray::tensorOps::copy< 3 >( boxAverageMin, m_boxAverageMin );
  real64 boxAverageMax[3];
  LvArray::tensorOps::copy< 3 >( boxAverageMax, m_boxAverageMax );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 const > const particleReferenceVolume = subRegion.getField< fields::mpm::particleReferenceVolume >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView1d< real64 > const particleKineticEnergy = subRegion.getField< fields::mpm::particleKineticEnergy >();

    arrayView1d< real64 const > particleInternalEnergy;
    if( computeInternalEnergyAndTemperature == 1 )
    {
      particleInternalEnergy = subRegion.getField< fields::mpm::particleInternalEnergy >();
    }

    // Accumulate values
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    // Need to convert to parallelDevicePolicy and incorporate RAJA reductions
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &boxMass, &boxParticleReferenceVolume, &boxStress, &boxPlasticStrain, &boxDamage, &boxMatVolume, &boxInternalEnergy, & boxKineticEnergy] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];
      real64 x = particlePosition[p][0];
      real64 y = particlePosition[p][1];
      real64 z = particlePosition[p][2];
      if( x > boxAverageMin[0] && x < boxAverageMax[0] &&
          y > boxAverageMin[1] && y < boxAverageMax[1] &&
          z > boxAverageMin[2] && z < boxAverageMax[2] )
      {
        boxMass += particleMass[p];
        boxMatVolume += particleVolume[p];
        boxParticleReferenceVolume += particleReferenceVolume[p];
        boxKineticEnergy += particleKineticEnergy[p] * particleVolume[p];
        
        if( computeInternalEnergyAndTemperature == 1 )
        {
          boxInternalEnergy += particleInternalEnergy[p] * particleVolume[p];
        }

        for( int i=0; i<6; ++i )
        {
          boxStress[i] += particleStress[p][i] * particleVolume[p]; // volume weighted average, will normalize later.
          boxPlasticStrain[i] += particlePlasticStrain[p][i] * particleVolume[p];
        }
        boxDamage += particleDamage[p] * particleReferenceVolume[p]; // reference volume weighted average, will normalize later.
      }
    } );
  } );

  // Additive sync: sxx, syy, szz, sxy, syz, sxz, mass, particle volume, damage
  // Check the voigt indexing of stress
  real64 boxSums[18];
  boxSums[0] = boxStress[0].get();       // sig_xx * volume
  boxSums[1] = boxStress[1].get();       // sig_yy * volume
  boxSums[2] = boxStress[2].get();       // sig_zz * volume
  boxSums[3] = boxStress[3].get();       // sig_yz * volume
  boxSums[4] = boxStress[4].get();       // sig_xz * volume
  boxSums[5] = boxStress[5].get();       // sig_xy * volume
  boxSums[6] = boxMass.get();            // total mass in box
  boxSums[7] = boxParticleReferenceVolume.get() > 0.0 ? boxParticleReferenceVolume.get() : 1.0;  // total particle reference volume in box; prevent div0 error
  boxSums[8] = boxDamage.get();          // damage * volume
  boxSums[9] = boxInternalEnergy.get();          // internal energy * volume
  boxSums[10] = boxKineticEnergy.get(); // kinetic energy * volume
  boxSums[11] = boxPlasticStrain[0].get(); // plasticStrain_xx
  boxSums[12] = boxPlasticStrain[1].get(); // plasticStrain_yy
  boxSums[13] = boxPlasticStrain[2].get(); // plasticStrain_zz
  boxSums[14] = boxPlasticStrain[3].get(); // plasticStrain_yz
  boxSums[15] = boxPlasticStrain[4].get(); // plasticStrain_xz
  boxSums[16] = boxPlasticStrain[5].get(); // plasticStrain_xy
  boxSums[17] = boxMatVolume.get();

  // Do an MPI sync to total these values and write from proc0 to a file.  Also compute global F
  // so file is directly plottable in excel as CSV or something.
  for( localIndex i = 0; i < 18; ++i )
  {
    real64 localSum = boxSums[i];
    real64 globalSum;
    MPI_Allreduce( &localSum,
                   &globalSum,
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   MPI_COMM_GEOS );
    boxSums[i] = globalSum;
  }

  int rank;
  MPI_Comm_rank( MPI_COMM_GEOS, &rank );
  if( rank == 0 )
  {
    // Calculate the box volume
    real64 boxVolume = m_domainExtent[0] * m_domainExtent[1] * m_domainExtent[2];

    // Write to file
    std::ofstream file;
    file.open( "boxAverageHistory.csv", std::ios::out | std::ios::app );
    if( file.fail() )
    {
      throw std::ios_base::failure( std::strerror( errno ) );
    }
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
    // time | sig_xx | sig_yy | sig_zz | sig_xy | sig_yz | sig_zx | density | damage | internal energy | kinetic energy | epxx | epyy | epzz | epyz | epxz | epxy | total particle volume
    file << time_n + dt
         << ","
         << boxSums[0] / boxVolume
         << ","
         << boxSums[1] / boxVolume
         << ","
         << boxSums[2] / boxVolume
         << ","
         << boxSums[3] / boxVolume
         << ","
         << boxSums[4] / boxVolume
         << ","
         << boxSums[5] / boxVolume
         << ","
         << boxSums[6] / boxVolume
         << ","
         << boxSums[8] / boxSums[7] // We normalize by total particle reference volume because this should equal one if all the material is
                                    // damaged
         << ", "
         << boxSums[9] / boxVolume
         << ", "
         << boxSums[10] / boxVolume
         << ", "
         << boxSums[11] / boxVolume
         << ", "
         << boxSums[12] / boxVolume
         << ", "
         << boxSums[13] / boxVolume
         << ", "
         << boxSums[14] / boxVolume
         << ", "
         << boxSums[15] / boxVolume
         << ", "
         << boxSums[16] / boxVolume
         << ", "
         << boxSums[17] << std::endl;
    file.close();
  }
}

void SolidMechanicsMPM::writeParticleData( const real64 GEOS_UNUSED_PARAM( time_n ),
                                           ParticleManager & GEOS_UNUSED_PARAM( particleManager ) )
{
  GEOS_LOG_RANK_0( "Particle data writing to file is not currently implemented!" );
  // // Each particle writes to separate file to avoid race condition (only one rank should possess the particle as a master)
  // particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  // {
  //   // Get fields
  //   arrayView1d< real64 const > const particleID = subRegion.getParticleID();
  //   arrayView1d< real64 const > const particlePosition = subRegion.getParticleCenter();
  //   arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();

  //   SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
  //   for( int pp = 0; pp < activeParticleIndices.size(); ++pp)
  //   {
  //     localIndex const p = activeParticleIndices[pp];

  //     std::ofstream file;
  //     file.open( "reactionHistory.csv", std::ios::out); // | std::ios::app );
  //     if( file.fail() )
  //       throw std::ios_base::failure( std::strerror( errno ) );
  //     //make sure write fails with exception if something is wrong
  //     file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );
  //     file << "time, F00, F11, F22, length_x, length_y, length_z, Rx-, Rx+, Ry-, Ry+, Rz-, Rz+, L00, L11, L22" << std::endl;
  //     file << std::setprecision( std::numeric_limits< long double >::digits10 )
  //         << 0.0 << ","
  //         << 1.0 << "," << 1.0 << "," << 1.0 << ","
  //         << m_domainExtent[0] << "," << m_domainExtent[1] << "," << m_domainExtent[2] << ","
  //         << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ","
  //         << 0.0 << "," << 0.0 << "," << 0.0
  //         << std::endl;
  //   }
  // } );
}

void SolidMechanicsMPM::computeBoxMetrics( ParticleManager & particleManager,
                                           arrayView1d< real64 > currentStress,
                                           real64 & boxMaterialVolume )
{
  // real64 boxStress[3]; // we sum stress * volume in particles, additive sync, then divide by box volume.
  // boxMaterialVolume = 0.0;

  RAJA::MultiReduceSum< RAJA::seq_multi_reduce, real64 > boxStress ( 3, 0.0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > boxMaterialVolumeReduce( 0.0 );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get fields
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();

    // Accumulate values
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      boxMaterialVolumeReduce += particleVolume[p];

      // Only need diagonal components
      for( int i=0; i < 3; ++i )
      {
        boxStress[i] += particleStress[p][i] * particleVolume[p]; // volume weighted average, will normalize later.
      }
    } );
  } );

  // Additive sync: sxx, syy, szz, sxy, syz, sxz, mass, particle volume, damage
  real64 boxSums[4];
  boxSums[0] = boxStress[0].get();       // sig_xx * volume
  boxSums[1] = boxStress[1].get();       // sig_yy * volume
  boxSums[2] = boxStress[2].get();       // sig_zz * volume
  boxSums[3] = boxMaterialVolumeReduce.get();

  // Do an MPI sync to total these values and write from proc0 to a file.  Also compute global F
  // so file is directly plottable in excel as CSV or something.
  for( localIndex i = 0; i < 4; ++i )
  {
    real64 localSum = boxSums[i];
    real64 globalSum;
    MPI_Allreduce( &localSum,
                   &globalSum,
                   1,
                   MPI_DOUBLE,
                   MPI_SUM,
                   MPI_COMM_GEOS );
    boxSums[i] = globalSum;
  }

  real64 boxVolume = m_domainExtent[0] * m_domainExtent[1] * m_domainExtent[2];

  currentStress[0] = boxSums[0] / boxVolume;
  currentStress[1] = boxSums[1] / boxVolume;
  currentStress[2] = boxSums[2] / boxVolume;

  boxMaterialVolume = boxSums[3];
}

void SolidMechanicsMPM::stressControl( real64 dt,
                                       ParticleManager & particleManager,
                                       SpatialPartition & GEOS_UNUSED_PARAM( partition ) )
{
  GEOS_MARK_FUNCTION;

  // arrayView1d< int const > const periodic = partition.getPeriodic();

  real64 targetStress[3];
  LvArray::tensorOps::copy< 3 >(targetStress, m_domainStress);

  array1d< real64 > currentStress;
  currentStress.resize( 3 );
  LvArray::tensorOps::fill< 3 >( currentStress, 0.0 );
  // if( periodic[0] == 1 || periodic[1] == 1 || periodic[2] == 1 )
  // {
  real64 boxMaterialVolume;
  computeBoxMetrics( particleManager,
                     currentStress,
                     boxMaterialVolume );
  // }

  // CC: TODO Still use box stress but enfore Lmax for each direction

  // // Non periodic directions should use boundary reactions instead of box stresses
  // for( int i=0; i<m_numDims; ++i)
  // {
  //   if( !periodic[i] )
  //   {
  //     real64 area = 1;
  //     for( int j=0; j< m_numDims; ++j)
  //     {
  //       if( j == i )
  //       {
  //         continue;
  //       }
  //       area *= m_domainExtent[j];
  //     }
  //     // x-, x+, y-, y+, z-, z+
  //     currentStress[i] = ( m_globalFaceReactions[2*i + 1] * m_xGlobalMax[i] + m_globalFaceReactions[ 2*i ] * m_xGlobalMin[i] ) / ( area * m_domainExtent[i] );
  //   }
  // }

  // Uses maximum bulk modulus ( lowest effective PID gains ) of all materials
  real64 maximumBulkModulus = 0.0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    const ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    array1d< real64 > bulkModulus;
    string constitutiveModelName = constitutiveModel.getCatalogName();
    if( constitutiveModelName == "Hyperelastic" ){
      const Hyperelastic & hyperelastic = dynamic_cast< const Hyperelastic & >( constitutiveModel );
      bulkModulus = hyperelastic.bulkModulus();
    }

    if( constitutiveModelName == "HyperelasticMMS" || constitutiveModelName == "Chiumenti" ){
      const HyperelasticMMS & hyperelasticMMS = dynamic_cast< const HyperelasticMMS & >( constitutiveModel );
      arrayView1d< real64 const > const lambda = hyperelasticMMS.lambda();
      arrayView1d< real64 const > const shearModulus = hyperelasticMMS.shearModulus();
      bulkModulus.resize(lambda.size());
      forAll< serialPolicy >( activeParticleIndices.size(), [=, &bulkModulus] GEOS_HOST ( localIndex const pp ) // Could be reduction instead of a loop
      {
        localIndex const p = activeParticleIndices[pp];
        bulkModulus[p] = conversions::lameConstants::toBulkMod( lambda[p], shearModulus[p] );
      } );
    }

    if( constitutiveModelName == "ElasticIsotropic" || constitutiveModelName == "CeramicDamage" || constitutiveModelName == "StrainHardeningPolymer"  || constitutiveModelName == "VonMisesJ" ){
      const ElasticIsotropic & elasticIsotropic = dynamic_cast< const ElasticIsotropic & >( constitutiveModel );
      bulkModulus = elasticIsotropic.bulkModulus();
    }

    if( constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      const ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< const ElasticTransverseIsotropic & >( constitutiveModel );
      bulkModulus = elasticTransverseIsotropic.effectiveBulkModulus();
    }

    if( constitutiveModelName == "Graphite" )
    {
      const Graphite & graphite = dynamic_cast< const Graphite & >( constitutiveModel );
      bulkModulus = graphite.effectiveBulkModulus();
    }

    if( constitutiveModelName == "Geomechanics" )
    {
      const Geomechanics & geomechanics = dynamic_cast< const Geomechanics & >( constitutiveModel );
      bulkModulus = geomechanics.bulkModulus();
    }

    forAll< serialPolicy >( activeParticleIndices.size(), [=, &maximumBulkModulus] GEOS_HOST ( localIndex const pp ) // Could be reduction instead of a loop
    {
      localIndex const p = activeParticleIndices[pp];
      maximumBulkModulus = fmax( maximumBulkModulus, bulkModulus[p] );
    } );
  } );

  // Do a reduce of maximumBulkModulus between ranks so stress control calculations are consistent
  real64 maximumBulkModulusThisRank = maximumBulkModulus;
  real64 maximumBulkModulusAllRanks = 0;
  maximumBulkModulusAllRanks = MpiWrapper::allReduce( maximumBulkModulusThisRank,
                                                      MpiWrapper::Reduction::Max,
                                                      MPI_COMM_GEOS );
  maximumBulkModulus = maximumBulkModulusAllRanks;

  GEOS_ERROR_IF( maximumBulkModulus <= 0.0, "At least one material must have a positive bulk or effective bulk modulus for stress control!" );

  // Could some numerical artifact ever produce a negative relative density? Do I need a check for that or clip is to 0
  real64 relativeDensity = fmin( 1.0, boxMaterialVolume / ( m_domainExtent[0] * m_domainExtent[1] * m_domainExtent[2] ) );

	// This will drive the response towards the desired stress but may be
	// unstable.
  // CC: TODO Use 1- domain porosity to make uncompacted region more stable
	real64 stressControlKp = relativeDensity * m_stressControlKp / ( maximumBulkModulus * dt );
	real64 stressControlKd = relativeDensity * m_stressControlKd / ( maximumBulkModulus );
	real64 stressControlKi = relativeDensity * m_stressControlKi / ( maximumBulkModulus * dt * dt );

	real64 error[3];
  LvArray::tensorOps::copy< 3 >(error, targetStress);
  LvArray::tensorOps::subtract< 3 >(error, currentStress);

	real64 dedt[3];
  LvArray::tensorOps::copy< 3 >( dedt, error);
  LvArray::tensorOps::subtract< 3 >( dedt, m_stressControlLastError);
  LvArray::tensorOps::scale< 3 >( dedt, 1.0/dt);

	real64 stressControlPTerm[3];
  LvArray::tensorOps::copy< 3 >( stressControlPTerm, error);
  LvArray::tensorOps::scale< 3 >( stressControlPTerm, stressControlKp);

  real64 stressControlITerm[3];
  LvArray::tensorOps::copy< 3 >( stressControlITerm, error);
  LvArray::tensorOps::scale< 3 >( stressControlITerm, stressControlKi*dt);
  LvArray::tensorOps::add< 3 >( m_stressControlITerm, stressControlITerm );

  real64 stressControlDTerm[3];
  LvArray::tensorOps::copy< 3 >( stressControlDTerm, dedt);
  LvArray::tensorOps::scale< 3 >(stressControlDTerm, stressControlKd);

	real64 L_new[3];
  LvArray::tensorOps::copy< 3 >( L_new, stressControlPTerm );
  LvArray::tensorOps::add< 3 >( L_new, m_stressControlITerm );
  LvArray::tensorOps::add< 3 >( L_new, stressControlDTerm );

	LvArray::tensorOps::copy< 3 >(m_stressControlLastError, error);

	// Limit L so boundary motion doesn't exceed CFL limit
	// boundaryV = strainRate*(m_xGlobalMax[0] - m_xGlobalMin[0]) = CFL*m_hEl[0]/dt
	real64 Lxx_min,
         Lxx_max,
         Lyy_min,
         Lyy_max,
         Lzz_min,
         Lzz_max;

	Lxx_max = m_cflFactor * ( m_hEl[0] / dt) / ( m_xGlobalMax[0] - m_xGlobalMin[0]);
	Lxx_min = -1.0 * Lxx_max;
	L_new[0] = LvArray::math::min( Lxx_max , LvArray::math::max( L_new[0], Lxx_min ) );

	Lyy_max = m_cflFactor * ( m_hEl[1] / dt) / ( m_xGlobalMax[1] - m_xGlobalMin[1]);
	Lyy_min = -1.0 * Lyy_max;
	L_new[1] = LvArray::math::min( Lyy_max , LvArray::math::max( L_new[1], Lyy_min ) );

	Lzz_max = m_cflFactor * ( m_hEl[2] / dt) / ( m_xGlobalMax[2] - m_xGlobalMin[2]);
	Lzz_min = -1.0 * Lzz_max;
	L_new[2] = LvArray::math::min( Lzz_max , LvArray::math::max( L_new[2], Lzz_min ) );

	for(int i=0; i < m_numDims; ++i )
	{
		if( m_stressControl[i] == 1 )
		{
			m_domainL[i] = L_new[i];
			m_domainF[i] += m_domainL[i]*m_domainF[i]*dt;
		}
	}
}


void SolidMechanicsMPM::applySuperimposedVelocityGradient( const real64 dt,
                                                           ParticleManager & particleManager,
                                                           SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  int const numDims = m_numDims;
  real64 domainL[3];
  LvArray::tensorOps::copy< 3 >( domainL, m_domainL );
  arrayView1d< int const > const periodic = partition.getPeriodic();

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_DEVICE ( localIndex const pp )
       {
        localIndex const p = activeParticleIndices[pp];

        for(int i=0; i < numDims; ++i)
        {
          if( periodic[i] )
          {
            particleVelocityGradient[p][i][i] += domainL[i];
            particlePosition[p][i] += particlePosition[p][i] * domainL[i] * dt;
          }
        }
      } ); // particle loop

      particlePosition.move( LvArray::MemorySpace::host );
      particleVelocityGradient.move( LvArray::MemorySpace::host );

  } ); // subregion loop
}

void SolidMechanicsMPM::initializeGridFields( NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  int const numVelocityFields = m_numVelocityFields;;

  // Scalars
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView2d< real64 > const gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 > const gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView1d< real64 > const gridSurfaceMass = nodeManager.getReference< array1d< real64 > >( viewKeyStruct::gridSurfaceMassString() );
  arrayView2d< real64 > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
  arrayView1d< int > const gridCohesiveNode = nodeManager.getReference< array1d< int > >( viewKeyStruct::gridCohesiveNodeString() );

  // Vectors
  arrayView3d< real64 > const gridDisplacement = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDisplacementString() );
  arrayView3d< real64 > const gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 > const gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVelocityString() );
  arrayView3d< real64 > const gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );
  arrayView2d< real64 > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView3d< real64 > const gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 > const gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );
  arrayView3d< real64 > const gridSurfaceTensionForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceTensionForceString() );
  arrayView3d< real64 > const gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridContactForceString() );

  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView2d< real64 > const gridSurfaceNormalWeights = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightsString() );
  arrayView2d< real64 > const gridSurfaceNormalWeightNormalization = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceNormalWeightNormalizationString() );
  arrayView3d< real64 > const gridParticleSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridParticleMappedSurfaceNormalString() );
  arrayView3d< real64 > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );
  arrayView2d< real64 > const gridSurfaceArea = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceAreaString() );

  arrayView3d< real64 > const gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );

  arrayView2d< real64 > const gridExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridExplicitSurfaceNormalString() );
  arrayView2d< int > const gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );
  arrayView3d< real64 > const gridCohesiveArea = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCohesiveAreaString() );
  arrayView3d< real64 > const gridCohesiveForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCohesiveForceString() );

  arrayView2d< real64 > const gridPrincipalExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridPrincipalExplicitSurfaceNormalString() );

  arrayView1d< globalIndex > const gridMaxMappedParticleID = nodeManager.getReference< array1d< globalIndex > >( viewKeyStruct::gridMaxMappedParticleIDString() );

  arrayView2d< real64 > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  // Development
  arrayView3d< real64 > const gridLRSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfaceNormalString() );
  arrayView3d< real64 > const gridLRSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfacePositionString() );

  forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    gridCohesiveNode[g] = 0;
    gridSurfaceMass[g] = 0.0;

    gridMaxMappedParticleID[g] = -1;

    for( int i = 0; i < 3; ++i )
    {
      gridDamageGradient[g][i] = 0.0;
      gridExplicitSurfaceNormal[g][i] = 0.0;
      gridPrincipalExplicitSurfaceNormal[g][i] = 0.0;
    }

    LvArray::tensorOps::fill< 6 >( gridBoreholeStress[g], 0.0 );

    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      gridSurfaceFieldMass[g][fieldIndex] = 0.0;
      gridCohesiveFieldFlag[g][fieldIndex] = 0;

      gridMass[g][fieldIndex] = 0.0;
      gridMaterialVolume[g][fieldIndex] = 0.0;
      gridDamage[g][fieldIndex] = 0.0;
      gridMaxDamage[g][fieldIndex] = 0.0;
      gridMassWeightedDamage[g][fieldIndex] = 0.0;

      gridSurfaceNormalWeightNormalization[g][fieldIndex] = 0;
      gridSurfaceNormalWeights[g][fieldIndex] = 0.0;

      gridSurfaceArea[g][fieldIndex] = 0.0;

      LvArray::tensorOps::fill< 3 >( gridLRSurfaceNormal[g][fieldIndex], 0.0 );
      LvArray::tensorOps::fill< 3 >( gridLRSurfacePosition[g][fieldIndex], 0.0 );

      for( int i = 0; i < 3; ++i )
      {
        gridCenterOfMass[g][fieldIndex][i] = 0.0;
        gridDisplacement[g][fieldIndex][i] = 0.0;
        gridCenterOfVolume[g][fieldIndex][i] = 0.0;
        gridVelocity[g][fieldIndex][i] = 0.0;
        gridDVelocity[g][fieldIndex][i] = 0.0;
        gridMomentum[g][fieldIndex][i] = 0.0;
        gridAcceleration[g][fieldIndex][i] = 0.0;
        gridInternalForce[g][fieldIndex][i] = 0.0;
        gridExternalForce[g][fieldIndex][i] = 0.0;
        gridSurfaceTensionForce[g][fieldIndex][i] = 0.0;
        gridContactForce[g][fieldIndex][i] = 0.0;
        gridNormalStress[g][fieldIndex][i] = 0.0;
        gridCohesiveArea[g][fieldIndex][i] = 0.0;
        gridCohesiveForce[g][fieldIndex][i] = 0.0;
        gridParticleSurfaceNormal[g][fieldIndex][i] = 0.0;
        gridSurfaceNormal[g][fieldIndex][i] = 0.0;
        gridSurfacePosition[g][fieldIndex][i] = 0.0;
      }
    }
  } );
}

void SolidMechanicsMPM::boundaryConditionUpdate( real64 dt, real64 time_n )
{
  GEOS_MARK_FUNCTION;

  int bcInterval = 0;

  for( localIndex i = 0; i < m_bcTable.size( 0 ); ++i ) // Naive method for determining what part of BC table we're currently in, can
                                                        // optimize later (TODO)
  {
    if( time_n + 0.5 * dt > m_bcTable[i][0] )
    {
      bcInterval = i;
    }
  }

  for( int i = 0; i < 6; ++i )
  {
    m_boundaryConditionTypes[i] = m_bcTable[bcInterval][i+1];
  }
}

void SolidMechanicsMPM::projectParticleSurfaceNormalsToGrid( DomainPartition & domain,
                                                             ParticleManager & particleManager,
                                                             NodeManager & nodeManager,
                                                             MeshLevel & mesh )
{
  real64 const smallMass = m_smallMass;

  // Grid fields
  arrayView1d< globalIndex > const gridMaxMappedParticleID = nodeManager.getReference< array1d< globalIndex > >( viewKeyStruct::gridMaxMappedParticleIDString() );
  arrayView2d< real64 > const gridPrincipalExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridPrincipalExplicitSurfaceNormalString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get particle fields
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    // arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    // Get nodes this particle maps to
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( localIndex gg=0; gg < 8 * numberOfVerticesPerParticle; gg++ )
      {
        // TODO Change this to use on the fly mappings
        localIndex const & g = mappedNodes[pp][gg];
        RAJA::atomicMax( parallelDeviceAtomic{}, &gridMaxMappedParticleID[g], particleID[p] );
      }
    } );
    ++subRegionIndex;
  } );

  syncGridFields( { viewKeyStruct::gridMaxMappedParticleIDString() }, domain, nodeManager, mesh, MPI_MAX );

  arrayView2d< real64 > const gridExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridExplicitSurfaceNormalString() );

  // Normalize mapped normals by mapped mass
  localIndex const numNodes = nodeManager.size();
  // auto globalToLocalMap = nodeManager.globalToLocalMap();
  forAll< serialPolicy >( numNodes, [&] GEOS_HOST ( localIndex const g )
  {
    LvArray::tensorOps::fill< 3 >( gridExplicitSurfaceNormal[g], 0.0 );

    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      auto globalToLocalMap = subRegion.globalToLocalMap();
      auto iter = globalToLocalMap.find( gridMaxMappedParticleID[g] );
      if( iter == globalToLocalMap.end() )
      {
        return;
      }

      // Even if ghost particles are writting in non home partition it would only result in doubling of normal length which doesn't change its use (magnitude is irrelevant, but could normalize if we want)
      arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
      localIndex const p = globalToLocalMap.at( gridMaxMappedParticleID[g] );
      LvArray::tensorOps::copy< 3 >( gridExplicitSurfaceNormal[g], particleSurfaceNormal[p] );
    } );
  } );

  syncGridFields( { viewKeyStruct::gridPrincipalExplicitSurfaceNormalString() }, domain, nodeManager, mesh, MPI_SUM );

  arrayView1d< real64 > const gridSurfaceMass = nodeManager.getReference< array1d< real64 > >( viewKeyStruct::gridSurfaceMassString() );

  subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get particle fields
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    // arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    // Get nodes this particle maps to
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( localIndex gg = 0; gg < 8 * numberOfVerticesPerParticle; gg++ )
      {
        localIndex const & g = mappedNodes[pp][gg];

        real64 index = ( LvArray::tensorOps::AiBi< 3 >( gridPrincipalExplicitSurfaceNormal[g], particleSurfaceNormal[p] ) < 0.0 ) ? -1.0 : 1.0;

        if( LvArray::tensorOps::l2NormSquared< 3 >( particleSurfaceNormal[p] ) > 1e-12 )
        {
          gridSurfaceMass[g] += particleMass[p] * shapeFunctionValues[pp][gg];
          for(int i = 0; i < 3; ++i)
          {
            gridExplicitSurfaceNormal[g][i] += index * particleMass[p] * shapeFunctionValues[pp][gg] * particleSurfaceNormal[p][i];
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop

  syncGridFields( { viewKeyStruct::gridSurfaceMassString(), viewKeyStruct::gridExplicitSurfaceNormalString() }, domain, nodeManager, mesh, MPI_SUM );

  // Normalize mapped normals by mapped mass
  // int const numNodes = nodeManager.size();
  forAll< serialPolicy >( numNodes, [=] GEOS_HOST ( localIndex const g )
  {
    if( gridSurfaceMass[g] > smallMass )
    {
      LvArray::tensorOps::normalize< 3 >( gridExplicitSurfaceNormal[g] );
    }
  } );
}

void SolidMechanicsMPM::initializeCohesiveReferenceConfiguration( DomainPartition & domain,
                                                                  ParticleManager & particleManager,
                                                                  NodeManager & nodeManager,
                                                                  MeshLevel & mesh )
{
  real64 const smallMass = m_smallMass;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numVelocityFields = m_numVelocityFields;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  real64 xLocalMax[3];
  LvArray::tensorOps::copy< 3 >( xLocalMax, m_xLocalMax );
  real64 xGlobalMin[3];
  LvArray::tensorOps::copy< 3 >( xGlobalMin, m_xGlobalMin );
  real64 xGlobalMax[3];
  LvArray::tensorOps::copy< 3 >( xGlobalMax, m_xGlobalMax );

  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const gridExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridExplicitSurfaceNormalString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    // arrayView1d< int const > const particleRank = subRegion.getParticleRank();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const> const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    arrayView2d< int > const particleCohesiveFieldMapping = subRegion.getField< fields::mpm::particleCohesiveFieldMapping >();

    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      real64 particleContributionToGrid;

      if( static_cast< SurfaceFlag >( particleSurfaceFlag[p] ) == SurfaceFlag::Cohesive )
      {
        for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
        {
          localIndex const mappedNode = mappedNodes[pp][g];
          localIndex nodeFlag = 0;
          if( damageFieldPartitioning == 1 && LvArray::tensorOps::l2Norm< 3 >( gridExplicitSurfaceNormal[mappedNode] ) > 1e-12 )
          {
            nodeFlag = ( LvArray::tensorOps::AiBi< 3 >( gridExplicitSurfaceNormal[mappedNode], particleSurfaceNormal[p] ) < 0.0 ) ? 1 : 0; // 0 for "A" field, 1 for "B" field
          }

          localIndex const fieldIndex = nodeFlag * m_numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

          particleCohesiveFieldMapping[p][g] = fieldIndex;

          particleContributionToGrid = particleMass[p] * shapeFunctionValues[pp][g];
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMass[mappedNode][fieldIndex], particleContributionToGrid );
        }
      }
    });
    ++subRegionIndex;
  });

  // Grid sync
  syncGridFields( { viewKeyStruct::gridMassString() }, domain, nodeManager, mesh, MPI_SUM );

  // Any nodes with both gridMass for fields A and B are nodes that could belong to a cohesive interface
  array1d< globalIndex > interfaceGridNodes;

  localIndex const numNodes = nodeManager.size();
  arrayView1d< globalIndex > localToGlobalMap = nodeManager.localToGlobalMap();
  forAll< serialPolicy >( numNodes, [=, &interfaceGridNodes] GEOS_HOST ( localIndex const g )
  {

    for( localIndex A = 0; A < numVelocityFields - 1; ++A )
    {
      for( localIndex B = A + 1; B < numVelocityFields; ++B )
      {
        bool active = ( gridMass[g][A] > smallMass ) && ( gridMass[g][B] > smallMass );

        if( active )
        {
          interfaceGridNodes.emplace_back( localToGlobalMap[g] );
        }
      }
    }

  } );

  // Collect all the nodes that belong to cohesive interfaces and distribute them to all the processes for mapping particles to
  // reference grid when particle advection triggers repartitioning
  array1d< int > dataSizes( MpiWrapper::commSize() );
  MpiWrapper::allGather( LvArray::integerConversion< int >( interfaceGridNodes.size() ), dataSizes, MPI_COMM_GEOS );
  int const totalDataSize = std::accumulate( dataSizes.begin(), dataSizes.end(), 0 );

  // Once the MPI exchange is done, `allData` will contain all the data of all the MPI ranks.
  // We want all ranks to get all the data. But each rank may have a different size of information.
  // Therefore, we use `allgatherv` that does not impose the same size across ranks like `allgather` does.
  stdVector< globalIndex > allData( totalDataSize );
  // `displacements` is the offset (relative to the receive buffer) to store the data for each rank.
  stdVector< int > displacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( dataSizes.begin(), dataSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv( interfaceGridNodes.data(), interfaceGridNodes.size(), allData.data(), dataSizes.data(), displacements.data(), MPI_COMM_GEOS );

 // Sort the grid indices and move any duplicates to the end.
  std::ptrdiff_t const numUniqueValues = LvArray::sortedArrayManipulation::makeSortedUnique( allData.begin(),
                                                                                             allData.end() );
  // Move unique global grid node indices to member variable
  m_cohesiveNodeGlobalIndices.insert( allData.begin(), allData.begin() + numUniqueValues );

  int numCohesiveNodes = m_cohesiveNodeGlobalIndices.size();

  m_referenceCohesiveGridNodePositions.resize( numCohesiveNodes, 3 );
  m_referenceCohesiveGridNodePartitioningSurfaceNormals.resize( numCohesiveNodes, 3 );

  arrayView2d< real64 > const referenceCohesiveGridNodePositions = m_referenceCohesiveGridNodePositions;
  arrayView2d< real64 > const referenceCohesiveGridNodePartitioningSurfaceNormals = m_referenceCohesiveGridNodePartitioningSurfaceNormals;

  array1d< int > localCohesiveGridNodes;

  forAll< serialPolicy >( numNodes, [=, &localCohesiveGridNodes] GEOS_HOST ( localIndex const g )
  {
    globalIndex const mappedNode = localToGlobalMap[ g ];

    if( m_cohesiveNodeGlobalIndices.contains(  mappedNode ) )
    {
      // CC: TODO must be a better way to find index in temp arrays
      localIndex nodeIndex = 0;
      for( int n = 0; n < numCohesiveNodes; ++n )
      {
        if( m_cohesiveNodeGlobalIndices[n] == mappedNode )
        {
          nodeIndex = n;
          break;
        }
      }

      for( int i = 0; i < 3; ++i)
      {
        referenceCohesiveGridNodePositions[nodeIndex][i] = gridPosition[g][i];
        referenceCohesiveGridNodePartitioningSurfaceNormals[nodeIndex][i] = gridExplicitSurfaceNormal[g][i];
      }

      localCohesiveGridNodes.emplace_back( nodeIndex );
    }

    for(localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
    {
      gridMass[g][fieldIndex] = 0.0;
    }
  } );

  // Sync initial grid positions with rank 0
  int const numRanks = MpiWrapper::commSize();
  int rank = MpiWrapper::commRank( MPI_COMM_GEOS );

  array1d< MPI_Request > mpiRequestIndices( numRanks );
  array1d< MPI_Request > mpiRequestPositions( numRanks );
  array1d< MPI_Request > mpiRequestSurfaceNormals( numRanks );

  array1d< MPI_Status > mpiStatusIndices( numRanks );
  array1d< MPI_Status > mpiStatusPositions( numRanks );
  array1d< MPI_Status > mpiStatusSurfaceNormals( numRanks );

  array1d< int > gridNodeIndices;
  if( rank != 0 )
  {
    // Send list of indices to overwrite nodal positions
    mpiRequestIndices[rank] = MPI_REQUEST_NULL;
    MpiWrapper::iSend( localCohesiveGridNodes.data(),
                       localCohesiveGridNodes.size(),
                       0,
                       0,
                       MPI_COMM_GEOS,
                       &mpiRequestIndices[rank] );

    // Send nodal positions
    mpiRequestPositions[rank] = MPI_REQUEST_NULL;
    MpiWrapper::iSend( m_referenceCohesiveGridNodePositions.data(),
                       m_referenceCohesiveGridNodePositions.size(),
                       0,
                       1,
                       MPI_COMM_GEOS,
                       &mpiRequestPositions[rank] );

    // Send nodal positions
    mpiRequestSurfaceNormals[rank] = MPI_REQUEST_NULL;
    MpiWrapper::iSend( m_referenceCohesiveGridNodePartitioningSurfaceNormals.data(),
                       m_referenceCohesiveGridNodePartitioningSurfaceNormals.size(),
                       0,
                       1,
                       MPI_COMM_GEOS,
                       &mpiRequestSurfaceNormals[rank] );
  }
  else
  {
    for( int r = 1; r < numRanks; r++)
    {
      MpiWrapper::recv( gridNodeIndices,
                        r,
                        0,
                        MPI_COMM_GEOS,
                        &mpiStatusIndices[r] );

      array2d< real64 > gridNodePositions( numCohesiveNodes, 3 );
      mpiRequestPositions[r] = MPI_REQUEST_NULL;
      MpiWrapper::iRecv( gridNodePositions.data(),
                         gridNodePositions.size(),
                         r,
                         1,
                         MPI_COMM_GEOS,
                         &mpiRequestPositions[r] );

      MpiWrapper::wait( &mpiRequestPositions[r], &mpiStatusPositions[r]);

      // Send nodal surface normals
      array2d< real64 > gridNodeSurfaceNormals( numCohesiveNodes, 3 );
      mpiRequestSurfaceNormals[r] = MPI_REQUEST_NULL;
      MpiWrapper::iRecv( gridNodeSurfaceNormals.data(),
                         gridNodeSurfaceNormals.size(),
                         r,
                         1,
                         MPI_COMM_GEOS,
                         &mpiRequestSurfaceNormals[r] );

      MpiWrapper::wait( &mpiRequestSurfaceNormals[r], &mpiStatusSurfaceNormals[r] );

      // Combine grid positions
      for( int g = 0; g < gridNodeIndices.size(); ++g )
      {
        for(int k = 0; k < 3; ++k){
          referenceCohesiveGridNodePositions[gridNodeIndices[g]][k] = gridNodePositions[gridNodeIndices[g]][k];
          referenceCohesiveGridNodePartitioningSurfaceNormals[gridNodeIndices[g]][k] = gridNodeSurfaceNormals[gridNodeIndices[g]][k];
        }
      }
    }
  }

  MpiWrapper::barrier();

  // Scatter referenceCohesiveGridNodePositions to all other ranks
  MpiWrapper::bcast( m_referenceCohesiveGridNodePositions.data(),
                     m_referenceCohesiveGridNodePositions.size(),
                     0,
                     MPI_COMM_GEOS );

  // Scatter referenceCohesiveGridNodePositions to all other ranks
  MpiWrapper::bcast( m_referenceCohesiveGridNodePartitioningSurfaceNormals.data(),
                     m_referenceCohesiveGridNodePartitioningSurfaceNormals.size(),
                     0,
                     MPI_COMM_GEOS );

  // Using mapping back to particles to identify particles involved in cohesive zone calculations and precompute shape functions from initial conditions
  subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();

    arrayView1d< int > particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    arrayView2d< globalIndex > particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();
    arrayView2d< real64 > particleReferenceShapeFunctionValues = subRegion.getField< fields::mpm::particleReferenceShapeFunctionValues >();

    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();

    arrayView2d< real64 > const particleCohesiveReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleCohesiveReferenceSurfaceNormal >();
    arrayView2d< real64 > const particleCohesiveReferenceSurfacePosition = subRegion.getField< fields::mpm::particleCohesiveReferenceSurfacePosition >();

    arrayView3d< real64 > const particleCohesiveReferenceDeformationGradient = subRegion.getField< fields::mpm::particleCohesiveReferenceDeformationGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];


      // if( static_cast< SurfaceFlag >( particleSurfaceFlag[p] ) == SurfaceFlag::Cohesive )
      // {
        // Store reference deformation gradient of particles to update surface normals and positions
        LvArray::tensorOps::copy< 3, 3 >( particleCohesiveReferenceDeformationGradient[p], particleDeformationGradient[p] );
        LvArray::tensorOps::copy< 3 >( particleCohesiveReferenceSurfaceNormal[p], particleSurfaceNormal[p] );
        LvArray::tensorOps::copy< 3 >( particleCohesiveReferenceSurfacePosition[p], particleSurfacePosition[p] );

        for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
        {
          globalIndex const mappedNode = localToGlobalMap[ mappedNodes[pp][g] ];

          if( m_cohesiveNodeGlobalIndices.contains( mappedNode ) )
          {
            if( static_cast< SurfaceFlag >( particleSurfaceFlag[p] ) == SurfaceFlag::Cohesive )
            {
              particleCohesiveZoneFlag[p] = 1;
            }

            // Store particle shape function data for future cohesive computations
            for( int gg = 0; gg < numberOfVerticesPerParticle * 8; gg++ )
            {
              particleReferenceMappedNodes[p][gg] = localToGlobalMap[ mappedNodes[pp][gg] ];
              particleReferenceShapeFunctionValues[p][gg] = shapeFunctionValues[pp][gg];
            }

            // If at least one corner is involved no need to look further
            break;
          }
        }
      // }
    } );
    ++subRegionIndex;
  } );

  // Now we need to compute the grid area at cohesive initialization
  // We map the surface position of each particle ( vector from particle center to interface surface ), this is also the particle surface normal direction
  array2d< real64 > tempGridMassLocal( numCohesiveNodes, m_numVelocityFields );
  array1d< real64 > tempGridVolumeLocal( numCohesiveNodes );
  array3d< real64 > tempGridParticleSurfaceNormalLocal( numCohesiveNodes, m_numVelocityFields, 3 );
  array3d< real64 > tempGridSurfacePositionLocal( numCohesiveNodes, m_numVelocityFields, 3 );

  // Initialize temporary grid fields to zero
  for( int g  = 0; g < numCohesiveNodes; ++g)
  {
    tempGridVolumeLocal[g] = 0.0;
    for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex )
    {
      tempGridMassLocal[g][fieldIndex] = 0.0;
      for( int i = 0; i < m_numDims; ++i)
      {
        tempGridParticleSurfaceNormalLocal[g][fieldIndex][i] = 0.0;
        tempGridSurfacePositionLocal[g][fieldIndex][i] = 0.0;
      }
    }
  }

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();

    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    arrayView2d< globalIndex const > const particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();
    arrayView2d< real64 const > const particleReferenceShapeFunctionValues = subRegion.getField< fields::mpm::particleReferenceShapeFunctionValues >();
    arrayView2d< int const > const particleCohesiveFieldMapping = subRegion.getField< fields::mpm::particleCohesiveFieldMapping >();

    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    auto globalToLocalMap = nodeManager.globalToLocalMap();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &tempGridMassLocal, &tempGridVolumeLocal, &tempGridParticleSurfaceNormalLocal, &tempGridSurfacePositionLocal] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Map particle displacement to grid
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        globalIndex const mappedNode = particleReferenceMappedNodes[p][g];

        if( m_cohesiveNodeGlobalIndices.contains( mappedNode ) )
        {
          real64 shapeFunctionValue = particleReferenceShapeFunctionValues[p][g];

          // CC: TODO must be a better way to find index in temp arrays
          localIndex nodeIndex = 0;
          for( int n = 0; n < numCohesiveNodes; ++n )
          {
            if( m_cohesiveNodeGlobalIndices[n] == mappedNode )
            {
              nodeIndex = n;
              break;
            }
          }

          localIndex const fieldIndex = particleCohesiveFieldMapping[p][g];

          if( particleCohesiveZoneFlag[p] == 1 )
          {
            tempGridMassLocal[nodeIndex][fieldIndex] += particleMass[p] * shapeFunctionValue;
            // tempGridVolumeLocal[nodeIndex] += particleVolume[p] * shapeFunctionValue;

            for( int i  = 0; i < m_numDims; ++i ){
              tempGridParticleSurfaceNormalLocal[nodeIndex][fieldIndex][i] += particleMass[p] * particleSurfaceNormal[p][i] * shapeFunctionValue;
              tempGridSurfacePositionLocal[nodeIndex][fieldIndex][i] += particleMass[p] * (  particlePosition[p][i] - referenceCohesiveGridNodePositions[nodeIndex][i] + particleSurfacePosition[p][i] ) * shapeFunctionValue;
            }

          }
          tempGridVolumeLocal[nodeIndex] += particleVolume[p] * shapeFunctionValue;
        }
      }
    });
  });

  // Sync temporary grid fields
  array2d< real64 > tempGridMassGlobal( numCohesiveNodes, m_numVelocityFields );
  array1d< real64 > tempGridVolumeGlobal( numCohesiveNodes );
  array3d< real64 > tempGridParticleSurfaceNormalGlobal( numCohesiveNodes, m_numVelocityFields, 3 );
  array3d< real64 > tempGridSurfacePositionGlobal( numCohesiveNodes, m_numVelocityFields, 3 );

  MpiWrapper::allReduce( tempGridMassLocal,
                         tempGridMassGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridVolumeLocal,
                         tempGridVolumeGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridParticleSurfaceNormalLocal,
                         tempGridParticleSurfaceNormalGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridSurfacePositionLocal,
                         tempGridSurfacePositionGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  forAll< serialPolicy >( numCohesiveNodes, [=, &tempGridParticleSurfaceNormalGlobal, & tempGridSurfacePositionGlobal] GEOS_HOST ( localIndex const g )
  {
    for(localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
    {
      if ( tempGridMassGlobal[g][fieldIndex] > m_smallMass )
      {
        for( int i = 0; i < m_numDims; ++i ){
          tempGridParticleSurfaceNormalGlobal[g][fieldIndex][i] /= tempGridMassGlobal[g][fieldIndex];
          tempGridSurfacePositionGlobal[g][fieldIndex][i] /= tempGridMassGlobal[g][fieldIndex];
        }
      }
    }
  } );

  m_referenceCohesiveGridNodeAreas.resize( numCohesiveNodes, m_numVelocityFields ); // This needs to be solver field
  m_cohesiveGridNodeDamages.resize( numCohesiveNodes, m_numVelocityFields ); // TODO: eventually this should account for pairings of velocity fields (e.g. multiple cohesive laws defined mapping to different field pairings)
  m_referenceCohesiveGridNodeSurfaceNormals.resize( numCohesiveNodes, m_numVelocityFields, 3 );

  m_maxCohesiveGridNodeNormalDisplacement.resize( numCohesiveNodes, m_numVelocityFields, m_numVelocityFields );
  m_maxCohesiveGridNodeTangentialDisplacement.resize( numCohesiveNodes, m_numVelocityFields, m_numVelocityFields );

  for( int g = 0; g < numCohesiveNodes; ++g )
  {
    for ( int a = 0; a < m_numVelocityFields; ++a )
    {
      for( int b = 0; b < m_numVelocityFields; ++b )
      {
        m_maxCohesiveGridNodeNormalDisplacement[g][a][b] = 0.0;
        m_maxCohesiveGridNodeTangentialDisplacement[g][a][b] = 0.0;
      }
    }
  }

  arrayView2d< real64 > const referenceCohesiveGridNodeAreas = m_referenceCohesiveGridNodeAreas;
  arrayView2d< real64 > const cohesiveGridNodeDamages = m_cohesiveGridNodeDamages;
  arrayView3d< real64 > const referenceCohesiveGridNodeSurfaceNormals = m_referenceCohesiveGridNodeSurfaceNormals;

  real64 L = LvArray::tensorOps::l2Norm< 3 >( hEl );
  real64 dA = pow( 2 * L / m_numSurfaceIntegrationPoints, 2 );

  // Integrate shapefunctions along cohesive interface (approximated as plane from mapped normals) to determine grid area
  // forAll< serialPolicy >( numCohesiveNodes, [=, &referenceCohesiveGridNodeAreas, &referenceCohesiveGridNodeSurfaceNormals ] GEOS_HOST ( localIndex const g )
  real64 totalCohesiveSurfaceArea = 0.0;
  for( localIndex g = 0; g < numCohesiveNodes; ++g)
  {
    for(localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
    {
      referenceCohesiveGridNodeAreas[g][fieldIndex] = 0.0; //Initialize this to zero
      cohesiveGridNodeDamages[g][fieldIndex] = 0.0; //Initialize this to zero

      if( LvArray::tensorOps::l2Norm< 3 >( tempGridParticleSurfaceNormalGlobal[g][fieldIndex] ) < 1e-16)
      {
        continue;
      }

      LvArray::tensorOps::copy< 3 >( referenceCohesiveGridNodeSurfaceNormals[g][fieldIndex], tempGridParticleSurfaceNormalGlobal[g][fieldIndex] );

      // Construct basis for integration of surface plane
      real64 n[3] = {};
      LvArray::tensorOps::copy< 3 >( n, tempGridParticleSurfaceNormalGlobal[g][fieldIndex] );
      LvArray::tensorOps::normalize< 3 >( n );

      real64 distanceToSurface = fabs(LvArray::tensorOps::AiBi< 3 >( tempGridSurfacePositionGlobal[g][fieldIndex], n ) );

      real64 s1[3] = {};
      real64 s2[3] = {};
      computeOrthonormalBasis( n,
                               s1,
                               s2 );

      for( int i=0; i < m_numSurfaceIntegrationPoints; ++i )
      {
        real64 eta = 2.0 * i / ( m_numSurfaceIntegrationPoints - 1 ) - 1.0; // Natural coordinate of surface
        for( int j=0; j < m_numSurfaceIntegrationPoints; ++j )
        {
          real64 xi = 2.0 * j / ( m_numSurfaceIntegrationPoints - 1 ) - 1.0; // Natural coordinate of surface

          real64 surfacePoint[3] = {};
          surfacePoint[0] = distanceToSurface * n[0] + L * eta * s1[0] + L * xi * s2[0];
          surfacePoint[1] = distanceToSurface * n[1] + L * eta * s1[1] + L * xi * s2[1];
          surfacePoint[2] = distanceToSurface * n[2] + L * eta * s1[2] + L * xi * s2[2];

          if( ( -hEl[0] < surfacePoint[0] ) &&
              ( surfacePoint[0] < hEl[0] ) &&
              ( -hEl[1] < surfacePoint[1] ) &&
              ( surfacePoint[1] < hEl[1] ) &&
              ( -hEl[2] < surfacePoint[2] ) &&
              ( surfacePoint[2] < hEl[2] ) )
          {
            referenceCohesiveGridNodeAreas[g][fieldIndex] += ( 1 - fabs( surfacePoint[0] / hEl[0] ) ) * ( 1 - fabs( surfacePoint[1] / hEl[1]) ) * ( 1 - fabs( surfacePoint[2] / hEl[2]) );
          }
        }
      }
      referenceCohesiveGridNodeAreas[g][fieldIndex] *= dA * tempGridVolumeGlobal[g] / ( hEl[0] * hEl[1] * hEl[2] );
      totalCohesiveSurfaceArea += referenceCohesiveGridNodeAreas[g][fieldIndex];
    }
  }
  // } );

  totalCohesiveSurfaceArea /= m_numVelocityFields;
  m_polymerCZThickness = m_totalBinderVolume / totalCohesiveSurfaceArea;


  // CC: debug, for debugging temp grid variables to visualize in paraview
  arrayView3d< real64 > const gridReferenceAreaVector = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridReferenceAreaVectorString() );
  arrayView3d< real64 > const gridReferenceSurfacePosition = nodeManager.getReference< array3d < real64 > >( viewKeyStruct::gridReferenceSurfacePositionString() );
  arrayView1d< real64 > const gridReferenceMaterialVolume = nodeManager.getReference< array1d < real64 > >( viewKeyStruct::gridReferenceMaterialVolumeString() );

  forAll< serialPolicy >( nodeManager.size(), [=] GEOS_HOST ( localIndex const g )
  {
    bool isCohesive = false;
    for( int n = 0; n < numCohesiveNodes; ++n)
    {
      if( localToGlobalMap[g] == m_cohesiveNodeGlobalIndices[n] )
      {
        for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
        {
          gridReferenceMaterialVolume[g] = tempGridVolumeGlobal[n];
          for(int i = 0; i < m_numDims; ++i)
          {
            gridReferenceSurfacePosition[g][fieldIndex][i] = tempGridSurfacePositionGlobal[n][fieldIndex][i];
            gridReferenceAreaVector[g][fieldIndex][i] = referenceCohesiveGridNodeAreas[n][fieldIndex] * referenceCohesiveGridNodeSurfaceNormals[n][fieldIndex][i];
          }
        }
        isCohesive = true;
        break;
      }
    }

    if( !isCohesive )
    {
      gridReferenceMaterialVolume[g] = 0.0;
      for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
      {
        for(int i = 0; i < m_numDims; ++i)
        {
          gridReferenceSurfacePosition[g][fieldIndex][i] = 0.0;
          gridReferenceAreaVector[g][fieldIndex][i] = 0.0;
        }
      }
    }
  } );

}


// CC: Currently unused but may be useful if we perform any ray casting with particles in the future
bool SolidMechanicsMPM::interiorToParticleProjectedArea( ParticleManager & particleManager,
                                                         globalIndex const GEOS_UNUSED_PARAM( gridIndex ),
                                                         int const gridFieldIndex,
                                                         real64 const (& gridSurfaceNormal) [3], //Assumed normalized already
                                                         real64 const (& gridSurfacePoint) [3] )
{
  bool interior = false;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    // arrayView2d< globalIndex const > const particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Only particles with cohesive zone flags should need to be checked
      if( particleCohesiveZoneFlag[p] == 1 )
      {
        // Likewise, if the particle does not map to the grid node, then it is outside of its support and does not need to be checked
        // CC: this could be incorrect, CPDI domain scaling should prevent particles from growing larger than cells, but if they do a node completely interior to particle
        // would not integrate the area correctly

        // For no we'll check all the particles, but as soon as a ray is found to be interior we can exit the loop

        // Check if particle maps to current field
        // TODO this should check the particleCohesiveFieldMapping
        int const nodeFlag = 0; // ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNode], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0; // 0 undamaged or "A" field, 1 for "B" field
        localIndex const fieldIndex = nodeFlag * m_numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1

        if( gridFieldIndex == fieldIndex )
        {
          real64 surfaceToParticleCenter[3] = {};
          LvArray::tensorOps::copy< 3 >( surfaceToParticleCenter, gridSurfacePoint );
          LvArray::tensorOps::subtract< 3 >( surfaceToParticleCenter, particlePosition[p] );

          real64 projSurfaceToParticleCenter[3] = {};
          projectToPlane( surfaceToParticleCenter,
                          gridSurfaceNormal,
                          projSurfaceToParticleCenter );

          // Get projected particle RVectors normal to grid surface normal
          real64 projR1[3] = {};
          real64 R1[3] = {};
          LvArray::tensorOps::copy< 3 >( R1, particleRVectors[p][0] );
          projectToPlane( R1,
                          gridSurfaceNormal,
                          projR1 );

          real64 R2[3] = {};
          LvArray::tensorOps::copy< 3 >( R2, particleRVectors[p][1] );
          real64 projR2[3] = {};
          projectToPlane( R2,
                          gridSurfaceNormal,
                          projR2 );

          real64 projR3[3] = {};
          real64 R3[3] = {};
          LvArray::tensorOps::copy< 3 >( R3, particleRVectors[p][2] );
          projectToPlane( R3,
                          gridSurfaceNormal,
                          projR3 );

          real64 R1mag = LvArray::tensorOps::l2Norm< 3 >( projR1 );
          real64 projR1U[3] = {};
           real64 cprojR1 = 0.0;
          if ( R1mag > 1e-16 )
          {
            LvArray::tensorOps::copy< 3 >( projR1U, projR1 );
            LvArray::tensorOps::normalize< 3 >( projR1U );
            cprojR1 = LvArray::tensorOps::AiBi< 3 >( projSurfaceToParticleCenter, projR1U );
          }

          real64 R2mag = LvArray::tensorOps::l2Norm< 3 >( projR2 );
          real64 projR2U[3] = {};
          real64 cprojR2 = 0.0;
          if( R2mag  > 1e-16)
          {
            LvArray::tensorOps::copy< 3 >( projR2U, projR2 );
            LvArray::tensorOps::normalize< 3 >( projR2U );
            cprojR2 = LvArray::tensorOps::AiBi< 3 >( projSurfaceToParticleCenter, projR2U );
          }

          real64 R3mag = LvArray::tensorOps::l2Norm< 3 >( projR3 );
          real64 projR3U[3] = {};
          real64 cprojR3 = 0.0;
          if( R3mag > 1e-16 )
          {
            LvArray::tensorOps::copy< 3 >( projR3U, projR3 );
            LvArray::tensorOps::normalize< 3 >( projR3U );
            cprojR3 = LvArray::tensorOps::AiBi< 3 >( projSurfaceToParticleCenter, projR3U );
          }

          if( R1mag > 1e-16 && R2mag > 1e-16 )
          {
            // GEOS_LOG_RANK( "R1 and R2" );
            if( cprojR1 / R1mag > -1.0 && cprojR1 / R1mag < 1.0 &&
                cprojR2 / R2mag > -1.0 && cprojR2 / R2mag < 1.0 )
            {
              interior = true;
            }
          }

          if( R2mag > 1e-16 && R3mag > 1e-16 )
          {
            //  GEOS_LOG_RANK( "R2 and R3" );
            if( cprojR2 / R2mag > -1.0 && cprojR2 / R2mag < 1.0 &&
                cprojR3 / R3mag > -1.0 && cprojR3 / R3mag < 1.0 )
            {
              interior = true;
            }
          }

          if( R1mag > 1e-16 && R3mag > 1e-16 )
          {
            // GEOS_LOG_RANK( "R1 and R3" );
            if( cprojR1 / R1mag > -1.0 && cprojR1 / R1mag < 1.0 &&
                cprojR3 / R3mag > -1.0 && cprojR3 / R3mag < 1.0 )
            {
              interior = true;
            }
          }

        }

      }
    }

  } );

  return interior;
}

// Currently unused but might be useful in the future
void SolidMechanicsMPM::projectToPlane( real64 const (& vector)[3],
                                        real64 const (& normal)[3],
                                        real64 (& projection)[3] )
{
  real64 normalComponent = LvArray::tensorOps::AiBi< 3 >( vector, normal );

  real64 normalVector[3] = {};
  LvArray::tensorOps::copy< 3 >( normalVector, normal );
  LvArray::tensorOps::scale< 3 >( normalVector, normalComponent );

  LvArray::tensorOps::copy< 3 >( projection, vector );
  LvArray::tensorOps::subtract< 3 >( projection, normalVector );
}

// Used for mapping displacement to cohesive grid nodes
void SolidMechanicsMPM::computeDistanceToParticleSurface( real64 (& normal)[3],
                                                          arraySlice2d< real64 const > const rVectors,
                                                          real64 distanceToSurface )
{
  LvArray::tensorOps::normalize< 3 >( normal );

  real64 r1[3] = {};
  LvArray::tensorOps::copy< 3 >(  r1, rVectors[0] );
  real64 r2[3] = {};
  LvArray::tensorOps::copy< 3 >(  r2, rVectors[1] );
  real64 r3[3] = {};
  LvArray::tensorOps::copy< 3 >(  r3, rVectors[2] );

  // Compute particle face surface normals from RVectors
  real64 s1[3] = {};
  real64 s2[3] = {};
  real64 s3[3] = {};
  real64 s4[3] = {};
  real64 s5[3] = {};
  real64 s6[3] = {};

  LvArray::tensorOps::crossProduct( s3, r1, r2 );
  LvArray::tensorOps::crossProduct( s1, r2, r3 );
  LvArray::tensorOps::crossProduct( s2, r3, r1 );
  LvArray::tensorOps::crossProduct( s6, r2, r1 );
  LvArray::tensorOps::crossProduct( s4, r3, r2 );
  LvArray::tensorOps::crossProduct( s5, r1, r3 );

  LvArray::tensorOps::normalize< 3 >( s1 );
  LvArray::tensorOps::normalize< 3 >( s2 );
  LvArray::tensorOps::normalize< 3 >( s3 );
  LvArray::tensorOps::normalize< 3 >( s4 );
  LvArray::tensorOps::normalize< 3 >( s5 );
  LvArray::tensorOps::normalize< 3 >( s6 );

  real64 dS1 = LvArray::tensorOps::AiBi< 3 >( s1, r1);
  real64 dS2 = LvArray::tensorOps::AiBi< 3 >( s2, r2);
  real64 dS3 = LvArray::tensorOps::AiBi< 3 >( s3, r3);
  real64 dS4 = LvArray::tensorOps::AiBi< 3 >( s4, r1);
  real64 dS5 = LvArray::tensorOps::AiBi< 3 >( s5, r2);
  real64 dS6 = LvArray::tensorOps::AiBi< 3 >( s6, r3);

  real64 dN1 = LvArray::tensorOps::AiBi< 3 >( normal, s1 );
  real64 dN2 = LvArray::tensorOps::AiBi< 3 >( normal, s2 );
  real64 dN3 = LvArray::tensorOps::AiBi< 3 >( normal, s3 );
  real64 dN4 = LvArray::tensorOps::AiBi< 3 >( normal, s4 );
  real64 dN5 = LvArray::tensorOps::AiBi< 3 >( normal, s5 );
  real64 dN6 = LvArray::tensorOps::AiBi< 3 >( normal, s6 );

  real64 tolerance = 1e-16;

  distanceToSurface = DBL_MAX;
  if( abs( dN1 ) > tolerance && dS1 / dN1 > 0 )
  {
    distanceToSurface = fmin( distanceToSurface, dS1 / dN1 );
  }

  if( abs( dN2 ) > tolerance && dS2 / dN2 > 0 )
  {
    distanceToSurface = fmin( distanceToSurface, dS2 / dN2 );
  }

  if( abs( dN3 ) > tolerance && dS3 / dN3 > 0 )
  {
    distanceToSurface = fmin( distanceToSurface, dS3 / dN3 );
  }

  if( abs( dN4 ) > tolerance && dS4 / dN4 > 0 )
  {
    distanceToSurface = fmin( distanceToSurface, dS4 / dN4 );
  }

  if( abs( dN5 ) > tolerance && dS5 / dN5 > 0 )
  {
    distanceToSurface = fmin( distanceToSurface, dS5 / dN5 );
  }

  if( abs( dN6 ) > tolerance && dS6 / dN6 > 0 )
  {
    distanceToSurface = fmin( distanceToSurface, dS6 / dN6 );
  }

}


void SolidMechanicsMPM::enforceCohesiveLaw( ParticleManager & particleManager,
                                            NodeManager & nodeManager )
{
  int const numDims = m_numDims;
  int numCohesiveNodes = m_cohesiveNodeGlobalIndices.size();

  array2d< real64 > tempGridMassLocal( numCohesiveNodes, m_numVelocityFields );
  // array2d< real64 > tempGridVolumeLocal( numCohesiveNodes, m_numVelocityFields );
  // array3d< real64 > tempGridCenterOfVolumeLocal( numCohesiveNodes, m_numVelocityFields, 3 );
  array3d< real64 > tempGridDisplacementLocal( numCohesiveNodes, m_numVelocityFields, 3 );
  array3d< real64 > tempGridParticleSurfaceNormalLocal( numCohesiveNodes, m_numVelocityFields, 3 );
  array4d< real64 > tempGridDeformationGradientCofactorLocal( numCohesiveNodes, m_numVelocityFields, 3, 3 );

  array3d< real64 > tempGridCohesiveTraction( numCohesiveNodes, m_numVelocityFields, 3 ); // Should zero this first

  // Initialize temporary grid fields to zero
  for( int g  = 0; g < numCohesiveNodes; ++g)
  {
    for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex )
    {
      tempGridMassLocal[g][fieldIndex] = 0.0;
      // tempGridVolumeLocal[g][fieldIndex] = 0.0;

      for( int i = 0; i < m_numDims; ++i)
      {
        tempGridDisplacementLocal[g][fieldIndex][i] = 0.0;
        // tempGridCenterOfVolumeLocal[g][fieldIndex][i] = 0.0;
        tempGridParticleSurfaceNormalLocal[g][fieldIndex][i] = 0.0;

        tempGridCohesiveTraction[g][fieldIndex][i] = 0.0;

        for( int j = 0; j < m_numDims; ++j)
        {
          tempGridDeformationGradientCofactorLocal[g][fieldIndex][i][j] = 0.0;
        }
      }
    }
  }

  // Map relevant particle fields to grid
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
    arrayView1d< int > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

    arrayView3d< real64 const > const particleReferenceRVectors = subRegion.getField< fields::mpm::particleReferenceRVectors >();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();

    arrayView2d< globalIndex const > const particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();
    arrayView2d< real64 const > const particleReferenceShapeFunctionValues = subRegion.getField< fields::mpm::particleReferenceShapeFunctionValues >();
    arrayView2d< int const > const particleCohesiveFieldMapping = subRegion.getField< fields::mpm::particleCohesiveFieldMapping >();

    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=, &tempGridMassLocal, &tempGridDisplacementLocal, &tempGridParticleSurfaceNormalLocal, &tempGridDeformationGradientCofactorLocal] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      if( particleCohesiveZoneFlag[p] == 1 )
      {
        particleCohesiveZoneFlag[p] = 0;
        if( particleDamage[p] < 1.0 )
        {
          // Map particle displacement to grid
          for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
          {
            globalIndex const mappedNode = particleReferenceMappedNodes[p][g];

            if( m_cohesiveNodeGlobalIndices.contains( mappedNode ) )
            {
              // CC: TODO must be a better way to find index in temp arrays
              localIndex nodeIndex = 0;
              for( localIndex n = 0; n < numCohesiveNodes; ++n )
              {
                if( m_cohesiveNodeGlobalIndices[n] == mappedNode )
                {
                  nodeIndex = n;
                  break;
                }
              }

              localIndex const fieldIndex = particleCohesiveFieldMapping[p][g];
              
              if( m_enableCohesiveFailure == 0 || m_cohesiveGridNodeDamages[nodeIndex][fieldIndex] < 1.0 )
              {
                particleCohesiveZoneFlag[p] = 1; // Reenable particle cohesive flag if any of the cohesive nodes are undamaged
              }

              real64 shapeFunctionValue = particleReferenceShapeFunctionValues[p][g];

              real64 deformationGradient[3][3] = {};
              LvArray::tensorOps::copy< 3, 3 >( deformationGradient, particleDeformationGradient[p] );

              real64 initialDistanceToSurface = 0.0;
              real64 referenceSurfaceNormal[3] = {};
              LvArray::tensorOps::copy< 3 >( referenceSurfaceNormal, particleReferenceSurfaceNormal[p] );
              computeDistanceToParticleSurface( referenceSurfaceNormal,
                                                particleReferenceRVectors[p],
                                                initialDistanceToSurface );

              real64 initialSurfacePoint[3] = {};
              LvArray::tensorOps::copy< 3 >( initialSurfacePoint, particleReferenceSurfaceNormal[p] );
              LvArray::tensorOps::scale< 3 >( initialSurfacePoint,  initialDistanceToSurface );

              real64 deformedSurfacePoint[3] = {};
              LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( deformedSurfacePoint, particleDeformationGradient[p], initialSurfacePoint );

              real64 deformationGradientCofactor[3][3] = {};
              LvArray::tensorOps::cofactor< 3 >( deformationGradientCofactor, deformationGradient );

              tempGridMassLocal[nodeIndex][fieldIndex] += particleMass[p] * shapeFunctionValue;
              // tempGridVolumeLocal[nodeIndex][fieldIndex] += particleVolume[p] * shapeFunctionValue;

              for( int i  = 0; i < numDims; ++i )
              {
                tempGridDisplacementLocal[nodeIndex][fieldIndex][i] += particleMass[p] * ( particlePosition[p][i] - particleReferencePosition[p][i] + deformedSurfacePoint[i] - initialSurfacePoint[i] ) * shapeFunctionValue;
                // tempGridCenterOfVolumeLocal[nodeIndex][fieldIndex][i] += particleVolume[p] * (particlePosition[p][i] - m_referenceCohesiveGridNodePositions[nodeIndex][i])* shapeFunctionValue;
                tempGridParticleSurfaceNormalLocal[nodeIndex][fieldIndex][i] += particleMass[p] * particleSurfaceNormal[p][i] * shapeFunctionValue;
              }

              for( int i = 0; i < 3; ++i )
              {
                for( int j = 0; j < 3; ++j )
                {
                  tempGridDeformationGradientCofactorLocal[nodeIndex][fieldIndex][i][j] += particleMass[p] * shapeFunctionValue * deformationGradientCofactor[i][j];
                }
              }
            }
          }
        }

        // If flag was disabled, update particle surface flag to 4 from 3
        if( particleCohesiveZoneFlag[p] == 0 )
        {
          particleSurfaceFlag[p] = 4;
        }
      }
    });
  });

  // Sync temporary grid fields
  array2d< real64 > tempGridMassGlobal( numCohesiveNodes, m_numVelocityFields );
  // array2d< real64 > tempGridVolumeGlobal( numCohesiveNodes, m_numVelocityFields );
  // array3d< real64 > tempGridCenterOfVolumeGlobal( numCohesiveNodes, m_numVelocityFields, 3 );
  array3d< real64 > tempGridDisplacementGlobal( numCohesiveNodes, m_numVelocityFields, 3 );
  array3d< real64 > tempGridParticleSurfaceNormalGlobal( numCohesiveNodes, m_numVelocityFields, 3 );
  array4d< real64 > tempGridDeformationGradientCofactorGlobal( numCohesiveNodes, m_numVelocityFields, 3, 3 );

  MpiWrapper::allReduce( tempGridMassLocal,
                         tempGridMassGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // MpiWrapper::allReduce( tempGridVolumeLocal,
  //                        tempGridVolumeGlobal,
  //                        MpiWrapper::Reduction::Sum,
  //                        MPI_COMM_GEOS );

  // MpiWrapper::allReduce( tempGridCenterOfVolumeLocal,
  //                        tempGridCenterOfVolumeGlobal,
  //                        MpiWrapper::Reduction::Sum,
  //                        MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridDisplacementLocal,
                         tempGridDisplacementGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridParticleSurfaceNormalLocal,
                         tempGridParticleSurfaceNormalGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridDeformationGradientCofactorLocal,
                         tempGridDeformationGradientCofactorGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // Normalize fields after global sync
  forAll< serialPolicy >( numCohesiveNodes, [=, &tempGridDisplacementGlobal, &tempGridParticleSurfaceNormalGlobal, &tempGridDeformationGradientCofactorGlobal] GEOS_HOST ( localIndex const g )
  {
    for(localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
    {
      if ( tempGridMassGlobal[g][fieldIndex] > m_smallMass )
      {
        for( int i = 0; i < m_numDims; ++i ){
          // tempGridCenterOfVolumeGlobal[g][fieldIndex][i] /= tempGridVolumeGlobal[g][fieldIndex];
          tempGridDisplacementGlobal[g][fieldIndex][i] /= tempGridMassGlobal[g][fieldIndex];
          tempGridParticleSurfaceNormalGlobal[g][fieldIndex][i] /= tempGridMassGlobal[g][fieldIndex];

          for( int j = 0; j < m_numDims; ++j)
          {
            tempGridDeformationGradientCofactorGlobal[g][fieldIndex][i][j] /=tempGridMassGlobal[g][fieldIndex];
          }
        }
      }
    }
  } );

  // Compute traction law on grid
  forAll< serialPolicy >( numCohesiveNodes, [&, tempGridMassGlobal, tempGridDisplacementGlobal, tempGridParticleSurfaceNormalGlobal, tempGridDeformationGradientCofactorGlobal] GEOS_HOST ( localIndex const gg )
  {
    for( localIndex A = 0; A < m_numVelocityFields - 1; A++ )
    {
      for( localIndex B = A + 1; B < m_numVelocityFields; B++ )
      {
        bool active = ( tempGridMassGlobal[gg][A] > m_smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( tempGridParticleSurfaceNormalGlobal[gg][A] ) > 1.0e-16 )
                      and
                      ( tempGridMassGlobal[gg][B] > m_smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( tempGridParticleSurfaceNormalGlobal[gg][B] ) > 1.0e-16 );


        if( active )
        {
          real64 initialAreaVectorA[3] = {};
          LvArray::tensorOps::copy< 3 >( initialAreaVectorA, m_referenceCohesiveGridNodeSurfaceNormals[gg][A] );

          real64 areaA = m_referenceCohesiveGridNodeAreas[gg][A];
          LvArray::tensorOps::scale< 3 >( initialAreaVectorA, areaA );

          if(m_numDims < 3 )
          {
            tempGridDeformationGradientCofactorGlobal[gg][A][2][2] = 1.0;
            tempGridDeformationGradientCofactorGlobal[gg][B][2][2] = 1.0;
          }

          real64 areaVectorA[3] = {};
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( areaVectorA, tempGridDeformationGradientCofactorGlobal[gg][A], initialAreaVectorA );

          real64 initialAreaVectorB[3] = {};
          LvArray::tensorOps::copy< 3 >( initialAreaVectorB, m_referenceCohesiveGridNodeSurfaceNormals[gg][B] );

          real64 areaB = m_referenceCohesiveGridNodeAreas[gg][B];
          LvArray::tensorOps::scale< 3 >( initialAreaVectorB, areaB );

          real64 areaVectorB[3] = {};
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( areaVectorB, tempGridDeformationGradientCofactorGlobal[gg][B], initialAreaVectorB );

          computeCohesiveTraction( gg,
                                   A,
                                   B,
                                   tempGridMassGlobal[gg][A],
                                   tempGridMassGlobal[gg][B],
                                   tempGridDisplacementGlobal[gg][A],
                                   tempGridDisplacementGlobal[gg][B],
                                   areaVectorA,
                                   areaVectorB,
                                   tempGridParticleSurfaceNormalGlobal[gg][A],
                                   tempGridParticleSurfaceNormalGlobal[gg][B],
                                   tempGridCohesiveTraction[gg][A],
                                   tempGridCohesiveTraction[gg][B] );
        }
        else
        {
          m_cohesiveGridNodeDamages[gg][A] = 1.0;
          m_cohesiveGridNodeDamages[gg][B] = 1.0;
        }
      }
    }
  } );

  // // For debugging print to console
  // forAll< serialPolicy >( numCohesiveNodes, [=, &tempGridDisplacementGlobal, &tempGridParticleSurfaceNormalGlobal] GEOS_HOST ( localIndex const g )
  // {
  //   for(localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
  //   {
  //    GEOS_LOG_RANK_0("g_cz: " << g << ", " << 
  //                    "fieldIndex: " << fieldIndex << ", " << 
  //                    "u: {" << tempGridDisplacementGlobal[g][fieldIndex][0] << ", " << tempGridDisplacementGlobal[g][fieldIndex][1] << ", " << tempGridDisplacementGlobal[g][fieldIndex][2] << "}, " << 
  //                    "area0: " << m_referenceCohesiveGridNodeAreas[g][fieldIndex] << ", " << 
  //                    "normal: {" << tempGridParticleSurfaceNormalGlobal[g][fieldIndex][0] << ", " << tempGridParticleSurfaceNormalGlobal[g][fieldIndex][1] << ", " << tempGridParticleSurfaceNormalGlobal[g][fieldIndex][2] << "}");
  //   }
  // } );

  // // Sync cohesive damage ( check if this is the reason some particles don't lose there cohesive flag after cohesive node is damaged)
  // MpiWrapper::allReduce( m_cohesiveGridNodeDamages.data(),
  //                        m_cohesiveGridNodeDamages.data(),
  //                        m_cohesiveGridNodeDamages.size(),
  //                        MpiWrapper::getMpiOp( MpiWrapper::Reduction::Max ),
  //                        MPI_COMM_GEOS );

  // for( int n = 0; n < numCohesiveNodes; ++n)
  // {
  //   for( int f = 0; f < m_numVelocityFields; f++ )
  //   {
  //     // GEOS_LOG_RANK_0( "n: " << n << ", f: " << f << ", c_dmg: " << m_cohesiveGridNodeDamages[n][f] << ", dn_max: " << << ", dt_max: " << );
  //   }
  // }

  // CC: debug, for debugging temp grid variables to visualize in paraview
  // arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView1d< int > const gridCohesiveNode = nodeManager.getReference< array1d< int > >( viewKeyStruct::gridCohesiveNodeString() );
  arrayView3d< real64 > const gridDisplacement = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDisplacementString() );
  // arrayView3d< real64 >  const gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );
  // arrayView3d< real64 > const gridParticleSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridParticleMappedSurfaceNormalStringString() );
  arrayView3d< real64 > const gridCohesiveArea = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCohesiveAreaString() );
  arrayView3d< real64 > const gridCohesiveForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCohesiveForceString() );

  arrayView1d< globalIndex > localToGlobalMap = nodeManager.localToGlobalMap();
  forAll< serialPolicy >( nodeManager.size(), [=] GEOS_HOST ( localIndex const g )
  {
    for( int n = 0; n < numCohesiveNodes; ++n)
      if( localToGlobalMap[g] == m_cohesiveNodeGlobalIndices[n] )
      {
        gridCohesiveNode[g] = 1;
        for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex)
        {
          for(int i = 0; i < m_numDims; ++i)
          {
            // gridCenterOfVolume[g][fieldIndex][i] = tempGridCenterOfVolumeGlobal[n][fieldIndex][i];
            // gridParticleSurfaceNormal[g][fieldIndex][i] = tempGridParticleSurfaceNormalGlobal[n][fieldIndex][i];
            gridCohesiveArea[g][fieldIndex][i] = m_referenceCohesiveGridNodeSurfaceNormals[n][fieldIndex][i];
            gridDisplacement[g][fieldIndex][i] = tempGridDisplacementGlobal[n][fieldIndex][i];
            gridCohesiveForce[g][fieldIndex][i] = tempGridCohesiveTraction[n][fieldIndex][i];
          }
        }
      }
  } );

  // Map cohesive law back to particle
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView1d< int > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    arrayView2d< globalIndex const > const particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();
    arrayView2d< real64 const > const particleReferenceShapeFunctionValues = subRegion.getField< fields::mpm::particleReferenceShapeFunctionValues >();
    arrayView2d< int const > const particleCohesiveFieldMapping = subRegion.getField< fields::mpm::particleCohesiveFieldMapping >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Zero particle fields for fresh mapping
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( int i=0; i< m_numDims; ++i )
      {
        particleCohesiveForce[p][i] = 0.0;
      }
    } );

    // Map to particles
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      if( particleCohesiveZoneFlag[p] == 1 )
      {
        for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
        {
          globalIndex const mappedNode = particleReferenceMappedNodes[p][g];

          if( m_cohesiveNodeGlobalIndices.contains( mappedNode ) )
          {

            real64 shapeFunctionValue = particleReferenceShapeFunctionValues[p][g];

            // CC: TODO must be a better way to find index in temp arrays
            localIndex nodeIndex = 0;
            for( localIndex n = 0; n < numCohesiveNodes; ++n )
            {
              if( m_cohesiveNodeGlobalIndices[n] == mappedNode )
              {
                nodeIndex = n;
                break;
              }
            }

            localIndex const fieldIndex = particleCohesiveFieldMapping[p][g];

            // TODO: if a particle cohesive zone flags are turned off because they become fully damaged and no undamaged particles are still mapping to a cohesive grid node,
            // Should that node be marke as fully damaged?
            if( tempGridMassGlobal[nodeIndex][fieldIndex] > m_smallMass) // 1e-20 )
            {
              if( m_enableCohesiveFailure == 0 || m_cohesiveGridNodeDamages[nodeIndex][fieldIndex] < 1.0 )
              {
                for( localIndex i = 0; i < m_numDims; ++i )
                {
                  particleCohesiveForce[p][i] += tempGridCohesiveTraction[nodeIndex][fieldIndex][i] * particleMass[p] * shapeFunctionValue / tempGridMassGlobal[nodeIndex][fieldIndex];
                }
              }
              else
              {
                // If cohesive node is fully damaged then any particles mapping to it should have their cohesive zone flags turned off.
                particleCohesiveZoneFlag[p] = 0;
              }
            }
            else
            {
              // If cohesive node is fully damaged then any particles mapping to it should have their cohesive zone flags turned off.
              particleCohesiveZoneFlag[p] = 0;
            }
          }
        }
      }

    } );
  } );
}

void SolidMechanicsMPM::computeCohesiveTraction( localIndex g,
                                                 localIndex A,
                                                 localIndex B,
                                                 real64 mA,
                                                 real64 mB,
                                                 arraySlice1d< real64 const > const dA,
                                                 arraySlice1d< real64 const > const dB,
                                                 real64 const (& sA )[3], // rename this to something other than s, to be consistent with other nomenclature
                                                 real64 const (& sB )[3], // Particle surface area along surface normal
                                                 arraySlice1d< real64 const > const nA,
                                                 arraySlice1d< real64 const > const nB,
                                                 arraySlice1d< real64 > const tA,
                                                 arraySlice1d< real64 > const tB )
{
  LvArray::tensorOps::fill< 3 >( tA, 0.0);
  LvArray::tensorOps::fill< 3 >( tB, 0.0);

  // Total mass for the contact pair.
  real64 mAB = mA + mB;

  // Outward normal of field A with respect to field B.
  real64 nAB[3];

  // Mass-weighted average of the field normals
  for( localIndex i = 0; i < 3; ++i )
  {
    nAB[i] = nA[i] * mA - nB[i] * mB;
  }

  // Normalize the effective surface normal
  if( m_planeStrain == 1 )
  {
    nAB[2] = 0.0;
  }

  real64 norm = LvArray::tensorOps::l2Norm< 3 >( nAB );

  // If normal magnitude is zero for any reason just skip (e.g. no traction from cohesive law)
  if ( norm < 1e-20 )
  {
    return;
  }

  nAB[0] /= norm;
  nAB[1] /= norm;
  nAB[2] /= norm;

  // Flip normal because positive displacement is away from interface
  nAB[0] *= -1.0;
  nAB[1] *= -1.0;
  nAB[2] *= -1.0;

  // Compute tangent direction to normal
  real64 dANorm = LvArray::tensorOps::l2Norm< 3 >( dA );
  real64 normalDisplacementAVec[3] = {};
  real64 tangentialDisplacementAVec[3] = {};
  real64 normalDisplacementA = 0.0;
  if( dANorm > 1e-20 )
  {
    normalDisplacementA = LvArray::tensorOps::AiBi< 3 >( nAB, dA );
    LvArray::tensorOps::copy< 3 >( normalDisplacementAVec, nAB );
    LvArray::tensorOps::scale< 3 >( normalDisplacementAVec, normalDisplacementA );

    LvArray::tensorOps::copy< 3 >( tangentialDisplacementAVec, dA );
    LvArray::tensorOps::subtract< 3 >( tangentialDisplacementAVec, normalDisplacementAVec );
  }

  // Invert surface normal for field B
  nAB[0] *= -1.0;
  nAB[1] *= -1.0;
  nAB[2] *= -1.0;

  real64 dBNorm = LvArray::tensorOps::l2Norm< 3 >( dB );
  real64 normalDisplacementBVec[3] = {};
  real64 tangentialDisplacementBVec[3] = {};
  real64 normalDisplacementB = 0.0;
  if( dBNorm  > 1e-20 )
  {
    normalDisplacementB = LvArray::tensorOps::AiBi< 3 >( nAB, dB );
    LvArray::tensorOps::copy< 3 >( normalDisplacementBVec, nAB );
    LvArray::tensorOps::scale< 3 >( normalDisplacementBVec, normalDisplacementB);

    LvArray::tensorOps::copy< 3 >( tangentialDisplacementBVec, dB);
    LvArray::tensorOps::subtract< 3 >( tangentialDisplacementBVec, normalDisplacementBVec );
  }

  // Revert surface normal
  nAB[0] *= -1.0;
  nAB[1] *= -1.0;
  nAB[2] *= -1.0;

  real64 totalNormalDisplacement = normalDisplacementA + normalDisplacementB;
  m_maxCohesiveGridNodeNormalDisplacement[g][A][B] = fmax( m_maxCohesiveGridNodeNormalDisplacement[g][A][B], totalNormalDisplacement);
  m_maxCohesiveGridNodeNormalDisplacement[g][B][A] = fmax( m_maxCohesiveGridNodeNormalDisplacement[g][B][A], totalNormalDisplacement);

  real64 tangentialInterfaceDisplacement[3]  = {};
  LvArray::tensorOps::copy< 3 >( tangentialInterfaceDisplacement, tangentialDisplacementBVec );
  LvArray::tensorOps::subtract< 3 >( tangentialInterfaceDisplacement, tangentialDisplacementAVec );

  real64 totalTangentialDisplacement = LvArray::tensorOps::l2Norm< 3 >( tangentialInterfaceDisplacement );
  m_maxCohesiveGridNodeTangentialDisplacement[g][A][B] = fmax( m_maxCohesiveGridNodeTangentialDisplacement[g][A][B], totalTangentialDisplacement);
  m_maxCohesiveGridNodeTangentialDisplacement[g][B][A] = fmax( m_maxCohesiveGridNodeTangentialDisplacement[g][B][A], totalTangentialDisplacement);

  real64 normalStress = 0.0;
  real64 shearStress = 0.0;
  real64 dmg = fmax( m_cohesiveGridNodeDamages[g][A], m_cohesiveGridNodeDamages[g][B]);
  switch( m_cohesiveLaw )
  {
    case CohesiveLawOption::Uncoupled:
      uncoupledCohesiveLaw( totalNormalDisplacement,
                            totalTangentialDisplacement,
                            normalStress,
                            shearStress,
                            dmg );
      break;
    case CohesiveLawOption::NeedlemanXu:
      needlemanXuCohesiveLaw( totalNormalDisplacement,
                              totalTangentialDisplacement,
                              normalStress,
                              shearStress,
                              dmg );
      break;
    case CohesiveLawOption::Polymer:
      polymerCohesiveLaw( totalNormalDisplacement,
                          totalTangentialDisplacement,
                          normalStress,
                          shearStress,
                          dmg );
      break;
    default:
      GEOS_ERROR( "No cohesive law of the kind specified!" );
      break;
  }

  // Mass-weighted average of the projected area
  // Do we want to add choice of weighting as we do in contact calculations for normals?
  real64 areaA[3] = {};
  LvArray::tensorOps::copy< 3 >( areaA, sA );
  LvArray::tensorOps::scale< 3 >( areaA, mA );

  real64 areaB[3] = {};
  LvArray::tensorOps::copy< 3 >( areaB, sB );
  LvArray::tensorOps::scale< 3 >( areaB, -mB );

  LvArray::tensorOps::add< 3 >( areaA, areaB );
  LvArray::tensorOps::scale< 3 >( areaA, 1 / mAB );

  real64 surfaceArea = LvArray::tensorOps::AiBi< 3 >( nAB, areaA );

  real64 normalForce = normalStress * surfaceArea;
  real64 shearForce = shearStress * surfaceArea;

  if ( m_preventCZInterpentration == 1 && normalForce < 0 )
  {
    normalForce = 0.0;
  }

  // If cohesive failure is enabled only compute tractions if cohesive node has not failed
  // (e.g. max normal and tangential displacements have not been surpassed)
  // Mark node as failed for particle mapping
  if( m_enableCohesiveFailure == 1)
  {
    m_cohesiveGridNodeDamages[g][A] = fmax(m_cohesiveGridNodeDamages[g][A], dmg);
    m_cohesiveGridNodeDamages[g][B] = fmax(m_cohesiveGridNodeDamages[g][B], dmg);
  }
  else
  {
    m_cohesiveGridNodeDamages[g][A] = 0.0;
    m_cohesiveGridNodeDamages[g][B] = 0.0;
  }

  for( localIndex i = 0; i < m_numDims; ++i )
  {
    tA[i] += normalForce * -nAB[i] * (1.0 - m_cohesiveGridNodeDamages[g][A]);
    tB[i] += normalForce * nAB[i] * (1.0 - m_cohesiveGridNodeDamages[g][B]);
  }

  if ( fabs(totalTangentialDisplacement) > 1e-20 )
  {
    real64 tAB[3] = {}; // Tangent unit vector
    LvArray::tensorOps::copy< 3 >( tAB, tangentialInterfaceDisplacement );
    LvArray::tensorOps::scale< 3 >( tAB, 1 / totalTangentialDisplacement );

    for( localIndex i = 0; i < m_numDims; ++i)
    {
      tA[i] += shearForce * tAB[i] * (1.0 - m_cohesiveGridNodeDamages[g][A]);
      tB[i] += shearForce * -tAB[i] * (1.0 - m_cohesiveGridNodeDamages[g][B]);
    }
  }
}

void SolidMechanicsMPM::uncoupledCohesiveLaw( real64 normalDisplacement,
                                              real64 tangentialDisplacement,
                                              real64 & normalStress,
                                              real64 & shearStress,
                                              real64 & damage )
{
  if( ( normalDisplacement >= m_maxCohesiveNormalDisplacement ) ||
      ( tangentialDisplacement >= m_maxCohesiveTangentialDisplacement ) )
  {
    damage = fmax(1.0, damage);
  }

  normalStress = -m_normalForceConstant * normalDisplacement;
  shearStress = -m_shearForceConstant * tangentialDisplacement;
}

void SolidMechanicsMPM::needlemanXuCohesiveLaw( real64 normalDisplacement,
                                                real64 tangentialDisplacement,
                                                real64 & normalStress,
                                                real64 & shearStress,
                                                real64 & damage )
{
  if( ( normalDisplacement >= m_maxCohesiveNormalDisplacement ) ||
      ( tangentialDisplacement >= m_maxCohesiveTangentialDisplacement ) )
  {
    damage = fmax(1.0, damage);
  }

  real64 e = 2.7182818284590452353602874713526624977572470936999; // Euler's number
  real64 r = 0.0;
  real64 q = 1.0;
  real64 maxNormalStress = m_maxCohesiveNormalStress;
  real64 maxShearStress = m_maxCohesiveShearStress;
  real64 normalCharacteristicDisplacement = m_characteristicNormalDisplacement;
  real64 tangentialCharacteristicDisplacement = m_characteristicTangentialDisplacement;
  real64 normalWorkOfSeparation = e * maxNormalStress * normalCharacteristicDisplacement;
  real64 shearWorkOfSeparation = sqrt( e / 2 ) * maxShearStress * tangentialCharacteristicDisplacement;

  real64 normalizedNormalDisplacement = normalDisplacement / normalCharacteristicDisplacement;
  real64 normalizedTangentialDisplacement = tangentialDisplacement/ tangentialCharacteristicDisplacement;

  normalStress = -( normalWorkOfSeparation / normalCharacteristicDisplacement ) * exp( -normalizedNormalDisplacement ) *
                        ( normalizedNormalDisplacement * exp( -pow( normalizedTangentialDisplacement, 2) ) + ( 1.0 - q ) / ( r - 1.0 ) * ( 1.0 - exp( -pow( normalizedTangentialDisplacement, 2 ) ) ) * ( r - normalizedNormalDisplacement ) );
  shearStress = -( shearWorkOfSeparation / tangentialCharacteristicDisplacement ) * ( 2 * normalCharacteristicDisplacement / tangentialCharacteristicDisplacement ) * normalizedTangentialDisplacement *
                        ( q + ( r - q ) / ( r - 1.0 ) * normalizedNormalDisplacement ) *
                        exp( -normalizedNormalDisplacement ) * exp( -pow( normalizedTangentialDisplacement, 2) );
}

void SolidMechanicsMPM::polymerCohesiveLaw( real64 normalDisplacement,
                                            real64 tangentialDisplacement,
                                            real64 & normalStress,
                                            real64 & shearStress,
                                            real64 & damage )
{
  real64 t = m_polymerCZThickness; // m_characteristicNormalDisplacement;
  real64 tSqr = t * t;

  real64 normalDisplacementSqr = normalDisplacement * normalDisplacement;
  real64 shearDisplacementSqr = tangentialDisplacement * tangentialDisplacement;

  // These should eventually be user inputs
  real64 k = m_polymerCZBulkModulus; // 0.8; // Polymer bulk modulus
  real64 g = m_polymerCZShearModulus; // 0.2; // Polymer shear modulus
  real64 gSqr = g * g;
  real64 yieldStrength0 = m_polymerCZYieldStrength0; // 0.003;
  real64 r0 = m_polymerCZR0; // 0.0225;
  real64 r1 = m_polymerCZR1; // 0.0108;
  real64 r2 = m_polymerCZR2; // 0.177;
  real64 Gr = m_polymerCZGr; // 0.00174;
  real64 lambdaMax = m_polymerCZMaxStretch; // 2.14; // Polymer max stretch

  real64 lambda = sqrt(pow(normalDisplacement + t,2) + shearDisplacementSqr)/t; // Stretch

  if( lambda > lambdaMax )
  {
    damage = 1.0;
  }

  real64 sigma_H = Gr*pow(lambda * lambda - 1/lambda,2);

  real64 tau = sqrt(gSqr*(4*t*normalDisplacement*(normalDisplacementSqr + shearDisplacementSqr)+pow(normalDisplacementSqr+shearDisplacementSqr,2)+tSqr * (4*normalDisplacementSqr+3*shearDisplacementSqr))/(pow(t*(t+normalDisplacement),2)));

  real64 gamma_p = fmax(0.0, (tau-(yieldStrength0+sigma_H))/(k + (4.0/3.0)*g) );

  real64 R_gamma = r0*exp(-pow(gamma_p/r1, r2));

  real64 yieldStrength = yieldStrength0 + R_gamma + sigma_H;

  real64 scale  = 1.0;
  if( tau > yieldStrength )
  {
    scale = yieldStrength / tau;
  }

  normalStress = -scale * (normalDisplacement*(k*(t+normalDisplacement)+g*(2*t+normalDisplacement)))/(t*(t+normalDisplacement));
  shearStress = -g * tangentialDisplacement / t;
}

void SolidMechanicsMPM::updateGridBoreholeStress( NodeManager & nodeManager )
{
  real64 boreholeRadiusSqr = m_boreholeRadius * m_boreholeRadius;
  real64 boreholeStress[6];
  LvArray::tensorOps::copy< 6 >( boreholeStress, m_boreholeStress );

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  forAll< serialPolicy >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    if ( LvArray::math::square( gridPosition[g][0] ) + LvArray::math::square( gridPosition[g][1] ) < boreholeRadiusSqr )
    {
      LvArray::tensorOps::copy< 6 >( gridBoreholeStress[g], boreholeStress );
    }
  } );
}

void SolidMechanicsMPM::particleToGrid( real64 const time_n,
                                        integer const cycleNumber,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // On-the-fly shape function computations
  int const hasContact = m_hasContact;
  localIndex const numDims = m_numDims;
  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const useArtificialViscosity = m_useArtificialViscosity;
  int const computeXProfile = m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) );
  int const explicitSurfaceNormalInfluence = m_explicitSurfaceNormalInfluence;
  int const enableSurfaceTension = m_enableSurfaceTension;

  localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView2d< real64 > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );

  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView2d< int > const & gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );

  arrayView2d< real64 const > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  arrayView3d< real64 > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView2d< real64 > const gridSurfaceFieldMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridSurfaceFieldMassString() );
  arrayView3d< real64 > const gridSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfacePositionString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > particleArtificialViscosity;
    if( useArtificialViscosity == 1 )
    {
      particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    }
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

    arrayView1d< real64 const > particleSurfaceCurvature;
    if( enableSurfaceTension == 1 )
    {
      particleSurfaceCurvature = subRegion.getField< fields::mpm::particleSurfaceCurvature >();
    }

    // Get views to mapping arrays
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    ParticleType const particleType = subRegion.getParticleType();

    #ifndef GEOS_USE_DEVICE
      arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
      arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
      arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];
      GEOS_UNUSED_VAR( particleRVectors );
      GEOS_UNUSED_VAR( particleType );
      GEOS_UNUSED_VAR( ijkMap );
    #endif

    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // On-the-fly shape function computation
      #ifdef GEOS_USE_DEVICE
        localIndex mappedNodes[64];
        real64 shapeFunctionValues[64];
        real64 shapeFunctionGradientValues[64][3];
        mapNodesAndComputeShapeFunctions( ijkMap,
                                          xLocalMin,
                                          hEl,
                                          particleType,
                                          particlePosition[p],
                                          particleRVectors[p],
                                          gridPosition,
                                          mappedNodes,
                                          shapeFunctionValues,
                                          shapeFunctionGradientValues);          
      #endif

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        real64 particleContributionToGrid;

        #ifdef GEOS_USE_DEVICE
          localIndex const mappedNode = mappedNodes[g];
          real64 const shapeFunctionValue = shapeFunctionValues[g];
          real64 shapeFunctionGradientValue [3];
          LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[g] );
        #else
          localIndex const mappedNode = mappedNodes[pp][g];
          real64 const shapeFunctionValue = shapeFunctionValues[pp][g];
          real64 shapeFunctionGradientValue [3];
          LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[pp][g] );
        #endif

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        RAJA::atomicOr( parallelDeviceAtomic{}, &gridCohesiveFieldFlag[mappedNode][fieldIndex], particleCohesiveZoneFlag[p]);

        if( computeXProfile == 1 )
        {
          for( localIndex i = 0 ; i < 3 ; ++i )
          { // map the volume-weighted normal stress to nodes for profile output.
            particleContributionToGrid = shapeFunctionValue * particleStress[p][i] * particleVolume[p];
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridNormalStress[mappedNode][fieldIndex][i], particleContributionToGrid );
          }

          particleContributionToGrid = shapeFunctionValue * particleDamage[p] * particleMass[p];
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMassWeightedDamage[mappedNode][fieldIndex], particleContributionToGrid );
        }

        particleContributionToGrid = particleMass[p] * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMass[mappedNode][fieldIndex], particleContributionToGrid );

        particleContributionToGrid = particleVolume[p] * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMaterialVolume[mappedNode][fieldIndex], particleContributionToGrid );

        // TODO: Normalizing by volume might be better
        particleContributionToGrid = particleMass[p] * ( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ) * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridDamage[mappedNode][fieldIndex], particleContributionToGrid );

        particleContributionToGrid = particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p];
        RAJA::atomicMax( parallelDeviceAtomic{}, &gridMaxDamage[mappedNode][fieldIndex], particleContributionToGrid );

        bool hasSurfaceNormal = LvArray::tensorOps::l2Norm< 3 >( particleSurfaceNormal[p] ) > 1e-16; // Change this to l2norm squared for computational efficiency
        real64 surfacePositionAlongNormal = 0.0;
        if( hasContact == 1 && hasSurfaceNormal )
        {
          // Surface positions
          particleContributionToGrid = shapeFunctionValue * particleMass[p];
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceFieldMass[mappedNode][fieldIndex], particleContributionToGrid );

          real64 surfacePositionRelativeToNode[3];
          LvArray::tensorOps::copy< 3 >( surfacePositionRelativeToNode, particleSurfacePosition[p] );
          LvArray::tensorOps::add< 3 >( surfacePositionRelativeToNode, particlePosition[p] );
          LvArray::tensorOps::subtract< 3 >( surfacePositionRelativeToNode, gridPosition[mappedNode] );

          // TODO: Check which of these produces better results
          // We need to use the surface position relative to the grid node
          // Two possible options map relative surface position either using particle surface normal or grid surface normal
          surfacePositionAlongNormal = LvArray::tensorOps::AiBi< 3 >( surfacePositionRelativeToNode, particleSurfaceNormal[p] );
        }

        for( localIndex i=0; i < numDims; ++i )
        {
          if( hasContact == 1 )
          {
            // Surface normals
            particleContributionToGrid = shapeFunctionGradientValue[i] * particleVolume[p];
            if( particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 )
            {
              // Also maps explicit particle surface normals which will dominate if m_explicitSurfaceNormalInfluence is large
              // If particle surface normal was disabled due to damage or CPDI domain scaling (e.g. zeroed) then the following does not add anything to gridSurfaceNormal
              particleContributionToGrid += explicitSurfaceNormalInfluence * particleSurfaceNormal[p][i] * shapeFunctionValue * particleVolume[p];
            }
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceNormal[mappedNode][fieldIndex][i], particleContributionToGrid );

            if( hasSurfaceNormal )
            {
              // Surface positions 
              particleContributionToGrid = shapeFunctionValue * particleMass[p] * surfacePositionAlongNormal * particleSurfaceNormal[p][i];
              RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfacePosition[mappedNode][fieldIndex][i], particleContributionToGrid );
            }
          }

          particleContributionToGrid = particleMass[p] * particleVelocity[p][i] * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMomentum[mappedNode][fieldIndex][i], particleContributionToGrid );

          particleContributionToGrid = particleBodyForce[p][i] * particleMass[p] + particleSurfaceTraction[p][i] + particleCohesiveForce[p][i];
          if( enableSurfaceTension == 1 )
          {
            particleContributionToGrid += particleSurfaceCurvature[p] * particleSurfaceNormal[p][i] * m_surfaceTensionCoefficient;
          }
          particleContributionToGrid *= shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridExternalForce[mappedNode][fieldIndex][i], particleContributionToGrid );

          particleContributionToGrid = particleVolume[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridCenterOfVolume[mappedNode][fieldIndex][i], particleContributionToGrid );

          particleContributionToGrid = particleMass[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridCenterOfMass[mappedNode][fieldIndex][i], particleContributionToGrid );

          // MH: This will modify the stress if there is a non-zero borehole pressure set by the boreholePressure
          // event.  This change is applied if the radius in the x-y plane, centered at origin, defining the extent of the borehole
          // pressure BC, should be bigger than the borehole but not near outer domain boundary.
          for( localIndex k = 0; k < numDims; ++k )
          {
            localIndex voigt = voigtMap[k][i];
            particleContributionToGrid = particleStress[p][voigt] - gridBoreholeStress[mappedNode][voigt];
            if( useArtificialViscosity == 1 )
            {
              particleContributionToGrid -= particleArtificialViscosity[p] * (k == i);
            }
            particleContributionToGrid *= -1 * shapeFunctionGradientValue[k] * particleVolume[p];
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridInternalForce[mappedNode][fieldIndex][i], particleContributionToGrid );
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop
}

void SolidMechanicsMPM::particleToGrid_noAtomics( real64 const time_n,
                                                  integer const cycleNumber,
                                                  ParticleManager & particleManager,
                                                  NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // On-the-fly shape function computations
  localIndex const numDims = m_numDims;
  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const useArtificialViscosity = m_useArtificialViscosity;
  int const computeXProfile = m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) );
  localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView2d< real64 > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );

  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView2d< int > const & gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );

  arrayView2d< real64 const > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

    // Get views to mapping arrays
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    ParticleType const particleType = subRegion.getParticleType();

    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // On-the-fly shape function computation
      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues);

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        localIndex const mappedNode = mappedNodes[g];
        real64 const shapeFunctionValue = shapeFunctionValues[g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );
        gridCohesiveFieldFlag[mappedNode][fieldIndex] |= particleCohesiveZoneFlag[p];

        if( computeXProfile == 1 )
        {
          for( localIndex i = 0 ; i < 3 ; ++i )
          { // map the volume-weighted normal stress to nodes for profile output.
            gridNormalStress[mappedNode][fieldIndex][i] += shapeFunctionValue * particleStress[p][i] * particleVolume[p];
          }

          gridMassWeightedDamage[mappedNode][fieldIndex] += shapeFunctionValue * particleDamage[p] * particleMass[p];
        }

        gridMass[mappedNode][fieldIndex] = particleMass[p] * shapeFunctionValue;

        gridMaterialVolume[mappedNode][fieldIndex] += particleVolume[p] * shapeFunctionValue;

        // TODO: Normalizing by volume might be better
        gridDamage[mappedNode][fieldIndex] += particleMass[p] * ( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ) * shapeFunctionValue;

        gridMaxDamage[mappedNode][fieldIndex] = LvArray::math::max( gridMaxDamage[mappedNode][fieldIndex], particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] );

        for( localIndex i = 0; i < numDims; ++i )
        {
          gridMomentum[mappedNode][fieldIndex][i] += particleMass[p] * particleVelocity[p][i] * shapeFunctionValues[g];

          gridExternalForce[mappedNode][fieldIndex][i] += ( particleBodyForce[p][i] * particleMass[p] + particleSurfaceTraction[p][i] + particleCohesiveForce[p][i] ) * shapeFunctionValue;

          gridCenterOfVolume[mappedNode][fieldIndex][i] += particleVolume[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;

          gridCenterOfMass[mappedNode][fieldIndex][i] += particleMass[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;

          // MH: This will modify the stress if there is a non-zero borehole pressure set by the boreholePressure
          // event.  This change is applied if the radius in the x-y plane, centered at origin, defining the extent of the borehole
          // pressure BC, should be bigger than the borehole but not near outer domain boundary.
          for( localIndex k=0; k < numDims; ++k )
          {
            localIndex voigt = voigtMap[k][i];
            gridInternalForce[mappedNode][fieldIndex][i] += -( ( particleStress[p][voigt] - gridBoreholeStress[mappedNode][voigt] ) - particleArtificialViscosity[p] * useArtificialViscosity * (k == i) ) * shapeFunctionGradientValues[g][k] * particleVolume[p];
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop
}

void SolidMechanicsMPM::particleToGrid_randomMix( real64 const time_n,
                                                  integer const cycleNumber,
                                                  ParticleManager & particleManager,
                                                  NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // On-the-fly shape function computations
  localIndex const numDims = m_numDims;
  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const useArtificialViscosity = m_useArtificialViscosity;
  int const computeXProfile = m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) );
  localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView2d< real64 > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );

  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView2d< int > const & gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );

  arrayView2d< real64 const > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

    // Get views to mapping arrays
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    ParticleType const particleType = subRegion.getParticleType();

    // Copy the activeParticleIndices list and shuffle it to reduce number of atomic hits
    // stdVector< localIndex > randomizer( activeParticleIndices.begin(), activeParticleIndices.end() );
    array1d< localIndex > randomizer( activeParticleIndices.size() );
    for( localIndex p = 0; p < activeParticleIndices.size(); ++p )
    {
      randomizer[p] = activeParticleIndices[p];
    }
    std::random_device rd;
    std::mt19937 generator(rd());
    std::shuffle( randomizer.begin(),  randomizer.end(), generator);
    arrayView1d< localIndex const > const randomized = randomizer;

    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = randomized[pp];

      // On-the-fly shape function computation
      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues);

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        real64 particleContributionToGrid;

        localIndex const mappedNode = mappedNodes[g];
        real64 const shapeFunctionValue = shapeFunctionValues[g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        RAJA::atomicOr( parallelDeviceAtomic{}, &gridCohesiveFieldFlag[mappedNode][fieldIndex], particleCohesiveZoneFlag[p]);

        if( computeXProfile == 1 )
        {
          for( int i = 0 ; i < 3 ; ++i )
          { // map the volume-weighted normal stress to nodes for profile output.
            particleContributionToGrid = shapeFunctionValue * particleStress[p][i] * particleVolume[p];
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridNormalStress[mappedNode][fieldIndex][i], particleContributionToGrid );
          }

          particleContributionToGrid = shapeFunctionValue * particleDamage[p] * particleMass[p];
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMassWeightedDamage[mappedNode][fieldIndex], particleContributionToGrid );
        }

        particleContributionToGrid = particleMass[p] * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMass[mappedNode][fieldIndex], particleContributionToGrid );

        particleContributionToGrid = particleVolume[p] * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMaterialVolume[mappedNode][fieldIndex], particleContributionToGrid );

        // TODO: Normalizing by volume might be better
        particleContributionToGrid = particleMass[p] * ( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ) * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridDamage[mappedNode][fieldIndex], particleContributionToGrid );

        particleContributionToGrid = particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p];
        RAJA::atomicMax( parallelDeviceAtomic{}, &gridMaxDamage[mappedNode][fieldIndex], particleContributionToGrid );

        for( int i=0; i < numDims; ++i )
        {
          particleContributionToGrid = particleMass[p] * particleVelocity[p][i] * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMomentum[mappedNode][fieldIndex][i], particleContributionToGrid );

          particleContributionToGrid = ( particleBodyForce[p][i] * particleMass[p] + particleSurfaceTraction[p][i] + particleCohesiveForce[p][i] ) * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridExternalForce[mappedNode][fieldIndex][i], particleContributionToGrid );

          particleContributionToGrid = particleVolume[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridCenterOfVolume[mappedNode][fieldIndex][i], particleContributionToGrid );

          particleContributionToGrid = particleMass[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridCenterOfMass[mappedNode][fieldIndex][i], particleContributionToGrid );

          // MH: This will modify the stress if there is a non-zero borehole pressure set by the boreholePressure
          // event.  This change is applied if the radius in the x-y plane, centered at origin, defining the extent of the borehole
          // pressure BC, should be bigger than the borehole but not near outer domain boundary.
          for( localIndex k=0; k < numDims; ++k )
          {
            localIndex voigt = voigtMap[k][i];
            particleContributionToGrid = -( ( particleStress[p][voigt] - gridBoreholeStress[mappedNode][voigt] ) - particleArtificialViscosity[p] * useArtificialViscosity * (k == i) ) * shapeFunctionGradientValues[g][k] * particleVolume[p];
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridInternalForce[mappedNode][fieldIndex][i], particleContributionToGrid );
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop
}

void SolidMechanicsMPM::particleToGrid_minimalAtomics( real64 const time_n,
                                                       integer const cycleNumber,
                                                       ParticleManager & particleManager,
                                                       NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // On-the-fly shape function computations
  localIndex const numDims = m_numDims;
  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const useArtificialViscosity = m_useArtificialViscosity;
  int const computeXProfile = m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) );
  localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView2d< real64 > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );

  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView2d< int > const & gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );

  arrayView2d< real64 const > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

    // Get views to mapping arrays
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    ParticleType const particleType = subRegion.getParticleType();

    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // printf("Started GPU loop\n");

      // Define temporary contributions and initialize them
      int gridCohesiveFieldFlagContribution = 0;
      real64 gridMassWeightedDamageContribution = 0.0;
      real64 gridMassContribution = 0.0;
      real64 gridMaterialVolumeContribution = 0.0;
      real64 gridDamageContribution = 0.0;
      real64 gridMaxDamageContribution = 0.0;

      real64 gridNormalStressContribution[3];
      real64 gridMomentumContribution[3];
      real64 gridExternalForceContribution[3];
      real64 gridCenterOfVolumeContribution[3];
      real64 gridCenterOfMassContribution[3];
      real64 gridInternalForceContribution[6];

      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];

      // On-the-fly shape function computation
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues);
      // printf("Mapped nodes and computed shapefunctions\n");

      // printf("Particle %d\n", p);
      // printf("Mapped Nodes:\n");
      // for( int i = 0; i < 64; ++i)
      // {
      //   printf("%d , ", mappedNodes[i]);
      // }
      // printf("\n");

      // Copy mappedNodes
      // Sort copied mappedNodes
      // Iterate over sorted mappedNodes
      // Only sum the nodal contributions after the mappedNodex index changes
      // array1d< localIndex > sortedMappedNodes( 64 );
      // TODO will need to clear the array each iteration since we capture it as a view from outside. Received error that we cannot create chai buffers on device
      // If we mix particles with different number of mapping vertices and sort then we could accidentally sort grid indices that don't belong to the particle
      // printf("Creating sorted arrays\n");
      localIndex sortedMappedNodes[64];
      // LvArray::tensorOps::copy< 64 >( sortedMappedNodes, mappedNodes );
      localIndex sortedIndices[64] = {};
      for( localIndex i = 0; i < 64; ++i ) // Not general for all particle types
      {
        sortedMappedNodes[i] = mappedNodes[i];
        sortedIndices[i] = i;
      }
      // printf("Starting quick sort\n");

      // quickSort( sortedMappedNodes, sortedIndices, 0, 63 );
      quickSortIterative( sortedMappedNodes, sortedIndices, 0, 63 );

      // printf("Quicksort complete\n");

      // Use sorted indices to sort the shape functions and shape function gradients
      real64 sortedShapeFunctionValues[64];
      real64 sortedShapeFunctionGradientValues[64][3] = { {0.0} };
      // printf("Sorted:\n");
      for( localIndex i = 0; i < 64; ++i)
      {
        // printf("%d (%d), ", sortedMappedNodes[i], sortedIndices[i]);
        sortedShapeFunctionValues[i] = shapeFunctionValues[sortedIndices[i]];
        LvArray::tensorOps::copy< 3 >( sortedShapeFunctionGradientValues[i], shapeFunctionGradientValues[sortedIndices[i]] );
      }
      // printf("\n");

      // printf("Sorted shape functions\n");

      localIndex mappedNode = sortedMappedNodes[0];
      localIndex fieldIndex = partitionField( numContactGroups,
                                       damageFieldPartitioning,
                                       particleGroup[p],
                                       particleDamageGradient[p],
                                       particleSurfaceNormal[p],
                                       gridDamageGradient[mappedNode] );

      // printf("Computed field index\n");

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {

        real64 const shapeFunctionValue = sortedShapeFunctionValues[g];
        const real64 * const shapeFunctionGradientValue = sortedShapeFunctionGradientValues[g];

        gridCohesiveFieldFlagContribution |= particleCohesiveZoneFlag[p];

        if( computeXProfile == 1 )
        {
          for( int i = 0 ; i < 3 ; ++i )
          { // map the volume-weighted normal stress to nodes for profile output.
            gridNormalStressContribution[i] += shapeFunctionValue * particleStress[p][i] * particleVolume[p];
          }

          gridMassWeightedDamageContribution += shapeFunctionValue * particleDamage[p] * particleMass[p];
        }

        gridMassContribution += particleMass[p] * shapeFunctionValue;
        gridMaterialVolumeContribution += particleVolume[p] * shapeFunctionValue;
        gridDamageContribution += particleMass[p] * ( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ) * shapeFunctionValue;
        gridMaxDamageContribution = fmax(gridMaxDamageContribution, particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p]);

        for( localIndex i = 0; i < numDims; ++i )
        {
          gridMomentumContribution[i] += particleMass[p] * particleVelocity[p][i] * shapeFunctionValue;
          gridExternalForceContribution[i] += ( particleBodyForce[p][i] * particleMass[p] + particleSurfaceTraction[p][i] + particleCohesiveForce[p][i] ) * shapeFunctionValue;
          gridCenterOfVolumeContribution[i] += particleVolume[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
          gridCenterOfMassContribution[i] += particleMass[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;

          // MH: This will modify the stress if there is a non-zero borehole pressure set by the boreholePressure
          // event.  This change is applied if the radius in the x-y plane, centered at origin, defining the extent of the borehole
          // pressure BC, should be bigger than the borehole but not near outer domain boundary.
          for( localIndex k = 0; k < numDims; ++k )
          {
            localIndex voigt = voigtMap[k][i];
            gridInternalForceContribution[voigt] += -( ( particleStress[p][voigt] - gridBoreholeStress[mappedNode][voigt] ) - particleArtificialViscosity[p] * useArtificialViscosity * (k == i) ) * shapeFunctionGradientValue[k] * particleVolume[p];
          }
        }

        if( sortedMappedNodes[g] != sortedMappedNodes[g+1] || g == 8 * numberOfVerticesPerParticle - 1 )
        {
          RAJA::atomicOr( parallelDeviceAtomic{}, &gridCohesiveFieldFlag[mappedNode][fieldIndex], gridCohesiveFieldFlagContribution );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMassWeightedDamage[mappedNode][fieldIndex], gridMassWeightedDamageContribution );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMass[mappedNode][fieldIndex], gridMassContribution );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMaterialVolume[mappedNode][fieldIndex], gridMaterialVolumeContribution );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridDamage[mappedNode][fieldIndex], gridDamageContribution );
          RAJA::atomicMax( parallelDeviceAtomic{}, &gridMaxDamage[mappedNode][fieldIndex], gridMaxDamageContribution );

          for( localIndex i = 0; i < numDims; ++i )
          {
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridNormalStress[mappedNode][fieldIndex][i], gridNormalStressContribution[i] );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMomentum[mappedNode][fieldIndex][i], gridMomentumContribution[i] );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridExternalForce[mappedNode][fieldIndex][i], gridExternalForceContribution[i] );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridCenterOfVolume[mappedNode][fieldIndex][i], gridCenterOfVolumeContribution[i] );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridCenterOfMass[mappedNode][fieldIndex][i], gridCenterOfMassContribution[i] );
          }

          for( localIndex i = 0; i < 6; ++i)
          {
            RAJA::atomicAdd( parallelDeviceAtomic{}, &gridInternalForce[mappedNode][fieldIndex][i], gridInternalForceContribution[i] );
          }

          if( g != 8 * numberOfVerticesPerParticle - 1 )
          {
            // Recompute fieldIndex
            mappedNode = sortedMappedNodes[g+1];
            fieldIndex = partitionField( numContactGroups,
                                         damageFieldPartitioning,
                                         particleGroup[p],
                                         particleDamageGradient[p],
                                         particleSurfaceNormal[p],
                                         gridDamageGradient[mappedNode] );

            // Zero grid node contributions
            gridCohesiveFieldFlagContribution = 0;
            gridMassWeightedDamageContribution = 0.0;
            gridMassContribution = 0.0;
            gridMaterialVolumeContribution = 0.0;
            gridDamageContribution = 0.0;
            gridMaxDamageContribution = 0.0;

            LvArray::tensorOps::fill< 3 >( gridNormalStressContribution, 0.0 );
            LvArray::tensorOps::fill< 3 >( gridMomentumContribution, 0.0 );
            LvArray::tensorOps::fill< 3 >( gridExternalForceContribution, 0.0 );
            LvArray::tensorOps::fill< 3 >( gridCenterOfVolumeContribution, 0.0 );
            LvArray::tensorOps::fill< 3 >( gridCenterOfMassContribution, 0.0 );
            LvArray::tensorOps::fill< 6 >( gridInternalForceContribution, 0.0 );
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop
}

void SolidMechanicsMPM::particleToGrid_reduction( real64 const time_n,
                                                  integer const cycleNumber,
                                                  ParticleManager & particleManager,
                                                  NodeManager & nodeManager )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( cycleNumber );
  GEOS_UNUSED_VAR( particleManager );
  GEOS_UNUSED_VAR( nodeManager );

  // GEOS_MARK_FUNCTION;

  // // On-the-fly shape function computations
  // localIndex const numDims = m_numDims;
  // localIndex const numContactGroups = m_numContactGroups;
  // localIndex const numVelocityFields = m_numVelocityFields;
  // int const damageFieldPartitioning = m_damageFieldPartitioning;
  // int const useArtificialViscosity = m_useArtificialViscosity;
  // int const computeXProfile = m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) );
  // localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  // real64 hEl[3];
  // LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  // real64 xLocalMin[3];
  // LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  // arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // // Grid fields
  // arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  // arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  // arrayView2d< real64 const > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  // // Grid fields to map particles to
  // arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  // arrayView2d< real64 > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  // arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  // arrayView2d< real64 > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  // arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );

  // arrayView2d< int > const & gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );

  // arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  // arrayView3d< real64 > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  // arrayView3d< real64 > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );

  // arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  // arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );
  // arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );

  // #if defined(GEOS_USE_CUDA)
  //   using reduce_policy = RAJA::cuda_multi_reduce_atomic;
  // #endif

  // #if defined(GEOS_USE_HIP)
  //   using reduce_policy = RAJA::hip_multi_reduce_atomic;
  // #endif

  // localIndex const numNodes = nodeManager.size();

  // // Scalar fields
  // localIndex count = numNodes * numVelocityFields;
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridMassReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridMaterialVolumeReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridDamageReduction( count, 0.0 );
  // RAJA::MultiReduceMax< reduce_policy, real64 > gridMaxDamageReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridMassWeightedDamageReduction( count, 0.0 );
  // RAJA::MultiReduceBitOr< reduce_policy, int > gridCohesiveFieldFlagReduction( count, 0.0 );

  // // Vector fields
  // count = numNodes * numVelocityFields * 3;
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridMomentumReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridExternalForceReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridCenterOfMassReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridCenterOfVolumeReduction( count, 0.0 );
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridNormalStressReduction( count, 0.0 );

  // count = numNodes * numVelocityFields * 6;
  // RAJA::MultiReduceSum< reduce_policy, real64 > gridInternalForceReduction( count, 0.0 );

  // localIndex subRegionIndex = 0;
  // particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  // {
  //   // Particle fields
  //   arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
  //   arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
  //   arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
  //   arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
  //   arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
  //   arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
  //   arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
  //   arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
  //   arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

  //   arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
  //   arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
  //   arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
  //   arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
  //   arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
  //   arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
  //   arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

  //   // Get views to mapping arrays
  //   int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

  //   // Map to grid
  //   SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
  //   ParticleType const particleType = subRegion.getParticleType();

  //   forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
  //   {
  //     localIndex const p = activeParticleIndices[pp];

  //     // On-the-fly shape function computation
  //     localIndex mappedNodes[64];
  //     real64 shapeFunctionValues[64];
  //     real64 shapeFunctionGradientValues[64][3];
  //     mapNodesAndComputeShapeFunctions( ijkMap,
  //                                       xLocalMin,
  //                                       hEl,
  //                                       particleType,
  //                                       particlePosition[p],
  //                                       particleRVectors[p],
  //                                       gridPosition,
  //                                       mappedNodes,
  //                                       shapeFunctionValues,
  //                                       shapeFunctionGradientValues);

  //     for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
  //     {
  //       localIndex const mappedNode = mappedNodes[g];
  //       real64 const shapeFunctionValue = shapeFunctionValues[g];
  //       localIndex const fieldIndex = partitionField( numContactGroups,
  //                                              damageFieldPartitioning,
  //                                              particleGroup[p],
  //                                              particleDamageGradient[p],
  //                                              particleSurfaceNormal[p],
  //                                              gridDamageGradient[mappedNode] );

  //       localIndex scalarIndex = numNodes * fieldIndex + mappedNode;
  //       localIndex vectorIndex;

  //       gridCohesiveFieldFlagReduction[scalarIndex] |= particleCohesiveZoneFlag[p];

  //       if( computeXProfile == 1 )
  //       {
  //         for( localIndex i = 0 ; i < 3 ; ++i )
  //         { // map the volume-weighted normal stress to nodes for profile output.
  //           vectorIndex = ( numNodes * 3 ) * fieldIndex + numNodes * i + mappedNode;
  //           gridNormalStressReduction[vectorIndex] += shapeFunctionValue * particleStress[p][i] * particleVolume[p];
  //         }

  //         gridMassWeightedDamageReduction[scalarIndex] += shapeFunctionValue * particleDamage[p] * particleMass[p];
  //       }

  //       gridMassReduction[scalarIndex] += particleMass[p] * shapeFunctionValue;
  //       gridMaterialVolumeReduction[scalarIndex] += particleVolume[p] * shapeFunctionValue;
  //       gridDamageReduction[scalarIndex] += particleMass[p] * ( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ) * shapeFunctionValue;
  //       gridMaxDamageReduction[scalarIndex].max( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ); // needs to be a max check

  //       for( localIndex i=0; i < numDims; ++i )
  //       {
  //         vectorIndex = ( numNodes * 3 ) * fieldIndex + numNodes * i + mappedNode;
  //         gridMomentumReduction[vectorIndex] += particleMass[p] * particleVelocity[p][i] * shapeFunctionValue;
  //         gridExternalForceReduction[vectorIndex] += ( particleBodyForce[p][i] * particleMass[p] + particleSurfaceTraction[p][i] + particleCohesiveForce[p][i] ) * shapeFunctionValue;
  //         gridCenterOfVolumeReduction[vectorIndex] += particleVolume[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
  //         gridCenterOfMassReduction[vectorIndex] += particleMass[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;

  //         // MH: This will modify the stress if there is a non-zero borehole pressure set by the boreholePressure
  //         // event.  This change is applied if the radius in the x-y plane, centered at origin, defining the extent of the borehole
  //         // pressure BC, should be bigger than the borehole but not near outer domain boundary.
  //         for( int k=0; k < numDims; ++k )
  //         {
  //           int voigt = voigtMap[k][i];
  //           vectorIndex = ( numNodes * 6 ) * fieldIndex + numNodes * voigt + mappedNode;
  //           gridInternalForceReduction[vectorIndex] += -( ( particleStress[p][voigt] - gridBoreholeStress[mappedNode][voigt] ) - particleArtificialViscosity[p] * useArtificialViscosity * (k == i) ) * shapeFunctionGradientValues[g][k] * particleVolume[p];
  //         }
  //       }
  //     }
  //   } ); // particle loop

  //   // Copy reduction values back to grid fields (in serial for now, but might need to try parallel on device)
  //   for( localIndex fieldIndex=0; fieldIndex < numVelocityFields; ++fieldIndex )
  //   {
  //     for( int nodeIndex=0; nodeIndex < numNodes; nodeIndex++ )
  //     {
  //       localIndex scalarIndex = numNodes * fieldIndex + nodeIndex;
  //       gridMass[nodeIndex][fieldIndex] = gridMassReduction[scalarIndex].get();
  //       gridMaterialVolume[nodeIndex][fieldIndex] = gridMaterialVolumeReduction[scalarIndex].get();
  //       gridDamage[nodeIndex][fieldIndex] =  gridDamageReduction[scalarIndex].get();
  //       gridMaxDamage[nodeIndex][fieldIndex] =  gridMaxDamageReduction[scalarIndex].get();
  //       gridMassWeightedDamage[nodeIndex][fieldIndex] = gridMassWeightedDamageReduction[scalarIndex].get();
  //       gridCohesiveFieldFlag[nodeIndex][fieldIndex] = gridCohesiveFieldFlagReduction[scalarIndex].get();

  //       localIndex vectorIndex;
  //       for(localIndex i = 0; i < numDims; ++i)
  //       {
  //         vectorIndex = ( numNodes * 3 ) * fieldIndex + numNodes * i + nodeIndex;
  //         gridMomentum[nodeIndex][fieldIndex][i] = gridMomentumReduction[vectorIndex].get();
  //         gridExternalForce[nodeIndex][fieldIndex][i] = gridExternalForceReduction[vectorIndex].get();
  //         gridCenterOfMass[nodeIndex][fieldIndex][i] = gridCenterOfMassReduction[vectorIndex].get();
  //         gridCenterOfVolume[nodeIndex][fieldIndex][i] = gridCenterOfVolumeReduction[vectorIndex].get();
  //         gridNormalStress[nodeIndex][fieldIndex][i] = gridNormalStressReduction[vectorIndex].get();
  //       }

  //       // Voigt notation
  //       for(localIndex i = 0; i < 6; ++i)
  //       {
  //         vectorIndex = ( numNodes * 6 ) * fieldIndex + numNodes * i + nodeIndex;
  //         gridInternalForce[nodeIndex][fieldIndex][i] = gridInternalForceReduction[vectorIndex].get();
  //       }
  //     }
  //   }

  //   // Increment subregion index
  //   ++subRegionIndex;
  // } ); // subregion loop
}

void SolidMechanicsMPM::particleColorSort( ParticleManager & particleManager )
{
  GEOS_UNUSED_VAR( particleManager );

  // real64 nEl[3];
  // LvArray::tensorOps::copy< 3 >( nEl, m_nEl );
  // real64 hEl[3];
  // LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  // real64 xLocalMin[3];
  // LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  // // arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // int colorDims[3];
  // colorDims[0] = 4;
  // colorDims[1] = 4;
  // colorDims[2] = 4;
  // int totalColors = colorDims[0]*colorDims[1]*colorDims[2];

  // array2d< array1d < int> > > idx( numSubRegions, totalColors );

  // localIndex subRegionIndex = 0;
  // particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  // {
  //   arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
  //   arrayView1d< int > const particleColor = subRegion.getField< fields::mpm::particleColor >();

  //   SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
  //   forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
  //   {
  //     localIndex const p = activeParticleIndices[pp];

  //     localIndex i = LvArray::math::floor( ( particlePosition[p][0] - xLocalMin[0] ) / hEl[0] );
  //     localIndex j = LvArray::math::floor( ( particlePosition[p][1] - xLocalMin[1] ) / hEl[1] );
  //     localIndex k = LvArray::math::floor( ( particlePosition[p][2] - xLocalMin[2] ) / hEl[2] );

  //     particleColor[p] = colorDims[0] * colorDims[1] * (k % colorDims[2]) + colorDims[0] * (j % colorDims[1]) + (i % colorDims[0]);
  //   });

  //   ++subRegionIndex;
  // });
}

void SolidMechanicsMPM::particleToGrid_colors( real64 const time_n,
                                                  integer const cycleNumber,
                                                  ParticleManager & particleManager,
                                                  NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( cycleNumber );
  GEOS_UNUSED_VAR( particleManager );
  GEOS_UNUSED_VAR( nodeManager );

  // On-the-fly shape function computations
  localIndex const numDims = m_numDims;
  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const useArtificialViscosity = m_useArtificialViscosity;
  int const computeXProfile = m_computeXProfile == 1 && ( ( m_nextXProfileWriteTime <= time_n ) || ( cycleNumber == 0 ) );
  localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );
  arrayView2d< real64 > const & gridMaxDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() );
  arrayView2d< real64 > const & gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );

  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView2d< int > const & gridCohesiveFieldFlag = nodeManager.getReference< array2d< int > >( viewKeyStruct::gridCohesiveFieldFlagString() );

  arrayView2d< real64 const > const gridBoreholeStress = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridBoreholeStressString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfaceTraction = subRegion.getParticleSurfaceTraction();

    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
    arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
    arrayView2d< real64 const > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 > const particleCohesiveForce = subRegion.getField< fields::mpm::particleCohesiveForce >();
    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Map to grid
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    ParticleType const particleType = subRegion.getParticleType();

    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // On-the-fly shape function computation
      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues);

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        localIndex const mappedNode = mappedNodes[g];
        real64 const shapeFunctionValue = shapeFunctionValues[g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );
        gridCohesiveFieldFlag[mappedNode][fieldIndex] |= particleCohesiveZoneFlag[p];

        if( computeXProfile == 1 )
        {
          for( localIndex i = 0 ; i < 3 ; ++i )
          { // map the volume-weighted normal stress to nodes for profile output.
            gridNormalStress[mappedNode][fieldIndex][i] += shapeFunctionValue * particleStress[p][i] * particleVolume[p];
          }

          gridMassWeightedDamage[mappedNode][fieldIndex] += shapeFunctionValue * particleDamage[p] * particleMass[p];
        }

        gridMass[mappedNode][fieldIndex] = particleMass[p] * shapeFunctionValue;

        gridMaterialVolume[mappedNode][fieldIndex] += particleVolume[p] * shapeFunctionValue;

        // TODO: Normalizing by volume might be better
        gridDamage[mappedNode][fieldIndex] += particleMass[p] * ( particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] ) * shapeFunctionValue;

        gridMaxDamage[mappedNode][fieldIndex] = LvArray::math::max( gridMaxDamage[mappedNode][fieldIndex], particleSurfaceFlag[p] == 1 || particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ? 1.0 : particleDamage[p] );

        for( localIndex i=0; i < numDims; ++i )
        {
          gridMomentum[mappedNode][fieldIndex][i] += particleMass[p] * particleVelocity[p][i] * shapeFunctionValues[g];

          gridExternalForce[mappedNode][fieldIndex][i] += ( particleBodyForce[p][i] * particleMass[p] + particleSurfaceTraction[p][i] + particleCohesiveForce[p][i] ) * shapeFunctionValue;

          gridCenterOfVolume[mappedNode][fieldIndex][i] += particleVolume[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;

          gridCenterOfMass[mappedNode][fieldIndex][i] += particleMass[p] * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;

          // MH: This will modify the stress if there is a non-zero borehole pressure set by the boreholePressure
          // event.  This change is applied if the radius in the x-y plane, centered at origin, defining the extent of the borehole
          // pressure BC, should be bigger than the borehole but not near outer domain boundary.
          for( localIndex k = 0; k < numDims; ++k )
          {
            localIndex voigt = voigtMap[k][i];
            gridInternalForce[mappedNode][fieldIndex][i] += -( ( particleStress[p][voigt] - gridBoreholeStress[mappedNode][voigt] ) - particleArtificialViscosity[p] * useArtificialViscosity * (k == i) ) * shapeFunctionGradientValues[g][k] * particleVolume[p];
          }
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop
}

// Computes the surface tension by first computing the surface curvate for each particle
void SolidMechanicsMPM::computeSPHSurfaceCurvature( ParticleManager & particleManager )
{
 GEOS_MARK_FUNCTION;

  int const planeStrain = m_planeStrain;
  real64 const neighborRadius = m_neighborRadius;

  // Get accessors for volume, position, damage, surface flag, and cohesive flag
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > const particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > const particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > const particleSurfaceNormalAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleSurfaceNormal" );
 
  // Get views of accessors
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > const particleVolumeView = particleVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particlePositionView = particlePositionAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > const particleSurfaceNormalView = particleSurfaceNormalAccessor.toNestedViewConst();

  // Perform neighbor operations
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();

    // Get const views for array of arrays
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and damage field gradient
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    // arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView1d< real64 > const particleSurfaceCurvature = subRegion.getField< fields::mpm::particleSurfaceCurvature >();

    // arrayView1d< real64 > const & particleVolume = subRegion.getParticleVolume();

    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Call kernel field gradient function
      computeKernelFieldDivergenceDevice( particlePosition[p],
                                          neighborRadius,
                                          planeStrain,
                                          numNeighborsAll[pp],
                                          neighborRegions[pp],
                                          neighborSubRegions[pp],
                                          neighborIndices[pp],
                                          particleVolumeView,
                                          particlePositionView,
                                          particleSurfaceNormalView,
                                          particleSurfaceCurvature[p] );
    } );
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::gridTrialUpdate( real64 dt,
                                         NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  localIndex const numDims = m_numDims;
  real64 const smallMass = m_smallMass;
  localIndex numVelocityFields = m_numVelocityFields;

  // Grid fields
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView3d< real64 > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 > const & gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVelocityString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );
  arrayView3d< real64 const > const & gridInternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridInternalForceString() );
  arrayView3d< real64 const > const & gridExternalForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridExternalForceString() );
  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView2d< real64 > const & gridDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() );

  arrayView2d< real64 const > const gridMaterialVolume = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() );
  arrayView3d< real64 > const & gridCenterOfVolume = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfVolumeString() );

  localIndex const numNodes = nodeManager.size();
  
  // RAJA::MultiReduceSum< RAJA::seq_multi_reduce, real64 > internalForce(3, 0.0);
  // forAll< parallelDevicePolicy<> >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
  forAll< serialPolicy >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    // Loop over velocity fields
    for( localIndex fieldIndex=0; fieldIndex< numVelocityFields; ++fieldIndex )
    {
      if( gridMass[g][fieldIndex] > smallMass ) // small mass threshold
      {
        gridDamage[g][fieldIndex] /= gridMass[g][fieldIndex];
        for( localIndex i=0; i < numDims; ++i )
        {
          real64 totalForce = gridInternalForce[g][fieldIndex][i] + gridExternalForce[g][fieldIndex][i];
          gridAcceleration[g][fieldIndex][i] = totalForce / gridMass[g][fieldIndex];
          gridDVelocity[g][fieldIndex][i] = gridAcceleration[g][fieldIndex][i] * dt;
          gridMomentum[g][fieldIndex][i] += totalForce * dt;
          gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i] / gridMass[g][fieldIndex];
          gridCenterOfMass[g][fieldIndex][i] /= gridMass[g][fieldIndex];
          gridCenterOfVolume[g][fieldIndex][i] /= gridMaterialVolume[g][fieldIndex];
        }
      }
      else
      {
        gridDamage[g][fieldIndex] = 0.0;
        for( localIndex i = 0; i < numDims; ++i )
        {
          gridAcceleration[g][fieldIndex][i] = 0.0;
          gridDVelocity[g][fieldIndex][i] = 0.0;
          gridVelocity[g][fieldIndex][i] = 0.0;
          gridMomentum[g][fieldIndex][i] = 0.0;
          gridCenterOfVolume[g][fieldIndex][i] = 0.0;
          gridCenterOfMass[g][fieldIndex][i] = 0.0;
        }
      }
    }
  } );
}

void SolidMechanicsMPM::enforceContact( real64 dt,
                                        // DomainPartition & domain,
                                        // ParticleManager & particleManager,
                                        NodeManager & nodeManager
                                        // MeshLevel & mesh
                                      )
{
  GEOS_MARK_FUNCTION;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView3d< real64 > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 > const & gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVelocityString() );
  arrayView3d< real64 > const & gridMomentum = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridMomentumString() );
  arrayView3d< real64 > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );
  arrayView3d< real64 > const & gridContactForce = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridContactForceString() );
  arrayView3d< real64 > const & gridCenterOfMass = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() );
  arrayView3d< real64 > const & gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );

  // Now computed in particleToGrid mapping
  // // Compute grid surface normals
  // computeGridSurfaceNormals( particleManager,
  //                            nodeManager );
  //
  // // Sync surface normals
  // syncGridFields( { viewKeyStruct::gridSurfaceNormalString() }, domain, nodeManager, mesh, MPI_SUM );

  // Apply symmetry boundary conditions to surface normals
  enforceGridVectorFieldSymmetryBC( gridSurfaceNormal, gridPosition, nodeManager.sets() );

  // // Normalize grid surface normals
  // normalizeGridSurfaceNormals( nodeManager );

  // // CC: Testing idea
  // // Compute grid surface normal weights
  // computeGridSurfaceNormalWeights( particleManager,
  //                                  nodeManager );

  // Apply symmetry boundary conditions to material position
  // TODO also enforce symmetry boundary conditions on grid center of volume
  enforceGridVectorFieldSymmetryBC( gridCenterOfMass, gridPosition, nodeManager.sets() );

  // Now computed in particleToGrid mapping
  // // Project particle surface positions to grid
  // computeGridSurfacePositions( particleManager, nodeManager );
  //
  // // Sync grid surface positions
  // syncGridFields( { viewKeyStruct::gridSurfacePositionString(), viewKeyStruct::gridSurfaceFieldMassString() }, domain, nodeManager, mesh, MPI_SUM );

  // Normalize grid surface normals and positions
  normalizeGridSurfaceNormalsAndPositions( nodeManager );

  // Compute contact forces
  computeContactForces( dt,
                        nodeManager );

  // Update grid momenta and velocities based on contact forces
  real64 const smallMass = m_smallMass;
  localIndex const numDims = m_numDims;
  localIndex const numVelocityFields = m_numVelocityFields;
  forAll< serialPolicy >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    for( localIndex fieldIndex=0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      if( gridMass[g][fieldIndex] > smallMass ) // small mass threshold
      {
        for( localIndex i=0; i < numDims; ++i )
        {
          real64 contactForce = gridContactForce[g][fieldIndex][i];
          gridAcceleration[g][fieldIndex][i] += contactForce / gridMass[g][fieldIndex];
          gridMomentum[g][fieldIndex][i] += contactForce * dt;
          gridDVelocity[g][fieldIndex][i] = contactForce / gridMass[g][fieldIndex] * dt; // Change in velocity from enforcing contact needed for XPIC update, CC: TODO clean up this code to avoid redundant operations
          gridVelocity[g][fieldIndex][i] = gridMomentum[g][fieldIndex][i] / gridMass[g][fieldIndex];
        }
      }
    }
  } );
}

void SolidMechanicsMPM::logisticRegression( ParticleManager & particleManager,
                                            NodeManager & nodeManager )
{
  GEOS_MARK_FUNCTION;

  real64 const smallMass = m_smallMass;
  localIndex const numContactGroups = m_numContactGroups;
  localIndex const numVelocityFields = m_numVelocityFields;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );

  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView3d< real64 const > const gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );
  arrayView3d< real64 > const gridLRSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfaceNormalString() );
  arrayView3d< real64 > const gridLRSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfacePositionString() );

  ParticleManager::ParticleViewAccessor< arrayView1d< localIndex const > > particleGroupAccessor = particleManager.constructArrayViewAccessor< localIndex, 1 >( "particleGroup" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particleSurfaceNormalAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleSurfaceNormal" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particleDamageGradientAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleDamageGradient" );

  ParticleManager::ParticleViewConst< arrayView1d< localIndex const > > particleGroupView =  particleGroupAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePositionView =  particlePositionAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleSurfaceNormalView =  particleSurfaceNormalAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleDamageGradientView =  particleDamageGradientAccessor.toNestedViewConst();

  arrayView1d< localIndex const > const numNeighborsAll = m_nodalNeighborList.m_numParticles.toViewConst();
  ArrayOfArraysView< localIndex const > const neighborRegions = m_nodalNeighborList.m_toParticleRegion.toViewConst();
  ArrayOfArraysView< localIndex const > const neighborSubRegions = m_nodalNeighborList.m_toParticleSubRegion.toViewConst();
  ArrayOfArraysView< localIndex const > const neighborIndices = m_nodalNeighborList.m_toParticleIndex.toViewConst();

  // Diagonal penalty matrix from paper
  real64 const lambda = 1e-7 * LvArray::math::square( ( m_hEl[0] + m_hEl[1] + m_hEl[2] ) / 3 ); // Paper doesn't say what equivalent grid size to use (min, max, or average)
  real64 const w_p = 1; // Particle weight, was 1 in paper

  // Iterate over nodes check if they have nonnegligible mass in more than one velocity field
  // If so we need to determine which particle is nearest the node and use their neighbor list to perform logistic regression
  localIndex const numNodes = nodeManager.size(); 
  forAll< serialPolicy >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
     for( localIndex a = 0; a < numVelocityFields-1; ++a )
     {
        if( gridMass[g][a] < smallMass )
        {
          continue;
        }
        
        for( localIndex b = a+1; b < numVelocityFields; ++b )
        {
          if( gridMass[g][b] < smallMass )
          {
            continue;
          }

          localIndex numNeighboringParticles = numNeighborsAll[g];

          real64 const mA = gridMass[g][a];
          real64 const mB = gridMass[g][b];

          // Use mass weighted normals as starting guess
          real64 n0[3];
          LvArray::tensorOps::scaledCopy< 3 >( n0, gridSurfaceNormal[g][b], mB );
          LvArray::tensorOps::scaledAdd< 3 >( n0, gridSurfaceNormal[g][a], -mA );
          LvArray::tensorOps::normalize< 3 >( n0 );
          
          real64 phi[4] = { n0[0], n0[1], n0[2], 0.0 }; // Initial guess for normal and offset phi={n_1, n_2, n_2, offset}

          real64 n_old[3]; // Eventually make the initial guess the mass weighted normals between the two fields
          real64 s_old[3];
          LvArray::tensorOps::copy< 3 >( n_old, n0 );
          LvArray::tensorOps::copy< 3 >( s_old, n0 );

          // Guess normal and offset of contact plane first (simple to use a weighting on the grid mapped normals from both fields and zero offset)
          bool converged = false;
          bool errored = false;
          localIndex const maxIters = 100; // Maximum number of iterations for logistic regression
          for(localIndex iter=0; iter < maxIters; ++iter)
          {
            real64 mat[4][4]; // J^T*W*J+Lambda
            for( localIndex i=0; i < 4; ++i )
            {
              for( localIndex j=0; j < 4; ++j )
              {
                mat[i][j] = lambda * ( i == j && i < 3 );
                for( localIndex n=0; n < numNeighboringParticles; ++n )
                {
                  localIndex const regionIndex = neighborRegions[g][n];
                  localIndex const subRegionIndex = neighborSubRegions[g][n];
                  localIndex const particleIndex = neighborIndices[g][n];

                  localIndex const fieldIndex = partitionField( numContactGroups,
                                                                damageFieldPartitioning,
                                                                particleGroupView[regionIndex][subRegionIndex][particleIndex],
                                                                particleDamageGradientView[regionIndex][subRegionIndex][particleIndex],
                                                                particleSurfaceNormalView[regionIndex][subRegionIndex][particleIndex],
                                                                gridDamageGradient[g] );
                  if( fieldIndex != a && fieldIndex != b )
                  {
                    continue;
                  }

                  real64 x_p[4] = { 1.0 };
                  x_p[0] = particlePositionView[regionIndex][subRegionIndex][particleIndex][0] - gridPosition[g][0]; // Unsure if position of material point is relative to grid node
                  x_p[1] = particlePositionView[regionIndex][subRegionIndex][particleIndex][1] - gridPosition[g][1];
                  x_p[2] = particlePositionView[regionIndex][subRegionIndex][particleIndex][2] - gridPosition[g][2];
                  x_p[3] = 1.0;
                  real64 exp_x_phi = LvArray::math::exp( -LvArray::tensorOps::AiBi< 4 >( x_p, phi ) );
                  real64 psi_sqr = LvArray::math::square( 2 * exp_x_phi / LvArray::math::square( 1.0 + exp_x_phi ) );

                  mat[i][j] += psi_sqr * w_p * x_p[i] * x_p[j];
                }
              }
            }

            real64 matInv[4][4];
            real64 matDet = LvArray::tensorOps::determinant< 4 >( mat );

            // printf("mat: {{ %.25f, %.25f, %.25f, %.25f}, { %.25f, %.25f, %.25f, %.25f}, { %.25f, %.25f, %.25f, %.25f}, { %.25f, %.25f, %.25f, %.25f}}, det: %.25f\n", mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0], mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1], mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2], mat[3][3], matDet);

            if( LvArray::math::abs(matDet) < 1e-25 )
            {
              GEOS_LOG_RANK("Matrix became singular!");
              errored = true;
              break;
            } 

            LvArray::tensorOps::invert< 4 >( matInv, mat );  

            real64 vec[4]; // J^T*W*(c-f(\phi^k))-Lambda*phi^k
            for( localIndex i = 0; i < 4; ++i )
            {
              vec[i] = -lambda * ( i < 3 ) * phi[i];
              for( localIndex n=0; n < numNeighboringParticles; ++n )
              {
                localIndex const regionIndex = neighborRegions[g][n];
                localIndex const subRegionIndex = neighborSubRegions[g][n];
                localIndex const particleIndex = neighborIndices[g][n];

                localIndex const fieldIndex = partitionField( numContactGroups,
                                                              damageFieldPartitioning,
                                                              particleGroupView[regionIndex][subRegionIndex][particleIndex],
                                                              particleDamageGradientView[regionIndex][subRegionIndex][particleIndex],
                                                              particleSurfaceNormalView[regionIndex][subRegionIndex][particleIndex],
                                                              gridDamageGradient[g] );

                if( fieldIndex != a && fieldIndex != b )
                {
                  continue;
                }

                real64 x_p[4];
                x_p[0] = particlePositionView[regionIndex][subRegionIndex][particleIndex][0] - gridPosition[g][0]; // Unsure if position of material point is relative to grid node
                x_p[1] = particlePositionView[regionIndex][subRegionIndex][particleIndex][1] - gridPosition[g][1];
                x_p[2] = particlePositionView[regionIndex][subRegionIndex][particleIndex][2] - gridPosition[g][2];
                x_p[3] = 1;
                
                real64 c_p = ( fieldIndex == a ) - ( fieldIndex == b ); // Flip this to be consistent with normal direction for average method (initial guess)

                real64 exp_x_phi = LvArray::math::exp( -LvArray::tensorOps::AiBi< 4 >( x_p, phi ) ); 
                real64 psi  = 2.0 * exp_x_phi / LvArray::math::square( 1.0 +  exp_x_phi ); 
                vec[i] += psi * w_p * ( c_p - ( 2.0 / ( 1.0 + exp_x_phi ) - 1.0 ) ) * x_p[i];
              }
            }

            // Compute increment in phi
            real64 dphi[4]; // Increment in normal and offset
            LvArray::tensorOps::Ri_eq_AijBj< 4, 4 >( dphi, matInv, vec );
            
            real64 oldPhiNorm = LvArray::tensorOps::l2Norm< 4 >( phi );

            LvArray::tensorOps::add< 4 >( phi, dphi ); // Update phi

            real64 newPhiNorm = LvArray::tensorOps::l2Norm< 4 >( phi );

            real64 n_new[3];
            n_new[0] = phi[0];
            n_new[1] = phi[1];
            n_new[2] = phi[2];
            real64 newNorm = LvArray::tensorOps::l2Norm< 3 >( n_new );
            if( newNorm > 1e-25 )
            {
              LvArray::tensorOps::scale< 3 >( n_new, 1 / newNorm );
            }
            else
            {
              GEOS_LOG_RANK("New normal had zero magnitude!");
              errored = true;
              break;
            }

            // Check convergence
            // real64 const epsilon = 1e-5; // Convergence criteria
            real64 const relativeChange = LvArray::math::abs( 1.0 - LvArray::tensorOps::AiBi< 3 >( n_new, n_old ) );
            real64 const relativeNormChange = LvArray::math::abs(newPhiNorm - oldPhiNorm)/ oldPhiNorm;
            real64 ds = ( -LvArray::math::log( 1.0 / 3.0 ) - phi[3] ) / newNorm; // Add check and constraint to make sure that surface position is never more than grid diagonal
            
            real64 s_new[3];
            LvArray::tensorOps::scaledCopy< 3 >( s_new, n_new, ds);
            
            LvArray::tensorOps::copy< 3 >( gridLRSurfacePosition[g][a], s_new );
            LvArray::tensorOps::copy< 3 >( gridLRSurfacePosition[g][b], s_new ); // Should flip direction? what do we do with multiple velocity field pairs
            LvArray::tensorOps::scaledCopy< 3 >( gridLRSurfaceNormal[g][a], n_new, -1.0 );
            LvArray::tensorOps::scaledCopy< 3 >( gridLRSurfaceNormal[g][b], n_new, 1.0 );

            real64 diff[3];
            LvArray::tensorOps::copy< 3 >( diff, s_new );
            LvArray::tensorOps::subtract< 3 >( diff, s_old );
            real64 step_change = LvArray::tensorOps::l2Norm< 3 >( diff ) / LvArray::tensorOps::l2Norm< 3 >( hEl );

             GEOS_LOG_RANK("g: " << g << ", " << 
                          "a: " << a << ", " << 
                          "b: " << b << ", " <<
                          "iter: " << iter << ", " << 
                          "phi: {" << phi[0] << ", " << phi[1] << ", " << phi[2] << ", " << phi[3] <<  "}" << ", " <<
                          "n_old: {" << n_old[0] << ", " << n_old[1] << ", " << n_old[2] << "}" << ", " <<
                          "n_new: {" << n_new[0] << ", " << n_new[1] << ", " << n_new[2] << "}" << ", " << 
                          "s_old: {" << s_old[0] << ", " << s_old[1] << ", " << s_old[2] << "}" << ", " << 
                          "s_new: {" << s_new[0] << ", " << s_new[1] << ", " << s_new[2] << "}" << ", " << 
                          "step change: " << step_change << ", " << 
                          "relChange: " << relativeChange << ", relNormChange: " << relativeNormChange  << ", ds: " << ds);
            
            if( step_change < 0.01 )
            // if( relativeChange < epsilon and relativeNormChange < 0.4 )
            {
              // Need to do contact calculations
              GEOS_LOG_RANK("Converged!");
              // LvArray::tensorOps::scaledCopy< 3 >( gridLRNormal[g][a], n_new, phi[3] );
              // LvArray::tensorOps::scaledCopy< 3 >( gridLRNormal[g][b], n_new, phi[3] ); // Should flip direction? what do we do with multiple velocity field pairs
              converged = true;
              break;
            }

            LvArray::tensorOps::copy< 3 >( n_old, n_new );
            LvArray::tensorOps::copy< 3 >( s_old, s_new );
          } 

          // Dump node info to test in MATLAB
          std::stringstream nstream;
          nstream << "g: " << g;
          std::stringstream sstream;
          sstream <<  "n0: {" << n0[0] << ", " << n0[1] << "," << n0[2] << "}\n";
          sstream << "x1, x2, x3, c\n";
          for( localIndex n=0; n < numNeighboringParticles; ++n )
          {
            localIndex const regionIndex = neighborRegions[g][n];
            localIndex const subRegionIndex = neighborSubRegions[g][n];
            localIndex const particleIndex = neighborIndices[g][n];

            localIndex const fieldIndex = partitionField( numContactGroups,
                                                          damageFieldPartitioning,
                                                          particleGroupView[regionIndex][subRegionIndex][particleIndex],
                                                          particleDamageGradientView[regionIndex][subRegionIndex][particleIndex],
                                                          particleSurfaceNormalView[regionIndex][subRegionIndex][particleIndex],
                                                          gridDamageGradient[g] );
            nstream << ", {" << n << ", " << regionIndex << ", " << subRegionIndex << ", " << particleIndex << ", " << fieldIndex << "}"  ;
            if( fieldIndex != a && fieldIndex != b )
            {
              continue;
            }

            real64 xp[3];
            LvArray::tensorOps::copy< 3 >( xp, particlePositionView[regionIndex][subRegionIndex][particleIndex] );
            LvArray::tensorOps::subtract< 3 >( xp, gridPosition[g] );
            sstream <<  xp[0] << " " << xp[1] << " " << xp[2] << " " << ((fieldIndex == a ) - ( fieldIndex == b )) << "\n";
          }
          GEOS_LOG_RANK(sstream.str());
          GEOS_LOG_RANK(nstream.str());

          if( !converged && !errored )
          {
            GEOS_LOG_RANK("Logistic regression did not converge! Using last iteration value.");
          }
        }
     } 
  } );
}


void SolidMechanicsMPM::mapSurfaceNormalsAndPositionsToParticles( ParticleManager & particleManager,
                                                                  NodeManager & nodeManager )
{
  // int const numContactGroups = m_numContactGroups;
  // int const damageFieldPartitioning = m_damageFieldPartitioning;

  // arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  // arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView3d< real64 const > const & gridLRSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfaceNormalString() );
  arrayView3d< real64 const > const & gridLRSurfacePosition = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridLRSurfacePositionString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {

    // arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    // arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    // arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    // arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();

    arrayView2d< real64 > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView2d< real64 > const particleReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
    arrayView2d< real64 > const particleReferenceSurfacePosition = subRegion.getField< fields::mpm::particleReferenceSurfacePosition >();

    // Get views to mapping arrays
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< localIndex const > const mappedFields = m_mappedFields[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];
    
    // Map to particles
    // ParticleType const particleType = subRegion.getParticleType();
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      if( particleSurfaceFlag[p] > 1 ) // Add check to determine if we want to override existing particle surface normal (e.g. normal mag > 0 )
      {
        if( LvArray::tensorOps::l2Norm< 3 >( particleSurfaceNormal[p] ) > 1e-20 && m_overwriteExistingNormalsAndPositions == 0 )
        {
          return;
        }

        LvArray::tensorOps::fill< 3 >( particleSurfaceNormal[p], 0.0 );
        LvArray::tensorOps::fill< 3 >( particleSurfacePosition[p], 0.0 );
        LvArray::tensorOps::fill< 3 >( particleReferenceSurfaceNormal[p], 0.0 );
        LvArray::tensorOps::fill< 3 >( particleReferenceSurfacePosition[p], 0.0 );

        real64 mappedMass = 0.0;
        for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
        {
          localIndex const mappedNode = mappedNodes[pp][g];
          localIndex const fieldIndex = mappedFields[pp][g];
          real64 const shapeFunctionValue = shapeFunctionValues[pp][g];
          real64 shapeFunctionGradientValue[3]; 
          LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[pp][g] );

          
          real64 norm = LvArray::tensorOps::l2Norm< 3 >( gridLRSurfaceNormal[mappedNode][fieldIndex] );
          if( norm <  1e-20 )
          {
            continue;
          }

          mappedMass += shapeFunctionValue;

          LvArray::tensorOps::scaledAdd< 3 >( particleSurfaceNormal[p], gridLRSurfaceNormal[mappedNode][fieldIndex], shapeFunctionValue );
          LvArray::tensorOps::scaledAdd< 3 >( particleSurfacePosition[p], gridLRSurfacePosition[mappedNode][fieldIndex], shapeFunctionValue );
        }

        if( mappedMass < 1e-20 )
        {
          return;
        }
        
        LvArray::tensorOps::normalize< 3 >( particleSurfaceNormal[p] ); // Should check for zero magnitude
        LvArray::tensorOps::scale< 3 >( particleSurfacePosition[p], mappedMass );

        // Update reference surface normals and positions so kinematic update will be correct based on particle deformation gradient
        real64 invF[3][3];
        LvArray::tensorOps::invert< 3 >( invF, particleDeformationGradient[p] );

        real64 cofF[3][3];
        LvArray::tensorOps::cofactor< 3 >( cofF, particleDeformationGradient[p] );

        real64 invCofF[3][3];
        LvArray::tensorOps::invert< 3 >( invCofF, cofF );

        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleReferenceSurfaceNormal[p], invCofF, particleSurfaceNormal[p] );
        LvArray::tensorOps::normalize< 3 >( particleReferenceSurfaceNormal[p] );
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleReferenceSurfacePosition[p], invF, particleSurfacePosition[p] );

        GEOS_LOG_RANK("p: " << p << ", " << 
                      "mappedMass: " << mappedMass << ", " << 
                      "pNormal: {"<< particleSurfaceNormal[p][0] << ", " << particleSurfaceNormal[p][1] << ", " << particleSurfaceNormal[p][2] << "}, " <<
                      "pPosition: {"<< particleSurfacePosition[p][0] << ", " << particleSurfacePosition[p][1] << ", " << particleSurfacePosition[p][2] << "}, " <<
                      "pRefNormal: {"<< particleReferenceSurfaceNormal[p][0] << ", " << particleReferenceSurfaceNormal[p][1] << ", " << particleReferenceSurfaceNormal[p][2] << "}, " <<
                      "pRefPosition: {"<< particleReferenceSurfacePosition[p][0] << ", " << particleReferenceSurfacePosition[p][1] << ", " << particleReferenceSurfacePosition[p][2] << "}");
      }
    } );
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::interpolateTable( real64 x,
                                          real64 dx,
                                          array2d< real64 > table,
                                          arrayView1d< real64 > output,
                                          SolidMechanicsMPM::InterpolationOption interpolationType )
{
  int numRows = table.size(0);
  int numColumns = table.size(1);
  int numDims = numColumns - 1;

  real64 xNext = x + dx;

  // If we preceed table use first row
  if( xNext <= table[0][0] )
  {
    for( localIndex i = 0; i < numDims; ++i )
    {
      output[i] = table[0][i + 1];
    }
    return;
  }

  // If we exceed table use last row
  if( xNext >= table[numRows - 1][0] )
  {
    for( localIndex i = 0; i < numDims; ++i )
    {
      output[i] = table[numRows - 1][i + 1];
    }
    return;
  }

  int tableInterval = 0;
  for( localIndex i = 0; i < numRows; ++i )
  {
    if( x + dx / 2 > table[i][0] )
    {
      tableInterval = i;
    }
  }

  real64 timeInterval = table[tableInterval + 1][0] - table[tableInterval][0]; // Time fInterval for current part of F table we're in
  real64 timePast = xNext - table[tableInterval][0]; // Time elapsed since switching intervals in F table
  real64 timeFrac = timePast / timeInterval;

  for( localIndex i = 0; i < numDims; ++i )
  {
      switch( interpolationType )
      {
        case SolidMechanicsMPM::InterpolationOption::Linear:
          // default linear interpolation
          output[i] = table[tableInterval][i + 1] * ( 1.0 - timeFrac ) + table[tableInterval + 1][i + 1] * timeFrac;
          break;
        case SolidMechanicsMPM::InterpolationOption::Cosine:
          // smooth-step interpolation with cosine, zero endpoint velocity
          output[i] = table[tableInterval][i + 1] - 0.5 * ( table[tableInterval + 1][i + 1] - table[tableInterval][i + 1] ) * ( cos( 3.141592653589793 * timeFrac ) - 1.0 );
          break;
        case SolidMechanicsMPM::InterpolationOption::Smoothstep:
          // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
          output[i] = table[tableInterval][i+1] + ( table[tableInterval+1][i+1] - table[tableInterval][i+1] ) *
                      ( 10.0 * pow( timeFrac, 3 ) - 15.0 * pow( timeFrac, 4 ) + 6.0 * pow( timeFrac, 5 ) );
          break;
        default:
          GEOS_ERROR( "No interpolation option of that type!" );
          break;
      }
  }
}

void SolidMechanicsMPM::interpolateValueInRange( real64 const & x,
                                                 real64 const & xmin,
                                                 real64 const & xmax,
                                                 real64 const & ymin,
                                                 real64 const & ymax,
                                                 real64 & output,
                                                 SolidMechanicsMPM::InterpolationOption interpolationType )
{ // Stripped down version of the table interpolation to be used when interpolating a 1D function
  // in a fixed and well-defined range (xmin<x<xmax).
  // Perhaps this should construct a table and call the tableInterpolation so we don't need to
  // keep both consistent if other inteprolation types are created.
  GEOS_MARK_FUNCTION;

  if( ( x < xmin ) || ( xmax <= xmin ) )
  {
   output = ymin;
  }
  else if ( x>xmax )
  {
    output = ymax;
  }
  else
  {

  real64 timeFrac = (x-xmin)/(xmax-xmin);
      switch( interpolationType )
      {
        case SolidMechanicsMPM::InterpolationOption::Linear:
          // default linear interpolation
          output = ymin * ( 1.0 - timeFrac ) + ymax * timeFrac;
          break;
        case SolidMechanicsMPM::InterpolationOption::Cosine:
          // smooth-step interpolation with cosine, zero endpoint velocity
          output = ymin - 0.5 * ( ymax - ymin ) * ( cos( 3.141592653589793 * timeFrac ) - 1.0 );
          break;
        case SolidMechanicsMPM::InterpolationOption::Smoothstep:
          // smooth-step interpolation with 5th order polynomial, zero endpoint velocity and acceleration
          output = ymin + ( ymax - ymin ) * ( 10.0 * pow( timeFrac, 3 ) - 15.0 * pow( timeFrac, 4 ) + 6.0 * pow( timeFrac, 5 ) );
          break;
        default:
          GEOS_ERROR( "No interpolation option of that type!" );
          break;
      }
  }
}

void SolidMechanicsMPM::interpolateFTable( real64 dt, real64 time_n )
{
  GEOS_MARK_FUNCTION;

  array1d< real64 > Fii_new(3);
  interpolateTable( time_n,
                    dt,
                    m_fTable,
                    Fii_new,
                    m_fTableInterpType );

  for( localIndex i = 0; i < m_numDims; ++i )
  {
    if( m_stressControl[i] != 1 )
    {
      real64 Fii_dot = ( Fii_new[i] - m_domainF[i] ) / dt;
      m_domainL[i] = Fii_dot / Fii_new[i]; // L = Fdot.Finv
      m_domainF[i] = Fii_new[i];
    }
  }
}

void SolidMechanicsMPM::interpolateStressTable( real64 dt,
                                                real64 time_n )
{
  GEOS_MARK_FUNCTION;

  interpolateTable( time_n,
                    dt,
                    m_stressTable,
                    m_domainStress,
                    m_fTableInterpType );
}

void SolidMechanicsMPM::gridToParticle( real64 dt,
                                        ParticleManager & particleManager,
                                        NodeManager & nodeManager,
                                        DomainPartition & domain,
                                        MeshLevel & mesh )
{
  GEOS_MARK_FUNCTION;

  // Grid-to-particle map
  switch(( SolidMechanicsMPM::UpdateMethodOption ) m_updateMethod)
  {
    case SolidMechanicsMPM::UpdateMethodOption::FLIP:
      performFLIPUpdate( dt,
                         particleManager,
                         nodeManager );
      break;
    case SolidMechanicsMPM::UpdateMethodOption::PIC:
      performPICUpdate( dt,
                        particleManager,
                        nodeManager );
      break;
    case SolidMechanicsMPM::UpdateMethodOption::XPIC:
      performXPICUpdate( dt,
                         particleManager,
                         nodeManager,
                         domain,
                         mesh );
      break;
    case SolidMechanicsMPM::UpdateMethodOption::FMPM:
      performFMPMUpdate( dt,
                         particleManager,
                         nodeManager,
                         domain,
                         mesh );
      break;
    default:
      GEOS_ERROR("SolidMechanicsMPM solver update method not recognized.");
      break;
  }
}

void SolidMechanicsMPM::performFLIPUpdate( real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager )
{
  // On-the-fly grid shape function computations
  real64 xLocalMin[3] = {};
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {

    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    // Registered by MPM solver
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    // arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    // arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    // arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const numContactGroups = m_numContactGroups;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    ParticleType const particleType = subRegion.getParticleType();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Initialize velocity gradient for this particle
      // if( m_directionalOverlapCorrection == 0 )
      // {
      for( localIndex i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }
      // }

      // On-the-fly grid shape function computations
      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues );

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        // localIndex const mappedNode = mappedNodes[pp][g];
        localIndex const mappedNode = mappedNodes[g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        for( localIndex i=0; i<numDims; ++i )
        {
          // // Full step update
          // // particlePosition[p][i] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionValues[pp][g] * dt;

          // // Half step update ( equivalent to old geos, but if there is only one cell in the simulation that is anomalous issues with the particle velocity gradient and by extension the stresses)
          // // If the simulation size is increased to 2x2x2 then the issue goes away and likewise if the rvectors are scaled to 80% for a single cell simulation
          // particlePosition[p][i] += ( gridVelocity[mappedNode][fieldIndex][i] - 0.5 * gridAcceleration[mappedNode][fieldIndex][i] * dt ) * shapeFunctionValues[pp][g] * dt;

          // particleVelocity[p][i] += gridAcceleration[mappedNode][fieldIndex][i] * dt * shapeFunctionValues[pp][g]; // FLIP
          // for( int j=0; j<numDims; ++j )
          // {
          //   particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j]; // Technically
          //                                                                                                                         // wrong,
          //                                                                                                                         // the
          //                                                                                                                         // best
          //                                                                                                                         // kind of
          //                                                                                                                         // wrong
          //                                                                                                                         // (end-of-step
          //                                                                                                                         // velocity
          //                                                                                                                         // with
          //                                                                                                                         // beginning-of-step
          //                                                                                                                         // gradient)
          // }

          // Half step update ( equivalent to old geos, but if there is only one cell in the simulation that is anomalous issues with the particle velocity gradient and by extension the stresses)
          // If the simulation size is increased to 2x2x2 then the issue goes away and likewise if the rvectors are scaled to 80% for a single cell simulation
          particlePosition[p][i] += ( gridVelocity[mappedNode][fieldIndex][i] - 0.5 * gridAcceleration[mappedNode][fieldIndex][i] * dt ) * shapeFunctionValues[g] * dt;

          particleVelocity[p][i] += gridAcceleration[mappedNode][fieldIndex][i] * dt * shapeFunctionValues[g]; // FLIP
          for( int j=0; j<numDims; ++j )
          {
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[g][j];
          }
        }
      }
    } );
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::performPICUpdate(  real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager )
{
  // On-the-fly grid shape function computations
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    // Registered by MPM solver
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    // Map to particles
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    int const numDims = m_numDims;
    int const damageFieldPartitioning = m_damageFieldPartitioning;
    int const numContactGroups = m_numContactGroups;
    ParticleType const particleType = subRegion.getParticleType();

    #ifndef GEOS_USE_DEVICE
      arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
      arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
      arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];
      GEOS_UNUSED_VAR( particleType );
      GEOS_UNUSED_VAR( particleRVectors );
      GEOS_UNUSED_VAR( gridPosition );
      GEOS_UNUSED_VAR( ijkMap );
    #endif
    
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      #ifdef GEOS_USE_DEVICE
        // On-the-fly grid shape function computations
        localIndex mappedNodes[64];
        real64 shapeFunctionValues[64];
        real64 shapeFunctionGradientValues[64][3];
        mapNodesAndComputeShapeFunctions( ijkMap,
                                          xLocalMin,
                                          hEl,
                                          particleType,
                                          particlePosition[p],
                                          particleRVectors[p],
                                          gridPosition,
                                          mappedNodes,
                                          shapeFunctionValues,
                                          shapeFunctionGradientValues );
      #endif                                  

      for( localIndex i = 0; i < 3; ++i )
      {
        particleVelocity[p][i] = 0.0;
        for( localIndex j = 0; j < 3; ++j )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        #ifdef GEOS_USE_DEVICE
          localIndex const mappedNode = mappedNodes[g];
          real64 const shapeFunctionValue = shapeFunctionValues[g];
          real64 shapeFunctionGradientValue[3]; 
          LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[g] );
        #else
          localIndex const mappedNode = mappedNodes[pp][g];
          real64 const shapeFunctionValue = shapeFunctionValues[pp][g];
          real64 shapeFunctionGradientValue[3]; 
          LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[pp][g] );
        #endif

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        for( localIndex i = 0; i < numDims; ++i )
        {
          // particlePosition[p][i] += ( gridVelocity[mappedNode][fieldIndex][i] - 0.5 * gridAcceleration[mappedNode][fieldIndex][i] * dt ) * shapeFunctionValue[g] * dt; // CC: position update doesn't seem consistent with old GEOS for FLIP and PIC, need to double check
          // particleVelocity[p][i] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionValue[g];
          // for( int j=0; j<numDims; ++j )
          // {
          //   particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValue[j];
          // }

          particlePosition[p][i] += ( gridVelocity[mappedNode][fieldIndex][i] - 0.5 * gridAcceleration[mappedNode][fieldIndex][i] * dt ) * shapeFunctionValue * dt; // CC: position update doesn't seem consistent with old GEOS for FLIP and PIC, need to double check
          particleVelocity[p][i] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionValue;
          for( int j=0; j<numDims; ++j )
          {
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValue[j];
          }
        }
      }
    } );
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::performXPICUpdate( real64 dt,
                                           ParticleManager & particleManager,
                                           NodeManager & nodeManager,
                                           DomainPartition & domain,
                                           MeshLevel & mesh )
{
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  real64 xLocalMax[3];
  LvArray::tensorOps::copy< 3 >( xLocalMax, m_xLocalMax);

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  // arrayView3d< real64 const > const & gridDVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVelocityString() ); // for multifield contact corrections
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 const > const & gridAcceleration = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridAccelerationString() );

  arrayView3d< real64 > const & gridVPlus = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVPlusString() );
  arrayView3d< real64 > const & gridDVPlus = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridDVPlusString() );

  int numNodes = nodeManager.size();
  int const numDims = m_numDims;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numContactGroups = m_numContactGroups;
  int const numVelocityFields = m_numVelocityFields;
  int const updateOrder = m_updateOrder;

  // For iterative XPIC solve
  array3d< real64 > vStar( numNodes, numVelocityFields, numDims );
  array3d< real64 > vMinus( numNodes, numVelocityFields, numDims );
  array3d< real64 > dVMinus( numNodes, numVelocityFields, numDims );

  // Zero out vStar for each order iteration
  for( int n = 0; n < numNodes; ++n)
  {
    for( int cg = 0; cg < numVelocityFields; ++cg)
    {
      for( localIndex i = 0; i < numDims; ++i)
      {
        vStar[n][cg][i] = 0.0;
        vMinus[n][cg][i] = gridVelocity[n][cg][i] - gridAcceleration[n][cg][i] * dt;

        gridDVPlus[n][cg][i] = 0.0; // gridDVelocity[n][cg][i]; // CC: this isn't working, currently disabled by writing 0 to it
        dVMinus[n][cg][i] = gridDVPlus[n][cg][i];
      }
    }
  }

  // Do XPIC iterations
  for(int r = 2; r <= updateOrder; ++r)
  {
    // Zero out vPlus for each order iteration
    for( int n = 0; n < numNodes; ++n)
    {
      for( int cg = 0; cg < numVelocityFields; ++cg)
      {
        for( localIndex i = 0; i < numDims; ++i)
        {
          gridVPlus[n][cg][i] = 0.0;
        }
      }
    }

    localIndex subRegionIndex = 0;
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      // Registered by subregion
      arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
      // arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
      arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
      arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
      arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

      // Registered by MPM solver
      arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
      arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
      arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

      // Get views to mapping arrays
      int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
      // arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
      // arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
      // arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

      ParticleType const particleType = subRegion.getParticleType();

      // Map to particles
      // Seems like LvArrays that aren't views need to be passed by reference to forAll loops (e.g. vPlus)
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        // On-the-fly shape function computation
        localIndex mappedNodes[64];
        real64 shapeFunctionValues[64];
        real64 shapeFunctionGradientValues[64][3];
        mapNodesAndComputeShapeFunctions( ijkMap,
                                          xLocalMin,
                                          hEl,
                                          particleType,
                                          particlePosition[p],
                                          particleRVectors[p],
                                          gridPosition,
                                          mappedNodes,
                                          shapeFunctionValues,
                                          shapeFunctionGradientValues );

        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; ++gi )
        {
          localIndex const mappedNodeI = mappedNodes[gi];
          // int const nodeFlagI = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeI], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
          // localIndex const fieldIndexI = nodeFlagI * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
          localIndex const fieldIndexI = partitionField( numContactGroups,
                                                  damageFieldPartitioning,
                                                  particleGroup[p],
                                                  particleDamageGradient[p],
                                                  particleSurfaceNormal[p],
                                                  gridDamageGradient[mappedNodeI] );


          if( gridMass[mappedNodeI][fieldIndexI] > m_smallMass )
          {
            for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; ++gj )
            {
              localIndex const mappedNodeJ = mappedNodes[gj];
              // int const nodeFlagJ = ( damageFieldPartitioning == 1 && LvArray::tensorOps::AiBi< 3 >( gridDamageGradient[mappedNodeJ], particleDamageGradient[p] ) < 0.0 ) ? 1 : 0;
              // localIndex const fieldIndexJ = nodeFlagJ * numContactGroups + particleGroup[p]; // This ranges from 0 to nMatFields-1
              localIndex const fieldIndexJ = partitionField( numContactGroups,
                                                      damageFieldPartitioning,
                                                      particleGroup[p],
                                                      particleDamageGradient[p],
                                                      particleSurfaceNormal[p],
                                                      gridDamageGradient[mappedNodeJ] );

              for (localIndex i = 0; i < numDims; ++i)
              {
                gridVPlus[mappedNodeI][fieldIndexI][i] += ( ( updateOrder - r + 1.0 ) / r ) *
                                                          ( particleMass[p] * shapeFunctionValues[gi] * shapeFunctionValues[gj] / gridMass[mappedNodeI][fieldIndexI] ) *
                                                            vMinus[mappedNodeJ][fieldIndexJ][i];
              }
            }

          }
        }
      } );
      ++subRegionIndex;
    } );

    // syncGridFields( { viewKeyStruct::gridVPlusString() }, domain, nodeManager, mesh, MPI_SUM );
    syncGridFields( { viewKeyStruct::gridVPlusString(), viewKeyStruct::gridDVPlusString() }, domain, nodeManager, mesh, MPI_SUM ); // Also need to sync dVPlus later for multifield contact correction

    // Update vStar
    for( int n=0; n < numNodes; ++n)
    {
      for( int cg=0; cg < numVelocityFields; ++cg)
      {
        for( localIndex i = 0; i < numDims; ++i)
        {
          vStar[n][cg][i] += LvArray::math::pow(-1.0, r) * ( gridVPlus[n][cg][i] ); // + ( updateOrder - 1.0 ) / updateOrder * gridDVPlus[n][cg][i] );
          vMinus[n][cg][i] = gridVPlus[n][cg][i];
        }
      }
    }

    // Add if statement? because this doesn't need to be performed for the last updateOrder
    // Zero out vPlus for each order iteration
    for( int n=0; n < numNodes; ++n)
    {
      for( int cg=0; cg < numVelocityFields; ++cg)
      {
        for( localIndex i = 0; i < numDims; ++i)
        {
          gridDVPlus[n][cg][i] = 0.0;
        }
      }
    }

    subRegionIndex = 0;
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      // Registered by subregion
      arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
      // arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
      arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
      arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
      arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

      // Registered by MPM solver
      arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
      arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
      arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

      // Get views to mapping arrays
      int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
      ParticleType const particleType = subRegion.getParticleType();

      // Update dVStar
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        // On-the-fly shape function computation
        localIndex mappedNodes[64];
        real64 shapeFunctionValues[64];
        real64 shapeFunctionGradientValues[64][3];
        mapNodesAndComputeShapeFunctions( ijkMap,
                                          xLocalMin,
                                          hEl,
                                          particleType,
                                          particlePosition[p],
                                          particleRVectors[p],
                                          gridPosition,
                                          mappedNodes,
                                          shapeFunctionValues,
                                          shapeFunctionGradientValues );

        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; ++gi )
        {
          localIndex const mappedNodeI = mappedNodes[gi];
          localIndex const fieldIndexI = partitionField( numContactGroups,
                                                  damageFieldPartitioning,
                                                  particleGroup[p],
                                                  particleDamageGradient[p],
                                                  particleSurfaceNormal[p],
                                                  gridDamageGradient[mappedNodeI] );

          if( gridMass[mappedNodeI][fieldIndexI] > m_smallMass )
          {
            for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; ++gj )
            {
              localIndex const mappedNodeJ = mappedNodes[gj];
              localIndex const fieldIndexJ = partitionField( numContactGroups,
                                                      damageFieldPartitioning,
                                                      particleGroup[p],
                                                      particleDamageGradient[p],
                                                      particleSurfaceNormal[p],
                                                      gridDamageGradient[mappedNodeJ] );


              for (localIndex i = 0; i < numDims; ++i)
              {
                gridDVPlus[mappedNodeI][fieldIndexI][i] += ( ( updateOrder - r ) / r ) *
                                                        ( particleMass[p] * shapeFunctionValues[gi] * shapeFunctionValues[gj] / gridMass[mappedNodeI][fieldIndexI] ) *
                                                          dVMinus[mappedNodeJ][fieldIndexJ][i];
              }
            }
          }
        }
      } );
      ++subRegionIndex;
    } );

    for( int n=0; n < numNodes; ++n)
    {
      for( int cg=0; cg < numVelocityFields; cg++)
      {
        for( localIndex i = 0; i < numDims; ++i)
        {
          dVMinus[n][cg][i] = gridDVPlus[n][cg][i];
        }
      }
    }
  } //End of updateOrder iterations


  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();

    // Registered by MPM solver
    arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Update particles position and velocities now
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Zero velocity gradient
      for( localIndex i=0; i < numDims; ++i )
      {
        particlePosition[p][i] -= particleVelocity[p][i] * dt / 2.0;
        particleVelocity[p][i] = 0.0;
        for( int j=0; j < numDims; ++j )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        localIndex const mappedNode = mappedNodes[pp][g];

        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        for( localIndex i=0; i<numDims; ++i )
        {
          real64 m = updateOrder;
          real64 S = shapeFunctionValues[pp][g];
          real64 gVPlus = gridVelocity[mappedNode][fieldIndex][i];
          real64 gA = gridAcceleration[mappedNode][fieldIndex][i];

          particlePosition[p][i] += S * gVPlus * dt  - ( S * gA * dt - m * S *( gVPlus - gA * dt ) + m * S * vStar[mappedNode][fieldIndex][i] ) * dt / 2.0;

          // particlePosition[p][i] += S * ( gVPlus * dt  - ( (1 + m ) * gA * dt + m * ( vStar[mappedNode][fieldIndex][i] - gVPlus ) ) * dt / 2.0 );

          particleVelocity[p][i] += S * ( m * ( gVPlus - vStar[mappedNode][fieldIndex][i] ) + ( 1 - m ) * gA * dt );

          // CC: What about update to velocity gradient?
          // Currently copy this from FLIP udpate with change from gridVelocity to vStar
          for( int j=0; j < numDims; ++j )
          {
            particleVelocityGradient[p][i][j] += gridVelocity[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[pp][g][j];
          }
        }
      }

    } );
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::performFMPMUpdate(  real64 dt,
                                            ParticleManager & particleManager,
                                            NodeManager & nodeManager,
                                            DomainPartition & domain,
                                            MeshLevel & mesh )
{
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  real64 xLocalMax[3];
  LvArray::tensorOps::copy< 3 >( xLocalMax, m_xLocalMax);
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  // Grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 const > const & gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );
  arrayView3d< real64 const > const & gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );
  arrayView3d< real64 > const & gridVPlus = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVPlusString() );

  int numNodes = nodeManager.size();
  int const numDims = m_numDims;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numContactGroups = m_numContactGroups;
  int const numVelocityFields = m_numVelocityFields;
  int const updateOrder = m_updateOrder;

  // Added these as fields to nodeManager to easily sync them in parallelization for each iteration (technically only vPlus needs to be synced )
  // For iterative FMPM solve
  array3d< real64 > vStar( numNodes, numVelocityFields, numDims );
  array3d< real64 > vMinus( numNodes, numVelocityFields, numDims );

  // Initialize FMPM variables
  for( int n=0; n < numNodes; ++n)
  {
    for( int cg=0; cg < numVelocityFields; cg++)
    {
      for( localIndex i = 0; i < numDims; ++i)
      {
        vStar[n][cg][i] = updateOrder * gridVelocity[n][cg][i];
        vMinus[n][cg][i] = vStar[n][cg][i];
      }
    }
  }

  // Perform FMPM order iterations
  for(int r=2; r <= updateOrder; ++r)
  {
    // Zero out vPlus for each order iteration
    for( int n=0; n < numNodes; ++n)
    {
      for( int cg=0; cg < numVelocityFields; cg++)
      {
        for( localIndex i = 0; i < numDims; ++i)
        {
          gridVPlus[n][cg][i] = 0.0;
        }
      }
    }

    localIndex subRegionIndex = 0;
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      // Registered by subregion
      arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
      arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
      arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
      arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

      arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
      arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

      // Get views to mapping arrays
      int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

      ParticleType const particleType = subRegion.getParticleType();

      // Map to particles
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        // On-the-fly shape function computation
        localIndex mappedNodes[64];
        real64 shapeFunctionValues[64];
        real64 shapeFunctionGradientValues[64][3];
        mapNodesAndComputeShapeFunctions( ijkMap,
                                          xLocalMin,
                                          hEl,
                                          particleType,
                                          particlePosition[p],
                                          particleRVectors[p],
                                          gridPosition,
                                          mappedNodes,
                                          shapeFunctionValues,
                                          shapeFunctionGradientValues );

        for( int gi = 0; gi < 8 * numberOfVerticesPerParticle; ++gi )
        {
          localIndex const mappedNodeI = mappedNodes[gi];
          localIndex const fieldIndexI = partitionField( numContactGroups,
                                                  damageFieldPartitioning,
                                                  particleGroup[p],
                                                  particleDamageGradient[p],
                                                  particleSurfaceNormal[p],
                                                  gridDamageGradient[mappedNodeI] );

          if( gridMass[mappedNodeI][fieldIndexI] > m_smallMass )
          {
            real64 Splus = particleMass[p] * shapeFunctionValues[gi] / gridMass[mappedNodeI][fieldIndexI];

            for(int gj = 0; gj < 8 * numberOfVerticesPerParticle; ++gj )
            {
              localIndex const mappedNodeJ = mappedNodes[gj];
              localIndex const fieldIndexJ = partitionField( numContactGroups,
                                                      damageFieldPartitioning,
                                                      particleGroup[p],
                                                      particleDamageGradient[p],
                                                      particleSurfaceNormal[p],
                                                      gridDamageGradient[mappedNodeJ] );

              for (localIndex i = 0; i < numDims; ++i)
              {
                gridVPlus[mappedNodeI][fieldIndexI][i] += Splus * shapeFunctionValues[gj] * vMinus[mappedNodeJ][fieldIndexJ][i];;
              }
            }
          }
        }
      } );
      ++subRegionIndex;
    } );

    syncGridFields( { viewKeyStruct::gridVPlusString() }, domain, nodeManager, mesh, MPI_SUM );

    // Update vStar
    real64 orderCoefficient = LvArray::math::pow(-1.0, 1.0+r) * ( updateOrder - r + 1.0 ) / r;
    for( int n=0; n < numNodes; ++n)
    {
      for( int cg=0; cg < numVelocityFields; cg++)
      {
        for( localIndex i = 0; i < numDims; ++i)
        {
          vStar[n][cg][i] += orderCoefficient * gridVPlus[n][cg][i];
          vMinus[n][cg][i] = gridVPlus[n][cg][i];
        }
      }
    }
  } //End of updateOrder iterations

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registered by subregion
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();

    // arrayView1d< real64 > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    // Update particles positions and velocities now
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    ParticleType const particleType = subRegion.getParticleType();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // On-the-fly shape function computation
      localIndex mappedNodes[64];
      real64 shapeFunctionValues[64];
      real64 shapeFunctionGradientValues[64][3];
      mapNodesAndComputeShapeFunctions( ijkMap,
                                        xLocalMin,
                                        hEl,
                                        particleType,
                                        particlePosition[p],
                                        particleRVectors[p],
                                        gridPosition,
                                        mappedNodes,
                                        shapeFunctionValues,
                                        shapeFunctionGradientValues );

      for(localIndex i = 0; i < numDims; ++i)
      {
        particlePosition[p][i] += 0.5 * dt * particleVelocity[p][i];
        particleVelocity[p][i] = 0.0;

        for( int j=0; j < numDims; ++j )
        {
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        localIndex const mappedNode = mappedNodes[g];
        localIndex const fieldIndex = partitionField( numContactGroups,
                                               damageFieldPartitioning,
                                               particleGroup[p],
                                               particleDamageGradient[p],
                                               particleSurfaceNormal[p],
                                               gridDamageGradient[mappedNode] );

        for( localIndex i=0; i < numDims; ++i )
        {
          particlePosition[p][i] += 0.5 * dt * shapeFunctionValues[g] * vStar[mappedNode][fieldIndex][i];
          particleVelocity[p][i] += shapeFunctionValues[g] * vStar[mappedNode][fieldIndex][i];

          for( int j=0; j < numDims; ++j )
          {
            particleVelocityGradient[p][i][j] += vStar[mappedNode][fieldIndex][i] * shapeFunctionGradientValues[g][j];
          }
        }
      }
    } );
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::updateSolverDependencies( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // Get particle damage values from constitutive model
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // // Particle fields
    // arrayView1d< real64 > const particleDamage = subRegion.getParticleDamage();
    // arrayView2d< real64 > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >();
    // arrayView1d< real64 > const particleTemperature = subRegion.getParticleTemperature();
    // arrayView1d< real64 > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();

    // Get constitutive model reference
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // arrayView2d< real64 const > constitutiveDamage;
    // arrayView3d< real64 const > constitutivePlasticStrain;
    // arrayView1d< real64 const > constitutiveTemperature;
    // arrayView2d< real64 const > constitutiveWavespeed;

    // int hasDamage = constitutiveModel.hasWrapper( "damage" );
    // int hasPlasticStrain = constitutiveModel.hasWrapper( "plasticStrain" );
    // int hasTemperature = constitutiveModel.hasWrapper( "temperature" );
    // int hasWavespeed = constitutiveModel.hasWrapper( "wavespeed" );


    // if( hasDamage )
    // {
    //   constitutiveDamage = constitutiveModel.getReference< array2d< real64 > >( "damage" );
    // }

    // if( hasPlasticStrain )
    // {
    //   constitutivePlasticStrain = constitutiveModel.getReference< array3d< real64 > >( "plasticStrain" );
    // }

    // if( hasTemperature )
    // {
    //   constitutiveTemperature = constitutiveModel.getReference< array1d< real64 > >( "temperature" );
    // }

    // if( hasWavespeed )
    // {
    //   constitutiveWavespeed = constitutiveModel.getReference< array2d< real64 > >( "wavespeed" );
    // }

    // forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    // {
    //   localIndex const p = activeParticleIndices[pp];

    //   if( hasDamage && constitutiveDamage[p][0] > particleDamage[p] ) // Damage can only increase - This will also preserve user-specified damage at
    //                                                      // initialization
    //   {
    //     particleDamage[p] = constitutiveDamage[p][0]; // TODO: Load any pre-damage into the constitutive model at initialization. Or,
    //                                                   // switch to using VTK input such that we can initialize any field we want without
    //                                                   // explicitly post-processing the input file.
    //   }

    //   if( hasPlasticStrain )
    //   {
    //     LvArray::tensorOps::copy< 6 >( particlePlasticStrain[p], constitutivePlasticStrain[p][0] );
    //   }

    //   if( hasTemperature )
    //   {
    //     particleTemperature[p] = constitutiveTemperature[p];
    //   }

    //   if( hasWavespeed )
    //   {
    //     particleWavespeed[p] = constitutiveWavespeed[p][0];
    //   }
    // } );

    if( constitutiveModel.hasWrapper( "damage" ) ) // Fragile code because someone could change the damage key without our knowledge. TODO: Make an
                                            // integrated test that checks this
    {
      arrayView1d< real64 > const particleDamage = subRegion.getParticleDamage();
      arrayView2d< real64 const > const constitutiveDamage = constitutiveModel.getReference< array2d< real64 > >( "damage" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        if( constitutiveDamage[p][0] > particleDamage[p] ) // Damage can only increase - This will also preserve user-specified damage at
                                                           // initialization
        {
          particleDamage[p] = constitutiveDamage[p][0]; // TODO: Load any pre-damage into the constitutive model at initialization. Or,
                                                        // switch to using VTK input such that we can initialize any field we want without
                                                        // explicitly post-processing the input file.
        }
      } );
    }

    if( constitutiveModel.hasWrapper( "plasticStrain" ) ) // Fragile code because someone could change the damage key without our knowledge. TODO: Make an
                                            // integrated test that checks this
    {
      arrayView2d< real64 > const particlePlasticStrain = subRegion.getField< fields::mpm::particlePlasticStrain >();
      arrayView3d< real64 const > const constitutivePlasticStrain = constitutiveModel.getReference< array3d< real64 > >( "plasticStrain" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        LvArray::tensorOps::copy< 6 >( particlePlasticStrain[p], constitutivePlasticStrain[p][0] );
      } );
    }

    if( constitutiveModel.hasWrapper( "temperature" ) )
    {
      arrayView1d< real64 > const particleTemperature = subRegion.getParticleTemperature();
      arrayView1d< real64 const > const constitutiveTemperature = constitutiveModel.getReference< array1d< real64 > >( "temperature" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        particleTemperature[p] = constitutiveTemperature[p];
      } );
    }


    if( constitutiveModel.hasWrapper( "wavespeed" ) )
    {
      arrayView1d< real64 > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();
      arrayView2d< real64 const > const constitutiveWavespeed = constitutiveModel.getReference< array2d< real64 > >( "wavespeed" );
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        particleWavespeed[p] = constitutiveWavespeed[p][0];
      } );
    }
  } );
}

// Does this require a MPI sync across partitions
// Seems to have run fine which would indicate the returned time step gets checked in a larger function call outside of solver
real64 SolidMechanicsMPM::getStableTimeStep( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxWavespeed( 0.0 );
  real64 length = m_planeStrain == 1 ? LvArray::math::min( m_hEl[0], m_hEl[1] ) : LvArray::math::min( m_hEl[0], LvArray::math::min( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();


    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp ) // would need reduction to parallelize
    {
      localIndex const p = activeParticleIndices[pp];
      maxWavespeed.max( particleWavespeed[p] + LvArray::tensorOps::l2Norm< 3 >( particleVelocity[p] ) );
    } );
  } );

  real64 dtReturn = maxWavespeed.get() > 1.0e-16 ? m_cflFactor * length / maxWavespeed.get() : DBL_MAX; // This partitions's dt, make it huge if wavespeed=0.0
                                                                                      // (this happens when there are no particles on this
                                                                                      // partition)

  // CC: TODO add check to make sure that next time step won't skip a MPM event
  return dtReturn;
}

void SolidMechanicsMPM::deleteBadParticles( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int size = 0;
  // Cases covered:
  // 1.) Particles that map outside the domain (including buffer cells)
  // 2.) Particles with unacceptable Jacobian (<0.1 or >10)
  // 3.) Particles that are removed by the machine sample event
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Move everything into the host memory space
    subRegion.forWrappers( [&]( WrapperBase & wrapper )
    {
      wrapper.move( LvArray::MemorySpace::host, true );
    } );

    // Get relevant particle arrays
    arrayView1d< int > const particleDeleteFlag = subRegion.getField< fields::mpm::particleDeleteFlag >();

    // Initialize the set of particles to delete
    std::set< localIndex > indicesToErase;
    forAll< serialPolicy >( subRegion.size(), [=, &indicesToErase] GEOS_HOST ( localIndex const p )
    {
      if( particleDeleteFlag[p] == 1 )
      {
        indicesToErase.insert( p );
      }
    } );
    subRegion.erase( indicesToErase );
    subRegion.setActiveParticleIndices();
    size = subRegion.activeParticleIndices().size();
  } );

  // particleManager.getRegion< ParticleRegion >( "ParticleRegion1" ).resize( size ) ;
}

void SolidMechanicsMPM::printProfilingResults()
{
  // Use MPI reduction to get the average elapsed time for each step on all partitions
  int rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  unsigned int numIntervals = m_profilingTimes.size() - 1;
  stdVector< real64 > timeIntervalsAllRanks( numIntervals );

  // Get total CPU time for the entire time step
  real64 totalStepTimeThisRank = m_profilingTimes[numIntervals] - m_profilingTimes[0];
  real64 totalStepTimeAllRanks;
  totalStepTimeAllRanks =  MpiWrapper::allReduce( totalStepTimeThisRank,
                                                  MpiWrapper::Reduction::Sum,
                                                  MPI_COMM_GEOS );

  // Get total CPU times for each queried time interval
  for( unsigned int i = 0; i < numIntervals; ++i )
  {
    real64 timeIntervalThisRank = ( m_profilingTimes[i+1] - m_profilingTimes[i] );
    real64 timeIntervalAllRanks = MpiWrapper::allReduce( timeIntervalThisRank,
                                                         MpiWrapper::Reduction::Sum,
                                                         MPI_COMM_GEOS );
    if( rank == 0 )
    {
      timeIntervalsAllRanks[i] = timeIntervalAllRanks;
    }
  }

  // Print out solver profiling
  MPI_Barrier( MPI_COMM_GEOS );
  if( rank == 0 )
  {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Fraction of total time for one step: " << std::endl;
    for( unsigned int i = 0; i < numIntervals; ++i )
    {
      std::cout << " (" << i << ") ";
      std::cout << std::fixed;
      std::cout << std::showpoint;
      std::cout << std::setprecision( 6 ) << timeIntervalsAllRanks[i] / totalStepTimeAllRanks;
      std::cout << ": " << m_profilingLabels[i] << std::endl;
    }
    std::cout << " ** Total step CPU time:  " << totalStepTimeAllRanks << " s **" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
  }

  // Reset profiling arrays
  m_profilingTimes.resize( 0 );
  m_profilingLabels.resize( 0 );
}

void SolidMechanicsMPM::computeSurfaceFlags( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  localIndex const numDims = m_numDims;
  real64 const neighborRadius = m_neighborRadius;

  // Get accessors for volume, position
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );

  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > particleVolumeView = particleVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePositionView = particlePositionAccessor.toNestedViewConst();

  // Perform neighbor operations
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and surface flags
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< int > const particleSurfaceFlag = subRegion.getField< fields::mpm::particleSurfaceFlag >();

    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Get number of neighbors and accessor indices
      localIndex numNeighbors = numNeighborsAll[p];
      arraySlice1d< localIndex const > const regionIndices = neighborRegions[pp];
      arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[pp];
      arraySlice1d< localIndex const > const particleIndices = neighborIndices[pp];

      // Calculate
      real64 totalVolume = 0.0;
      real64 centroid[3];
      for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
      {
        localIndex regionIndex = regionIndices[neighborIndex];
        localIndex subRegionIndex = subRegionIndices[neighborIndex];
        localIndex particleIndex = particleIndices[neighborIndex];
        real64 rSquared = 0.0;
        for( localIndex i = 0; i < numDims; ++i )
        {
          real64 relativePositionComponent = particlePositionView[regionIndex][subRegionIndex][particleIndex][i] - particlePosition[p][i];
          rSquared += relativePositionComponent * relativePositionComponent;
        }
        real64 kernelVal = kernel( LvArray::math::sqrt( rSquared ));
        totalVolume += kernelVal * particleVolumeView[regionIndex][subRegionIndex][particleIndex];
        for( localIndex i = 0; i < numDims; ++i )
        {
          centroid[i] += kernelVal * particleVolumeView[regionIndex][subRegionIndex][particleIndex] * particlePositionView[regionIndex][subRegionIndex][particleIndex][i];
        }
      }
      real64 positionDifference[3];
      for( localIndex i = 0; i< numDims; ++i )
      {
        centroid[i] /= totalVolume;
        positionDifference[i] = centroid[i] - particlePosition[p][i];
      }
      real64 deltaSquared = 0.0;
      for( localIndex i = 0; i < numDims; ++i )
      {
        deltaSquared += positionDifference[i] * positionDifference[i];
      }
      // 0.08 was hand-picked based on testing done in Mathematica. Works in 2D and 3D!
      if( particleSurfaceFlag[p] != 1 && deltaSquared >= 0.08 * 0.08 * neighborRadius * neighborRadius )
      {
        particleSurfaceFlag[p] = 1;
      }
    } );
} );
}

// Before recomputing these
// -heal all damage
// -compute damage gradient without particle surface normals (unless we plan not to overwrite previously exisiting ones)
// -reset F for cpdidomainscaled particles
// -reset cohesive surface flag damage
void SolidMechanicsMPM::computeParticleSurfaceNormalsAndPositions( ParticleManager & particleManager,
                                                                   NodeManager & nodeManager )
{
  // Map particle mass to grid
  GEOS_MARK_FUNCTION;

  int const numDims = m_numDims;
  real64 const explicitSurfaceNormalInfluence = m_explicitSurfaceNormalInfluence;

  // Grid fields
  arrayView2d< real64 > const & gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView3d< real64 > const & gridSurfaceNormal = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridSurfaceNormalString() );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    
    // Map to grid
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< localIndex const > const mappedFields = m_mappedFields[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        real64 particleContributionToGrid;

        localIndex const mappedNode = mappedNodes[pp][g];
        localIndex const fieldIndex = mappedFields[pp][g];
        real64 const shapeFunctionValue = shapeFunctionValues[pp][g];
        real64 shapeFunctionGradientValue [3];
        LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[pp][g] );

       
        particleContributionToGrid = particleMass[p] * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &gridMass[mappedNode][fieldIndex], particleContributionToGrid );

        for( int i = 0; i < numDims; ++i )
        {
          particleContributionToGrid += shapeFunctionGradientValue[i] * particleVolume[p];
          if( particleSurfaceFlag[p] == 2 || particleSurfaceFlag[p] == 3 || particleSurfaceFlag[p] == 4 ) // Update this with enum type for type safety to specifically only implement 2 and 3 for now (those with explicit surface normals)
          {
            // Also maps explicit particle surface normals which will dominate if m_explicitSurfaceNormalInfluence is large
            // If particle surface normal was disabled due to damage or CPDI domain scaling (e.g. zeroed) then the following does not add anything to gridSurfaceNormal
            particleContributionToGrid += explicitSurfaceNormalInfluence * particleSurfaceNormal[p][i] * shapeFunctionValue * particleVolume[p];
          }
          RAJA::atomicAdd( parallelDeviceAtomic{}, &gridSurfaceNormal[mappedNode][fieldIndex][i], particleContributionToGrid );
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop

  switch( m_normalAndPositionMethod )
  {
      case NormalsAndPositionsMethodOption::LR:
        {
          logisticRegression( particleManager,
                              nodeManager );
          break;
        }
      default:
        {
          GEOS_ERROR( "Unknown particle surface normal and position computation method." );
          break;
        }
  }

  mapSurfaceNormalsAndPositionsToParticles( particleManager,
                                            nodeManager );

  // Zero grid mass for now so we don't mass anything up downstream
  int const numNodes = nodeManager.size();
  forAll< parallelDevicePolicy<> >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const g )
  {
    for( integer f = 0; f < m_numVelocityFields; ++f )
    {
      gridMass[g][f] = 0.0;
    }
  } ); // node loop
}

void SolidMechanicsMPM::logisticRegression1D( ParticleManager & particleManager,
                                              NodeManager & nodeManager )
{
  localIndex const numDims = m_numDims;
  localIndex const numContactGroups = m_numContactGroups;
  int const damageFieldPartitioning = m_damageFieldPartitioning;
  int const numNodes = nodeManager.size();

  // Temporarily register fields here that we need
  array2d< real64 > nodeSurfaceMass( numNodes, m_numVelocityFields );
  array3d< real64 > nodeSurfaceCenterOfMass( numNodes, m_numVelocityFields, 3 );

  // Initialize temporary fields
  for( int n = 0; n < numNodes; ++n )
  {
    for( int f = 0; f < m_numVelocityFields; ++f )
    {
      nodeSurfaceMass[n][f] = 0.0;
      LvArray::tensorOps::fill< 3 >( nodeSurfaceCenterOfMass[n][f], 0.0 );
    }
  }

  // Get array views of relevant temporary fields
  arrayView3d< real64 > const nodeSurfaceCenterOfMassView = nodeSurfaceCenterOfMass.toView();

  // Get relevant grid fields
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  // Map relevant values from particles to grid
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< int const > const particleGroup = subRegion.getParticleGroup();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();

    // Get views to mapping arrays
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 const > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];
    
    // ParticleType const particleType = subRegion.getParticleType();
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( localIndex g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        real64 particleContributionToGrid;

        localIndex const mappedNode = mappedNodes[pp][g];
        real64 const shapeFunctionValue = shapeFunctionValues[pp][g];
        real64 shapeFunctionGradientValue [3];
        LvArray::tensorOps::copy< 3 >( shapeFunctionGradientValue, shapeFunctionGradientValues[pp][g] );

        localIndex const fieldIndex = partitionField( numContactGroups,
                                                      damageFieldPartitioning,
                                                      particleGroup[p],
                                                      particleDamageGradient[p],
                                                      particleSurfaceNormal[p],
                                                      gridDamageGradient[mappedNode] );

        particleContributionToGrid = particleMass[p] * ( particleSurfaceFlag[p] > 1 ) * shapeFunctionValue;
        RAJA::atomicAdd( parallelDeviceAtomic{}, &nodeSurfaceMass[mappedNode][fieldIndex], particleContributionToGrid );
        for( int i = 0; i < numDims; ++i )
        {
          particleContributionToGrid = particleMass[p] * ( particleSurfaceFlag[p] > 1 ) * ( particlePosition[p][i] - gridPosition[mappedNode][i] ) * shapeFunctionValue;
          RAJA::atomicAdd( parallelDeviceAtomic{}, &nodeSurfaceCenterOfMassView[mappedNode][fieldIndex][i], particleContributionToGrid );
        }
      }
    } ); // particle loop

    // Increment subregion index
    ++subRegionIndex;
  } ); // subregion loop

  // Normalize node surface center of mass by surface mas
  for( int n = 0; n < numNodes; ++n )
  {
    for( int f = 0; f < m_numVelocityFields; ++f )
    {
      if( nodeSurfaceMass[n][f] > m_smallMass)
      {
        LvArray::tensorOps::scale< 3 >( nodeSurfaceCenterOfMass[n][f], 1/nodeSurfaceMass[n][f] );
      }
    }
  }

  // Get accessors and views
  ParticleManager::ParticleViewAccessor< arrayView1d< localIndex const > > particleGroupAccessor = particleManager.constructArrayViewAccessor< localIndex, 1 >( "particleGroup" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particleSurfaceNormalAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleSurfaceNormal" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particleDamageGradientAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleDamageGradient" );

  ParticleManager::ParticleViewConst< arrayView1d< localIndex const > > particleGroupView =  particleGroupAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePositionView =  particlePositionAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleSurfaceNormalView =  particleSurfaceNormalAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleDamageGradientView =  particleDamageGradientAccessor.toNestedViewConst();

  // Get nodal neighbor list views
  arrayView1d< localIndex const > const numNeighborsAll = m_nodalNeighborList.m_numParticles.toViewConst();
  ArrayOfArraysView< localIndex const > const neighborRegions = m_nodalNeighborList.m_toParticleRegion.toViewConst();
  ArrayOfArraysView< localIndex const > const neighborSubRegions = m_nodalNeighborList.m_toParticleSubRegion.toViewConst();
  ArrayOfArraysView< localIndex const > const neighborIndices = m_nodalNeighborList.m_toParticleIndex.toViewConst();

  // Compute normals and positions on grid use logistic regression along shared normal determined by vector between center of masses from either velocity field
  for( int n = 0; n < numNodes; ++n )
  {
    for( int a = 0; a < m_numVelocityFields-1; ++a )
    {
      if( nodeSurfaceMass[n][a] <= m_smallMass )
      {
        continue;
      }

      for( int b = a+1; b < m_numVelocityFields; ++b )
      {
        if( nodeSurfaceMass[n][b] <= m_smallMass )
        {
          continue;
        }

        real64 xAB[3];
        LvArray::tensorOps::copy< 3 >( xAB, nodeSurfaceCenterOfMass[n][b] );
        LvArray::tensorOps::subtract< 3 >( xAB, nodeSurfaceCenterOfMass[n][a] );

        real64 nAB[3];
        LvArray::tensorOps::copy< 3 >( nAB, xAB );
        LvArray::tensorOps::normalize< 3 >( nAB );

        // Iterative solve of logistic regression in 1D to compute the middle of the interface along the shared normal
        // Requires the node-particle neighborlist

        // Compute the project coordinates of each particle
        int numNeighbors = numNeighborsAll[n];
        for( int p = 0; p < numNeighbors; ++p )
        {
          localIndex const ri = neighborRegions[n][p];
          localIndex const sri = neighborSubRegions[n][p];
          localIndex const pi = neighborIndices[n][p];

          // Check if particle maps to either velocity field currently under consideration
          localIndex fieldIndex = partitionField( numContactGroups,
                                                  damageFieldPartitioning,
                                                  particleGroupView[ri][sri][pi],
                                                  particleDamageGradientView[ri][sri][pi],
                                                  particleSurfaceNormalView[ri][sri][pi],
                                                  gridDamageGradient[n] ); 
          if( fieldIndex != a && fieldIndex != b )
          {
            continue;
          }

          real64 relPos[3];
          LvArray::tensorOps::copy< 3 >( relPos, particlePositionView[ri][sri][pi] );
          LvArray::tensorOps::subtract< 3 >( relPos, gridPosition[n] );
          real64 distAlongNorm = LvArray::tensorOps::AiBi< 3 >( relPos, nAB );
          real64 y = fieldIndex == a;
          
          GEOS_UNUSED_VAR(distAlongNorm);
          GEOS_UNUSED_VAR(y);
        }
      }
    }
  }

  // Map normals and positions back to particles

}

void SolidMechanicsMPM::computeSphF( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const planeStrain = m_planeStrain;
  real64 const neighborRadius = m_neighborRadius;

  // Get accessors for volume, position, reference position
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particleReferencePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleReferencePosition" );

  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > particleVolumeView = particleVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePositionView = particlePositionAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particleReferencePositionView = particleReferencePositionAccessor.toNestedViewConst();

  // Perform neighbor operations
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    // Get particle position and sphF
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView3d< real64 > const particleSPHF = subRegion.getField< fields::mpm::particleSPHF >();

    // Loop over neighbors
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Get number of neighbors and accessor indices
      localIndex numNeighbors = numNeighborsAll[p];
      arraySlice1d< localIndex const > const regionIndices = neighborRegions[pp];
      arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[pp];
      arraySlice1d< localIndex const > const particleIndices = neighborIndices[pp];

      // Declare and size neighbor data arrays - TODO: switch to std::array? But then we'd need to template computeKernelFieldGradient
      // stdVector< real64 > neighborVolumes( numNeighbors );
      // stdVector< stdVector< real64 > > neighborPositions;
      // neighborPositions.resize( numNeighbors, stdVector< real64 >( 3 ) );
      // stdVector< stdVector< real64 > > neighborDisplacements;
      // neighborDisplacements.resize( numNeighbors, stdVector< real64 >( 3 ) );
      //
      // // Populate neighbor data arrays
      // for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
      // {
      //   localIndex regionIndex = regionIndices[neighborIndex];
      //   localIndex subRegionIndex = subRegionIndices[neighborIndex];
      //   localIndex particleIndex = particleIndices[neighborIndex];
      //   real64 r0Squared = 0.0;
      //   for( int i=0; i<m_numDims; ++i )
      //   {
      //     real64 thisComponent = particleReferencePosition[p][i];
      //     real64 neighborComponent = particleReferencePositionView[regionIndex][subRegionIndex][particleIndex][i];
      //     r0Squared += ( thisComponent - neighborComponent ) * ( thisComponent - neighborComponent );
      //   }
      //   if( r0Squared <= m_neighborRadius * m_neighborRadius ) // We only consider neighbor particles that would have been neighbors in
      //                                                           // the reference configuration
      //   {
      //     neighborVolumes[neighborIndex] = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
      //   }
      //   else
      //   {
      //     neighborVolumes[neighborIndex] = 0.0;
      //   }
      //   for( int i=0; i<3; ++i )
      //   {
      //     neighborPositions[neighborIndex][i] = particlePositionView[regionIndex][subRegionIndex][particleIndex][i];
      //     neighborDisplacements[neighborIndex][i] = particlePositionView[regionIndex][subRegionIndex][particleIndex][i] -
      //                                               particleReferencePositionView[regionIndex][subRegionIndex][particleIndex][i];
      //   }
      // }

      // // Call kernel field gradient function
      // computeKernelVectorGradient( particlePosition[p],        // input
      //                              numNeighbors,
      //                              neighborPositions,          // input
      //                              neighborVolumes,            // input
      //                              neighborDisplacements,      // input
      //                              particleSPHF[p] );          // OUTPUT: Spatial displacement gradient

      computeKernelVectorGradientDevice( particlePosition[p],        // input
                                         neighborRadius,
                                         planeStrain,
                                         numNeighbors,
                                         regionIndices,
                                         subRegionIndices,
                                         particleIndices,
                                         particleVolumeView,
                                         particlePositionView,
                                         particleReferencePositionView,
                                         particleSPHF[p] );          // OUTPUT: Spatial displacement gradient

      // Convert spatial displacement gradient to deformation gradient
      real64 temp[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
      LvArray::tensorOps::scaledAdd< 3, 3 >( temp, particleSPHF[p], -1.0 );
      LvArray::tensorOps::invert< 3 >( temp );
      LvArray::tensorOps::copy< 3, 3 >( particleSPHF[p], temp );
    } );
  } );
}

// void SolidMechanicsMPM::directionalOverlapCorrection( real64 dt, ParticleManager & particleManager )
// {
//   // Get the SPH version of the deformation gradient
//   computeSphF( particleManager );

//   // If we're at a surface, induce corrective deformation normal to the surface
//   particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
//   {
//     arrayView1d< real64 > const particleOverlap = subRegion.getField< fields::mpm::particleOverlap >();
//     particleOverlap.zero();
//     arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();
//     arrayView3d< real64 > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
//     particleVelocityGradient.zero();
//     arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
//     arrayView3d< real64 const > const particleSPHF = subRegion.getField< fields::mpm::particleSPHF >();

//     SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
//     forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST ( localIndex const pp )
//     {
//       localIndex const p = activeParticleIndices[pp];
//       if( LvArray::tensorOps::l2NormSquared< 3 >( particleDamageGradient[p] ) > 0.0 )
//       {
//         real64 surfaceNormal[3];
//         LvArray::tensorOps::copy< 3 >( surfaceNormal, particleDamageGradient[p] );
//         if( m_planeStrain == 1 )
//         {
//           surfaceNormal[2] = 0.0; // Just in case
//         }
//         LvArray::tensorOps::normalize< 3 >( surfaceNormal );

//         // Polar decompositions
//         real64 Rp[3][3], Rsph[3][3], Vp[3][3], Vsph[3][3];
//         LvArray::tensorOps::polarDecomposition< 3 >( Rp, particleDeformationGradient[p] );
//         LvArray::tensorOps::polarDecomposition< 3 >( Rsph, particleSPHF[p] );
//         LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >(Vp, particleDeformationGradient[p], Rp );
//         LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 3 >(Vsph, particleSPHF[p], Rsph );

//         // Compute directional stretches
//         real64 stretch[3], normalStretch;
//         real64 stretchSph[3], normalStretchSph;
//         LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >(stretch, Vp, surfaceNormal );
//         normalStretch = LvArray::tensorOps::AiBi< 3 >( stretch, surfaceNormal );
//         LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >(stretchSph, Vsph, surfaceNormal );
//         normalStretchSph = LvArray::tensorOps::AiBi< 3 >( stretchSph, surfaceNormal );

//         // Detect overlap and correct. TODO: Implement a ramp?
//         real64 overlap = normalStretch / normalStretchSph;
//         particleOverlap[p] = overlap;
//         if( overlap > 1.2 && normalStretchSph < 1 )
//         {
//           // Get orthonormal basis to transform V into
//           real64 s1[3], s2[3];
//           computeOrthonormalBasis( surfaceNormal, s1, s2 );
//           real64 Q[3][3], Qtranspose[3][3], targetF[3][3], targetFDot[3][3], FInverse[3][3];
//           for( int i=0; i<3; ++i )
//           {
//             Q[i][0] = surfaceNormal[i];
//             Q[i][1] = s1[i];
//             Q[i][2] = s2[i];
//             Qtranspose[0][i] = surfaceNormal[i];
//             Qtranspose[1][i] = s1[i];
//             Qtranspose[2][i] = s2[i];
//           }

//           // Perform transformation, substitute in the directional stretch predicted by Fsph, transform back, re-assemble F
//           real64 VpPrime[3][3], temp[3][3];
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( temp, Qtranspose, Vp );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( VpPrime, temp, Q );
//           VpPrime[0][0] = normalStretchSph;
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( temp, Q, VpPrime );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( Vp, temp, Qtranspose );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( targetF, Vp, Rp );

//           // Deduce the L necessary to achieve the new F
//           LvArray::tensorOps::copy< 3, 3 >( targetFDot, targetF );
//           LvArray::tensorOps::scaledAdd< 3, 3 >( targetFDot, particleDeformationGradient[p], -1.0 );
//           LvArray::tensorOps::scale< 3, 3 >( targetFDot, 1.0 / dt );
//           LvArray::tensorOps::invert< 3 >( FInverse, particleDeformationGradient[p] );
//           LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleVelocityGradient[p], targetFDot, FInverse );
//         }
//       }
//     } );
//   } );

//   // TODO: If we're not at a surface, volumetric overlap correction should be performed after the F update
// }

int SolidMechanicsMPM::evaluateSeparabilityCriterion( int const & planeStrain,
                                                      int const & numContactGroups,
                                                      int const & treatFullyDamagedAsSingleField,
                                                      real64 const & separabilityMinDamage,
                                                      real64 const & thinFeatureDFGThreshold,
                                                      real64 const & neighborRadius,
                                                      localIndex const & A,
                                                      localIndex const & B,
                                                      real64 const & damageA,
                                                      real64 const & damageB,
                                                      real64 const & maxDamageA,
                                                      real64 const & maxDamageB,
                                                      arraySlice1d< real64 const > const dmgGrad,
                                                      arraySlice1d< real64 const > const xA,
                                                      arraySlice1d< real64 const > const xB )
// m_treatFullyDamagedAsSingleField makes fields inseparable if damageA = damageB = 1, so we aren't putting
// arbitrary separation planes (and potential surfaces for accumulated overlap that needs corrections)
// between fully damaged materials. There is a potential issue that approaching damaged bodies
// would get contact gaps locked in.
{
  int separable = 0;
  // At least one field is fully damaged and both fields have the minimum separable level of damage.
  // The "a%b" is the "mod(a,b)" command, and indicates whether materials are from same contact group.
  if( ( ( maxDamageA >= 0.9999 || maxDamageB >= 0.9999 ) &&
        ( damageA >= separabilityMinDamage && damageB >= separabilityMinDamage ) ) ||
      ( A % numContactGroups != B % numContactGroups ) )
  {
    real64 xi = 1e-3;
    if( planeStrain == 1 )
    {
      xi /= m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1];
    }
    else
    {
      xi /= m_hEl[0] * m_hEl[0] + m_hEl[1] * m_hEl[1] + m_hEl[2] * m_hEl[2];
    }

    // We artificially map damage for surface particles with explicit normals and positions as one to gridMaxDamage and gridDamage
    // If this flag is set, the contact fields are treated as inseparable which is the wrong behavior
    if( treatFullyDamagedAsSingleField == 1 && LvArray::tensorOps::l2NormSquared< 3 >( dmgGrad ) < xi )
    {
      // separable = ( fmin( damageA, damageB ) < 0.9999 ) ? 1 : 0; // Might not need this check with additional comparision of damage gradient in if above
      separable = 0;
    }
    else
    {
      separable = 1;
    }
  }

  // For a thin (relative to grid spacing) strip of material with surface flags or damaged surfaces, the material at opposing surfaces will have opposing
	// damage-field gradients, and this will create a spurious slip surfaces within the material.
	// We can detect this case by seeing if the distance between the position of the two fields is small relative to the grid cell spacing.  This will
	// only happen in the case of a thin strip of material when using DFG, so in this case we set separable = 0
	// The threshold spacing is normalized by the neighbor radius.
	if ( thinFeatureDFGThreshold < 10.0 ) // Threshold defaults to 1.e99 so we only do this calculation if a lower value is set.
	{                                       // After testing we may default to a lower value.
		real64 dx[3] = {};
    LvArray::tensorOps::copy< 3 >( dx, xA );
    LvArray::tensorOps::subtract< 3 >( dx, xB );

    real64 productOfSquares = ( dx[0] * dx[0] ) * ( dx[1] * dx[1] ) * ( dx[2] * dx[2] );

    if ( productOfSquares < ( thinFeatureDFGThreshold * neighborRadius ) * ( thinFeatureDFGThreshold * neighborRadius ) )
		{
			separable = 0;
		}
	}

  return separable;
}

// All master particles should have centers inside the domain if particleCenters are corrected correctly during repartitioning
// CPDI: Edge case, corner of large or long particle beyond ghost cells but center is still inside domain?
void SolidMechanicsMPM::flagOutOfRangeParticles( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Identify global domain bounds
    real64 globalMin[3]; // including buffer cells
    real64 globalMax[3]; // including buffer cells
    for( int i=0; i<3; ++i )
    {
      globalMin[i] = m_xGlobalMin[i] - m_hEl[i];
      globalMax[i] = m_xGlobalMax[i] + m_hEl[i];
    }

    // Get particle fields
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< int > const particleDeleteFlag = subRegion.getField< fields::mpm::particleDeleteFlag >();
    ParticleType particleType = subRegion.getParticleType();

    // Define tolerance
    real64 tolerance[3];
    for( int i=0; i<3; ++i )
    {
      tolerance[i] = m_hEl[i] * DBL_EPSILON;
    }

    switch( particleType )
    {
      case ParticleType::SinglePoint:
        {
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
          {
            localIndex const p = activeParticleIndices[pp];
            for( int i=0; i<3; ++i )
            {
              if( particlePosition[p][i] < globalMin[i] + tolerance[i] || globalMax[i] - tolerance[i] < particlePosition[p][i] )
              {
                particleDeleteFlag[p] = 1;
                break; // TODO: if this doesn't work, just modify "i"
              }
            }
          } );
          break;
        }
      case ParticleType::CPDI:
        {
          arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
          int const signs[8][3] = { { 1, 1, 1},
            { 1, 1, -1},
            { 1, -1, 1},
            { 1, -1, -1},
            {-1, 1, 1},
            {-1, 1, -1},
            {-1, -1, 1},
            {-1, -1, -1} };
          forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
          {
            localIndex const p = activeParticleIndices[pp];
            for( int cornerIndex=0; cornerIndex<8; cornerIndex++ )
            {
              for( int i=0; i<3; ++i )
              {
                real64 cornerPositionComponent = particlePosition[p][i] + signs[cornerIndex][0] * particleRVectors[p][0][i] + signs[cornerIndex][1] * particleRVectors[p][1][i] + signs[cornerIndex][2] *
                                                particleRVectors[p][2][i];
                if( cornerPositionComponent < globalMin[i] + tolerance[i] || globalMax[i] - tolerance[i] < cornerPositionComponent )
                {
                  particleDeleteFlag[p] = 1;
                  break;
                }
              }
              if( particleDeleteFlag[p] == 1 )
              {
                break;
              }
            }
          } );
          break;
        }
      default:
        {
          GEOS_ERROR( "Particle type \"" << particleType << "\" is not yet supported." );
          break;
        }
    }
  } );
}

void SolidMechanicsMPM::computeRVectors( ParticleManager & particleManager )
{
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    if( subRegion.hasRVectors() )
    {
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
      arrayView3d< real64 const > const particleReferenceRVectors = subRegion.getField< fields::mpm::particleReferenceRVectors >();
      arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        for( int i=0; i<3; ++i )
        {
          for( int j=0; j<3; ++j )
          {
            particleRVectors[p][i][j] = particleReferenceRVectors[p][i][0] * particleDeformationGradient[p][j][0] +
                                        particleReferenceRVectors[p][i][1] * particleDeformationGradient[p][j][1] +
                                        particleReferenceRVectors[p][i][2] * particleDeformationGradient[p][j][2];
          }
        }
      } );
    }
  } );
}

void SolidMechanicsMPM::cpdiDomainScaling( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const planeStrain = m_planeStrain;
  int const disableSurfaceNormalsAndPositionsOnCPDIScaling = m_disableSurfaceNormalsAndPositionsOnCPDIScaling;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    if( subRegion.getParticleType() == ParticleType::CPDI )
    {
      real64 const lCrit = planeStrain == 1 ? 0.49999 * fmin( hEl[0], hEl[1] ) : 0.49999 * fmin( hEl[0], fmin( hEl[1], hEl[2] ) );
      arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
      arrayView2d< real64 > const particleReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
      arrayView2d< real64 > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
      arrayView2d< real64 > const particleReferenceSurfacePosition = subRegion.getField< fields::mpm::particleReferenceSurfacePosition >();
      arrayView2d< real64 > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
      arrayView1d< int > const particleDomainScaledFlag = subRegion.getField< fields::mpm::particleDomainScaledFlag >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
      {
        // Reset cpdi domain scaled flag
        particleDomainScaledFlag[p] = 0;

        arraySlice1d< real64 > const r1 = particleRVectors[p][0];
        arraySlice1d< real64 > const r2 = particleRVectors[p][1];
        arraySlice1d< real64 > const r3 = particleRVectors[p][2];

        bool scale = false;
        if( planeStrain == 1 ) // 2D cpdi domain scaling
        {
          // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
          real64 l[2][3];
          for( int i=0; i<3; ++i )
          {
            l[0][i] = r1[i] + r2[i]; // la
            l[1][i] = r1[i] - r2[i]; // lb
          }

          // scale l-vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
          for( int i = 0; i < 2; ++i )
          {
            real64 lLength = LvArray::math::sqrt( l[i][0] * l[i][0] + l[i][1] * l[i][1] + l[i][2] * l[i][2] );
            if( lLength > lCrit )
            {
              l[i][0] *= lCrit / lLength;
              l[i][1] *= lCrit / lLength;
              l[i][2] *= lCrit / lLength;
              scale = true;
            }
          }

          // reconstruct r-vectors.  eq. 11 in the CPDI domain scaling paper.
          if( scale )
          {
            for( int i=0; i<3; ++i )
            {
              r1[i] = 0.5 * (l[0][i] + l[1][i]);
              r2[i] = 0.5 * (l[0][i] - l[1][i]);
            }
          }
        }
        else // 3D cpdi domain scaling
        {
          // Initialize l-vectors.  Eq. 8a-d in the CPDI domain scaling paper.
          real64 l[4][3];
          for( int i=0; i<3; ++i )
          {
            l[0][i] = r1[i] + r2[i] + r3[i]; // la
            l[1][i] = r1[i] - r2[i] + r3[i]; // lb
            l[2][i] = r2[i] - r1[i] + r3[i]; // lc
            l[3][i] = r3[i] - r1[i] - r2[i]; // ld
          }

          // scale l vectors if needed.  Eq. 9 in the CPDI domain scaling paper.
          for( int i = 0; i < 4; ++i )
          {
            real64 lLength = LvArray::math::sqrt( l[i][0] * l[i][0] + l[i][1] * l[i][1] + l[i][2] * l[i][2] );
            if( lLength > lCrit )
            {
              l[i][0] *= lCrit / lLength;
              l[i][1] *= lCrit / lLength;
              l[i][2] *= lCrit / lLength;
              scale = true;
            }
          }

          // reconstruct r vectors.  eq. 11 in the CPDI domain scaling paper.
          if( scale )
          {
            for( int i=0; i<3; ++i )
            {
              r1[i] = 0.25 * ( l[0][i] + l[1][i] - l[2][i] - l[3][i] );
              r2[i] = 0.25 * ( l[0][i] - l[1][i] + l[2][i] - l[3][i] );
              r3[i] = 0.25 * ( l[0][i] + l[1][i] + l[2][i] + l[3][i] );
            }
          }
        }

        if( scale )
        {
          particleDomainScaledFlag[p] = 1;
        }

        // Turn off particle explicit surface normals and posititons
        if( disableSurfaceNormalsAndPositionsOnCPDIScaling == 1 && particleDomainScaledFlag[p] == 1 )
        {
          LvArray::tensorOps::fill< 3 >( particleSurfaceNormal[p], 0.0 );
          LvArray::tensorOps::fill< 3 >( particleSurfacePosition[p], 0.0 );

          LvArray::tensorOps::fill<3>( particleReferenceSurfaceNormal[p], 0.0 );
          LvArray::tensorOps::fill<3>( particleReferenceSurfacePosition[p], 0.0 );
        }
      } );
    }
  } );
}

// Should only be an option for CPDI, CPTI, and CPDI2 particles, right?
void SolidMechanicsMPM::subdivideParticles( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const numDims = m_numDims;

  // Determine critical length for subdividing particles from problem dimensions
  real64 const lCritSqr = LvArray::math::pow( m_planeStrain == 1 ? 0.49999 * fmin( m_hEl[0], m_hEl[1] ) : 0.49999 * fmin( m_hEl[0], fmin( m_hEl[1], m_hEl[2] ) ), 2 );

  // Find max particle ID for assigning new particle IDs
  // Determine number of new particles per region
  RAJA::ReduceMax< parallelDeviceReduce, globalIndex > maxParticleIDRank( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, int > numDivisibleParticles( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, int > numNewParticles( 0 );
  array1d< int > numDivisibleParticlesPerSubRegion( m_numberOfSubRegions );
  array1d< int > numNewParticlesPerSubRegion( m_numberOfSubRegions );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
    arrayView1d< int > const particleSubdivideFlag = subRegion.getField< fields::mpm::particleSubdivideFlag >();

    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionNumDivisibleParticles( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionNumNewParticles( 0 );

    // Loop over all particles, because there should be no ghost particles yet
    forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
    {
      maxParticleIDRank.max( particleID[p] );

      // If particle rvector is beyond characteristic grid cell direction, divide along that direction
      // This prevents us from unnecessarily adding more particles than we need to
      int numberOfDivisions = 0;
      for(int i = 0; i < numDims; ++i)
      {
        if( LvArray::tensorOps::l2NormSquared< 3 >( particleRVectors[p][i] ) > lCritSqr )
        {
          ++numberOfDivisions;
        }
      }

      if( numberOfDivisions > 0 )
      {
        particleSubdivideFlag[p] = 1;
        subRegionNumDivisibleParticles += 1;
        subRegionNumNewParticles += LvArray::math::pow(2, numberOfDivisions) - 1;
        numDivisibleParticles += 1;
        numNewParticles += LvArray::math::pow(2, numberOfDivisions) - 1;
      }
    } );

    numDivisibleParticlesPerSubRegion[subRegionIndex] = subRegionNumDivisibleParticles;
    numNewParticlesPerSubRegion[subRegionIndex] = subRegionNumNewParticles;

    ++subRegionIndex;
  } );

  // Reduce global ID
  globalIndex maxParticleID = MpiWrapper::max( maxParticleIDRank.get() );

  // Gather array for number of new particles per rank for assigning global IDs
  array1d< int > numNewParticlesPerRank( MpiWrapper::commSize() );
  MpiWrapper::allGather( numNewParticles.get(), numNewParticlesPerRank, MPI_COMM_GEOS );

  int rank = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  globalIndex currGlobalIndex = maxParticleID+1;
  int totalNewParticles = 0;
  for( int i = 0; i < numNewParticlesPerRank.size(); ++i )
  {
    if( i < rank-1)
    {
      currGlobalIndex += numNewParticlesPerRank[i];
    }
    totalNewParticles += numNewParticlesPerRank[i];
  }

  // Subdivide particles
  subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    localIndex oldSubRegionSize = subRegion.size();
    localIndex newSubRegionSize = oldSubRegionSize + numNewParticlesPerSubRegion[subRegionIndex];
    subRegion.resize( newSubRegionSize );

    // It looks like there was an issue with getting a copy of the arrayview size before resizing subregion.
    // Calling after as reference fixed it, but it is not clear if calling an arrayview reference before subregion resize also works, NEED TO TEST!!!!!!!
    arrayView1d< globalIndex > const & particleID = subRegion.getParticleID();
    arrayView1d< real64 > const & particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const & particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 > const & particleReferenceVolume = subRegion.getField< fields::mpm::particleReferenceVolume >();
    arrayView2d< real64 > const & particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const & particleReferencePosition = subRegion.getField< fields::mpm::particleReferencePosition >();
    arrayView3d< real64 > const & particleRVectors = subRegion.getParticleRVectors();
    arrayView3d< real64 > const & particleReferenceRVectors = subRegion.getField< fields::mpm::particleReferenceRVectors >();
    arrayView1d< int > const & particleSubdivideFlag = subRegion.getField< fields::mpm::particleSubdivideFlag >();
    arrayView1d< int > const & particleCopyFlag = subRegion.getField< fields::mpm::particleCopyFlag >();

    localIndex subRegionNewParticleIndex = oldSubRegionSize;

    // Only iterate over old particles, assumes on resize they are all at the front of the arrays
    // TODO convert this to RAJA forAll for device parallelization
    // forAll< serialPolicy >( oldSubRegionSize, [=, &currGlobalIndex, &subRegionNewParticleIndex] GEOS_HOST_DEVICE ( localIndex const p )
    for( localIndex p = 0; p < oldSubRegionSize; ++p)
    {
      if( particleSubdivideFlag[p] == 1 )
      {
        // Determine which directions to subdivide along
        int numDivisions = 0;
        array1d< array1d< int > > rVectorDivisions( numDims );
        array1d< int > subdivideDirections( numDims );
        LvArray::tensorOps::fill< 3 >( subdivideDirections, 0 );
        for( int d = 0; d < numDims; ++d)
        {
          // Check if any of the RVectors have lengths beyond the critical (as determined by the grid cell)
          if( LvArray::tensorOps::l2NormSquared< 3 >( particleRVectors[p][d] ) > lCritSqr )
          {
            subdivideDirections[d] = 1;
            rVectorDivisions[d].resize(2);
            rVectorDivisions[d][0] = 1;
            rVectorDivisions[d][1] = -1;
            ++numDivisions;
          }
          else
          {
            rVectorDivisions[d].resize(1);
            rVectorDivisions[d][0] = 0;
          }
        }

        auto rVectorOffsetCombinations = generateCombinations( rVectorDivisions );

        //Subdivide particle and copy particle field data
        int subdivideFactor = LvArray::math::pow(2, numDivisions);
        real64 newMass = particleMass[p] / subdivideFactor;
        real64 newVolume = particleVolume[p] / subdivideFactor;
        real64 newReferenceVolume = particleReferenceVolume[p] / subdivideFactor;
        for(int np = 1; np < subdivideFactor; ++np )
        {
          // Update particle mass, volume, initial volume, centers, reference positions, initial R vectors and Rvectors
          particleID[subRegionNewParticleIndex] = currGlobalIndex++;
          particleMass[subRegionNewParticleIndex] = newMass;
          particleVolume[subRegionNewParticleIndex] = newVolume;
          particleReferenceVolume[subRegionNewParticleIndex] = newReferenceVolume;

          LvArray::tensorOps::copy< 3, 3 >( particleReferenceRVectors[subRegionNewParticleIndex], particleReferenceRVectors[p] );
          LvArray::tensorOps::copy< 3, 3 >( particleRVectors[subRegionNewParticleIndex], particleRVectors[p] );

          for( int d = 0; d < numDims; ++d )
          {
            if( subdivideDirections[d] == 1 )
            {
              LvArray::tensorOps::scale< 3 >( particleReferenceRVectors[subRegionNewParticleIndex][d], 0.5 );
              LvArray::tensorOps::scale< 3 >( particleRVectors[subRegionNewParticleIndex][d], 0.5 );
            }
          }

          LvArray::tensorOps::copy< 3 >( particlePosition[subRegionNewParticleIndex], particlePosition[p] );
          LvArray::tensorOps::copy< 3 >( particleReferencePosition[subRegionNewParticleIndex], particleReferencePosition[p] );
          for(int di =0; di < numDims; ++di)
          {
            for(int dj =0; dj < numDims; ++dj)
            {
              particlePosition[subRegionNewParticleIndex][dj] += rVectorOffsetCombinations[np][di] * particleRVectors[subRegionNewParticleIndex][di][dj];
              particleReferencePosition[subRegionNewParticleIndex][dj] += rVectorOffsetCombinations[np][di] * particleReferenceRVectors[subRegionNewParticleIndex][di][dj];
            }
          }

          // Probably need another flag since casting from localIndex (unsigned) to int could potentially overflow
          // but for now we need -1 to screen particles that should not copy
          particleCopyFlag[subRegionNewParticleIndex] = static_cast< int >( p );

          ++subRegionNewParticleIndex;
        }

        // Modifying original particle (globalID does not need updating)
        particleMass[p] = newMass;
        particleVolume[p] = newVolume;
        particleReferenceVolume[p] = newReferenceVolume;

        for( int d = 0; d < numDims; ++d )
        {
          if( subdivideDirections[d] == 1 )
          {
            LvArray::tensorOps::scale< 3 >( particleReferenceRVectors[p][d], 0.5 );
            LvArray::tensorOps::scale< 3 >( particleRVectors[p][d], 0.5 );
          }
        }

        for(int di =0; di < numDims; ++di)
        {
          for(int dj =0; dj < numDims; ++dj)
          {
            particlePosition[p][dj] += rVectorOffsetCombinations[0][di] * particleRVectors[p][di][dj];
            particleReferencePosition[p][dj] += rVectorOffsetCombinations[0][di] * particleReferenceRVectors[p][di][dj];
          }
        }

        // Turn off subdivide flag for particle before copying to new particles
        particleSubdivideFlag[p] = 0;
      }
    }
    // );

    // Copy all other fields that do not need modification
    std::set< std::string > ignoreFieldCopy( { "particleID",
                                               "particleMass",
                                               "particleVolume",
                                               "particleReferenceVolume",
                                               "particleCenter",
                                               "particleReferencePosition",
                                               "particleRVectors",
                                               "particleReferenceRVectors",
                                               "particleCopyFlag" } );

    subRegion.forWrappers( [&]( WrapperBase & fieldWrapper )
    {
      string const fieldName = fieldWrapper.getName();

      // Filter out only particle fields for copy by prefix
      if( fieldName.substr(0, 8) != "particle" || ignoreFieldCopy.count( fieldName ) > 0)
      {
        return;
      }

      types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
      {
        using ArrayType = camp::first< decltype( tupleOfTypes ) >;
        using T = typename ArrayType::ValueType;

        auto sourceArray = Wrapper< ArrayType >::cast( fieldWrapper ).reference().toView();

        forAll< serialPolicy >( sourceArray.size( 0 ), [&]( localIndex const pp )
        {
          if( particleCopyFlag[pp] >= 0 )
          {
            // For scalar particle fields need to do assignment manually
            if constexpr( ArrayType::NDIM == 1 )
            {
              // If wrapper name is among those that should be overriden by new sub region skip (e.g. mass, density, etc.)
              // We currently only overwrite scalar quantities, but may need to adjust if we overwrite nonscalar quantities
              sourceArray[pp] = sourceArray[static_cast< localIndex >( particleCopyFlag[pp] )];
            }
            else
            {
              auto destinationSlice = sourceArray[pp];
              auto sourceSlice = sourceArray[static_cast< localIndex >( particleCopyFlag[pp] )];
              LvArray::forValuesInSliceWithIndices( destinationSlice, [slice=sourceSlice] ( T & val, auto const ... indices )
              {
                val = slice( indices ... );
              } );
            }
          }
        } );

        }, fieldWrapper );
    } );

    // Copy constitutive model fields (e.g. stresses)
    string const & modelName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, modelName );

    constitutiveModel.forWrappers( [&]( WrapperBase & fieldWrapper )
    {
      // Do not copy default and reference scalar values from constitutive model
      if(fieldWrapper.numArrayDims() == 0)
      {
        return;
      }

      types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
      {
        using ArrayType = camp::first< decltype( tupleOfTypes ) >;
        using T = typename ArrayType::ValueType;

        auto sourceArray = Wrapper< ArrayType >::cast( fieldWrapper ).reference().toView();
        forAll< serialPolicy >( sourceArray.size( 0 ), [&]( localIndex const pp )
        {
          if( particleCopyFlag[pp] >= 0 )
          {
            if constexpr( ArrayType::NDIM == 1 )
            {
              sourceArray[pp] = sourceArray[static_cast< localIndex >( particleCopyFlag[pp] )];
            }
            else
            {
              auto destinationSlice = sourceArray[pp];
              auto sourceSlice = sourceArray[static_cast< localIndex >( particleCopyFlag[pp] )];
              LvArray::forValuesInSliceWithIndices( destinationSlice, [slice=sourceSlice] ( T & val, auto const ... indices )
              {
                val = slice( indices ... );
              } );
            }
          }
        } );

        }, fieldWrapper );
    } );

    //Erase copy flags after completing particle subdivision
    forAll< serialPolicy >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
    {
      particleCopyFlag[p] = -1;
    } );

    // Rebuild active particle indices to include new particles
    subRegion.setActiveParticleIndices();
    ++subRegionIndex;
  } );

  GEOS_LOG_RANK_IF( totalNewParticles > 0, "Generated " << totalNewParticles << " particles from subdividing overly deformed particles!"  );
}

void SolidMechanicsMPM::resizeMappingArrays( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  m_mappedNodes.resize( m_numberOfSubRegions );
  m_mappedFields.resize( m_numberOfSubRegions );
  m_shapeFunctionValues.resize( m_numberOfSubRegions );
  m_shapeFunctionGradientValues.resize( m_numberOfSubRegions );

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    int const numberOfActiveParticles = subRegion.activeParticleIndices().size();
    m_mappedNodes[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle() );
    m_mappedFields[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle() );
    m_shapeFunctionValues[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle() );
    m_shapeFunctionGradientValues[subRegionIndex].resize( numberOfActiveParticles, 8 * subRegion.numberOfVerticesPerParticle(), 3 );
    ++subRegionIndex;
  } );
}


void SolidMechanicsMPM::correctGhostParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                                             SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< int const > const periodic = partition.getPeriodic();
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMinNoGhost );
  real64 xLocalMax[3];
  LvArray::tensorOps::copy< 3 >( xLocalMax, m_xLocalMaxNoGhost );

  real64 const neighborRadius = m_neighborRadius;

  real64 xLocalMinWithGhost[3];
  LvArray::tensorOps::copy< 3 >( xLocalMinWithGhost, m_xLocalMin );
  real64 xLocalMaxWithGhost[3];
  LvArray::tensorOps::copy< 3 >( xLocalMaxWithGhost, m_xLocalMax );

  real64 xMin[3];
  LvArray::tensorOps::copy< 3 >( xMin, m_xGlobalMin );
  real64 xMax[3];
  LvArray::tensorOps::copy< 3 >( xMax, m_xGlobalMax );
  real64 xExtent[3];
  LvArray::tensorOps::copy< 3 >( xExtent, m_domainExtent );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    
    SortedArrayView< localIndex const > const inactiveParticleIndices = subRegion.inactiveParticleIndices();
    forAll< serialPolicy >( inactiveParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = inactiveParticleIndices[pp];

      for( int i=0; i < 3; ++i)
      {
        if(periodic[i])
        {
          real64 sign = 1.0;
          // Determine if partition is located on negative or positive side of domain boundary
          if( LvArray::math::abs( xLocalMin[i] - xMin[i] ) < LvArray::math::abs( xLocalMax[i] - xMax[i] ) )
          {
            sign *= -1.0;
          }

          // particlePosition[p][i] = Mod(particlePosition[p][i]-xGlobalMin[i], xGlobalMax[i]-xGlobalMin[i])+xGlobalMin[i];
          // particlePosition[p][i] = Mod(particlePosition[p][i] - xMin[i] - sign * hEl[i] * LvArray::math::ceil(m_neighborRadius/hEl[i]), xExtent[i]) + xMin[i] + sign * hEl[i] * LvArray::math::ceil(m_neighborRadius/hEl[i]);
          particlePosition[p][i] = Mod(particlePosition[p][i] - xMin[i] - sign * neighborRadius, xExtent[i]) + xMin[i] + sign * neighborRadius;
        }
      }
    });
  });
}


void SolidMechanicsMPM::correctParticleCentersAcrossPeriodicBoundaries( ParticleManager & particleManager,
                                                                        SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  // We make assumptions that the grid is uniformly laid out. This is dangerous and should be fixed in case we decide to have nonuniform grids for computational efficiency

  arrayView1d< int const > const periodic = partition.getPeriodic();
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xMin[3];
  LvArray::tensorOps::copy< 3 >( xMin, m_xGlobalMin );
  // real64 xMax[3];
  // LvArray::tensorOps::copy< 3 >( xMax, m_xGlobalMax);
  real64 xExtent[3];
  LvArray::tensorOps::copy< 3 >( xExtent, m_domainExtent );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    arrayView2d< real64 > particlePosition = subRegion.getParticleCenter();

    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      //Check over every dimension
      for( localIndex i=0; i < 3; ++i )
      {
        if(periodic[i])
        {
          // particlePosition[p][i] = Mod( particlePosition[p][i] - xMin[i], xMax[i] - xMin[i] ) + sxMin[i];
          particlePosition[p][i] = Mod( particlePosition[p][i] - xMin[i], xExtent[i]) + xMin[i];
        }
      }
    });
  });
}

real64 SolidMechanicsMPM::Mod(real64 num, real64 denom)
{
  if(isZero(denom))
  {
    return num;
  }
  return num - denom * LvArray::math::floor(num/denom);
}

int SolidMechanicsMPM::combinations( int n, int k )
{
  return factorial( n )/( factorial( k ) * factorial( n - k ) );
}

int SolidMechanicsMPM::factorial( int n )
{
  int val = 1;
  for(int i = 1; i <= n; ++i){
    val *= i;
  }
  return val;
}

void SolidMechanicsMPM::populateMappingArrays( ParticleManager & particleManager,
                                               NodeManager & nodeManager ) //,
                                              //  SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );
  real64 xLocalMin[3];
  LvArray::tensorOps::copy< 3 >( xLocalMin, m_xLocalMin );
  real64 xLocalMax[3];
  LvArray::tensorOps::copy< 3 >( xLocalMax, m_xLocalMax);
  arrayView3d< localIndex const > const ijkMap = m_ijkMap;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get views to mapping arrays
    arrayView2d< localIndex > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];
    arrayView3d< real64 > const shapeFunctionGradientValues = m_shapeFunctionGradientValues[subRegionIndex];

    // Get particle fields
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();

    // Populate mapping arrays based on particle type
    switch( subRegion.getParticleType() )
    {
      case ParticleType::SinglePoint:
        {
          forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
          {
            localIndex const p = activeParticleIndices[pp];

            computeSinglePointShapeFunctions( gridPosition,
                                              particlePosition[p],
                                              ijkMap,
                                              xLocalMin,
                                              hEl,
                                              mappedNodes[pp],
                                              shapeFunctionValues[pp],
                                              shapeFunctionGradientValues[pp] );
          } );
          break;
        }
      case ParticleType::CPDI:
        {
          arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
          forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
          {
            localIndex const p = activeParticleIndices[pp];

            computeCPDIShapeFunctions( gridPosition,
                                       particlePosition[p],
                                       particleRVectors[p],
                                       ijkMap,
                                       xLocalMin,
                                       hEl,
                                       mappedNodes[pp],
                                       shapeFunctionValues[pp],
                                       shapeFunctionGradientValues[pp] );
          } );
          break;
        }
      default:
        {
          GEOS_ERROR( "Particle type \"" << subRegion.getParticleType() << "\" is not yet supported." );
          break;
        }
    }

    // Increment the subRegion index
    ++subRegionIndex;
  } );
}

void SolidMechanicsMPM::generateNodalNeighborList( ParticleManager & particleManager,
                                                   NodeManager & nodeManager )
{
  // WARNING!!! This will most likely have a race condition running on GPU
  // TODO: Fix it in the future!
  // int const damageFieldPartitioning = m_damageFieldPartitioning;
  // localIndex const numContactGroups = m_numContactGroups;
  // localIndex const numVelocityFields = m_numVelocityFields;
  localIndex const numNodes = nodeManager.size();

  // arrayView2d< real64 const > const gridDamageGradient = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridDamageGradientString() );

  // localIndex numBins = numNodes * numVelocityFields;  
  localIndex numBins = numNodes;
  
  array1d< localIndex > nodeFieldCount( numBins );
  for( localIndex b = 0; b < numBins; ++b )
  {
    nodeFieldCount[b] = 0;
  }
  arrayView1d< localIndex > const nodeFieldCountView = nodeFieldCount;

  // Count number of particles mapped to each node-field pair
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle(); 
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];

    localIndex const numParticlesInSubRegion = subRegion.size(); // Should include ghosted particles
    forAll< serialPolicy >( numParticlesInSubRegion, [=] GEOS_HOST ( localIndex const p )
    {
      // Copy list of mapped nodes into another temporary array
      stdVector< localIndex > uniqueNodes( 8 * numberOfVerticesPerParticle, -1 );
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        uniqueNodes[g] = mappedNodes[p][g];
      }

      // Sort the grid indices and move any duplicates to the end.
      std::ptrdiff_t const numUniqueValues = LvArray::sortedArrayManipulation::makeSortedUnique( uniqueNodes.begin(),
                                                                                                 uniqueNodes.end() );

      for( int g = 0; g < numUniqueValues; ++g )
      {
        localIndex nodeIndex = uniqueNodes[g];
        ++nodeFieldCountView[nodeIndex];
      }
    });
    ++subRegionIndex;
  });

  // GEOS_LOG_RANK("Counted number of bins...");

  m_nodalNeighborList.freeOnDevice(); // just being careful
  m_nodalNeighborList.resizeFromCapacities( nodeFieldCount );

  // arrayView1d< localIndex > const numNeighborsAll = m_nodalNeighborList.m_numParticles.toView();
  ArrayOfArraysView< localIndex > const neighborRegions = m_nodalNeighborList.m_toParticleRegion.toView();
  ArrayOfArraysView< localIndex > const neighborSubRegions = m_nodalNeighborList.m_toParticleSubRegion.toView();
  ArrayOfArraysView< localIndex > const neighborIndices = m_nodalNeighborList.m_toParticleIndex.toView();

  // // Check size of neighborlist after resize
  // for(int b = 0; b < numBins; ++b )
  // {
  //   GEOS_LOG_RANK("b: " <<  b << ", " << 
  //                 "numNeighborsAll[b]: "<< numNeighborsAll[b] << ", " << 
  //                 "neighborRegions[b].size(): "<< neighborRegions[b].size() << ", " << 
  //                 "neighborSubRegions[b].size(): "<< neighborSubRegions[b].size() << ", " <<
  //                 "neighborIndices[b].size(): "<< neighborIndices[b].size() );
  // }

  // Initialize bin count
  array1d< localIndex > binCounts( numBins );
  for( localIndex b = 0; b < numBins; ++b )
  {
    binCounts[b] = 0;
  }
  arrayView1d< localIndex > const binCountsView = binCounts;

  // Populate node particle data into many particle relation (aka nodal particle neighborlist )
  subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    localIndex const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle(); 
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];

    // arrayView1d< localIndex const > const particleGroup = subRegion.getParticleGroup();
    // arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    // arrayView2d< real64 const > const particleDamageGradient = subRegion.getField< fields::mpm::particleDamageGradient >();

    ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
    localIndex regionIndexOfSubRegion = region.getIndexInParent();
    localIndex subRegionIndexInRegion = subRegion.getIndexInParent();

    localIndex const numParticlesInSubRegion = subRegion.size(); // Should include ghosted particles
    forAll< serialPolicy >( numParticlesInSubRegion, [=] GEOS_HOST ( localIndex const p )
    {
      // Copy list of mapped nodes into another temporary array
      stdVector< localIndex > uniqueNodes( 8 * numberOfVerticesPerParticle, -1 );
      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        uniqueNodes[g] = mappedNodes[p][g];
      }

      // Sort the grid indices and move any duplicates to the end.
      std::ptrdiff_t const numUniqueValues = LvArray::sortedArrayManipulation::makeSortedUnique( uniqueNodes.begin(),
                                                                                                 uniqueNodes.end() );

      for( int g = 0; g < numUniqueValues; ++g )
      {
        localIndex nodeIndex = uniqueNodes[g];
        // localIndex fieldIndex = partitionField( numContactGroups,
        //                                         damageFieldPartitioning,
        //                                         particleGroup[p],
        //                                         particleDamageGradient[p],
        //                                         particleSurfaceNormal[p],
        //                                         gridDamageGradient[nodeIndex] ); // Use dfg here
        // localIndex const binIndex = numNodes * fieldIndex + nodeIndex;
    
        // loclaIndex const binCount = binCountsView[binIndex];

        // ++numNeighborsAll[binIndex];
        // neighborRegions[binIndex][binCount] = regionIndexOfSubRegion;
        // neighborSubRegions[binIndex][binCount] = subRegionIndexInRegion;
        // neighborIndices[binIndex][binCount] = p;
        // ++binCountsView[binIndex]; // Will need RAJA atomics
        localIndex const binCount = binCountsView[nodeIndex];

        // ++numNeighborsAll[nodeIndex];
        neighborRegions[nodeIndex][binCount] = regionIndexOfSubRegion;
        neighborSubRegions[nodeIndex][binCount] = subRegionIndexInRegion;
        neighborIndices[nodeIndex][binCount] = p;
        ++binCountsView[nodeIndex]; // Will need RAJA atomics
      }
    });
    ++subRegionIndex;
  });

  // // Check size of neighborlist after resize
  // for(int b = 0; b < numBins; ++b )
  // {
  //   GEOS_LOG_RANK("b: " <<  b << ", " << 
  //                 "numNeighborsAll[b]: "<< numNeighborsAll[b] << ", " << 
  //                 "neighborRegions[b].size(): "<< neighborRegions[b].size() << ", " << 
  //                 "neighborSubRegions[b].size(): "<< neighborSubRegions[b].size() << ", " <<
  //                 "neighborIndices[b].size(): "<< neighborIndices[b].size() );
    
  //   std::stringstream strStream;
  //   strStream << "b: " << b << ", neighbors: ";
  //   for(int n = 0; n < numNeighborsAll[b]; ++n)
  //   {
  //     strStream << "{" << neighborRegions[b][n] << ", " << neighborSubRegions[b][n] << ", " << neighborIndices[b][n] << "}";
  //   }
  //   GEOS_LOG_RANK(strStream.str());
  // }
} 

void SolidMechanicsMPM::computeSinglePointShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                          arraySlice1d< real64 const > const particlePosition,
                                                          arrayView3d< int const > const ijkMap,
                                                          real64 const (& xLocalMin)[3],
                                                          real64 const (& hEl)[3],
                                                          arraySlice1d< int > const mappedNodes,
                                                          arraySlice1d< real64 > const shapeFunctionValues,
                                                          arraySlice2d< real64 > const shapeFunctionGradientValues  )
{
  // get IJK associated with particle center
  int centerIJK[3];
  for( int i=0; i<3; ++i )
  {
    centerIJK[i] = LvArray::math::floor( ( particlePosition[i] - xLocalMin[i] ) / hEl[i] );
  }

  // get node IDs, weights and grad weights
  int node = 0;
  int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
  auto corner_x = gridPosition[corner];

  real64 xRel = (particlePosition[0] - corner_x[0]) / hEl[0];
  real64 yRel = (particlePosition[1] - corner_x[1]) / hEl[1];
  real64 zRel = (particlePosition[2] - corner_x[2]) / hEl[2];

  for( int i=0; i<2; ++i )
  {
    real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
    real64 dxWeight = i/hEl[0] - (1-i)/hEl[0];
    for( int j=0; j<2; ++j )
    {
      real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
      real64 dyWeight = j/hEl[1] - (1-j)/hEl[1];
      for( int k=0; k<2; ++k )
      {
        real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
        real64 dzWeight = k/hEl[2] - (1-k)/hEl[2];
        mappedNodes[node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k];
        shapeFunctionValues[node] = xWeight * yWeight * zWeight;
        shapeFunctionGradientValues[node][0] = dxWeight * yWeight * zWeight;
        shapeFunctionGradientValues[node][1] = xWeight * dyWeight * zWeight;
        shapeFunctionGradientValues[node][2] = xWeight * yWeight * dzWeight;
        ++node;
      }
    }
  }
}

GEOS_DEVICE
GEOS_FORCE_INLINE
void SolidMechanicsMPM::computeCPDIShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                   arraySlice1d< real64 const > const particlePosition,
                                                   arraySlice2d< real64 const > const particleRVectors,
                                                   arrayView3d< int const > const ijkMap,
                                                   real64 const (& xLocalMin)[3],
                                                   real64 const (& hEl)[3],
                                                   arraySlice1d< int > const mappedNodes,
                                                   arraySlice1d< real64 > const shapeFunctionValues,
                                                   arraySlice2d< real64 > const shapeFunctionGradientValues )
{
  int const signs[8][3] = { { -1, -1, -1 },
                            {  1, -1, -1 },
                            {  1,  1, -1 },
                            { -1,  1, -1 },
                            { -1, -1,  1 },
                            {  1, -1,  1 },
                            {  1,  1,  1 },
                            { -1,  1,  1 } };

  real64 alpha[8][8];
  real64 cpdiVolume, oneOverV;
  real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

  for( int i=0; i<3; ++i )
  {
    p_r1[i] = particleRVectors[0][i];
    p_r2[i] = particleRVectors[1][i];
    p_r3[i] = particleRVectors[2][i];
  }

  // We need this because of CPDI domain scaling
  cpdiVolume = 8.0*(-(p_r1[2]*p_r2[1]*p_r3[0]) + p_r1[1]*p_r2[2]*p_r3[0] + p_r1[2]*p_r2[0]*p_r3[1] - p_r1[0]*p_r2[2]*p_r3[1] - p_r1[1]*p_r2[0]*p_r3[2] + p_r1[0]*p_r2[1]*p_r3[2]);
  oneOverV = 1.0 / cpdiVolume;

  alpha[0][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[0][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[0][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
  alpha[1][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[1][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[1][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[2][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[2][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[2][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[3][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[3][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[3][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
  alpha[4][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[4][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[4][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
  alpha[5][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[5][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[5][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[6][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[6][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[6][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[7][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[7][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[7][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);

  // get IJK associated with each corner
  int cornerIJK[8][3]; // CPDI can map to up to 8 cells
  for( int corner=0; corner < 8; ++corner )
  {
    for( int i=0; i<3; ++i )
    {
      real64 cornerPositionComponent = particlePosition[i] +
                                       signs[corner][0] * particleRVectors[0][i] +
                                       signs[corner][1] * particleRVectors[1][i] +
                                       signs[corner][2] * particleRVectors[2][i];

      cornerIJK[corner][i] = LvArray::math::floor( ( cornerPositionComponent - xLocalMin[i] ) / hEl[i] ); // TODO: Temporarily store the CPDI
                                                                                                // corners since they're re-used
                                                                                                // below?
    }
  }

  // get node IDs associated with each corner from IJK map, along with weights and grad weights
  // *** The order in which we access the IJK map must match the order we evaluate the shape functions! ***
  int node = 0;
  for( int corner=0; corner<8; ++corner )
  {
    int cornerNode = ijkMap[cornerIJK[corner][0]][cornerIJK[corner][1]][cornerIJK[corner][2]];
    auto cornerNodePosition = gridPosition[cornerNode];

    real64 x, y, z;
    x = particlePosition[0] + signs[corner][0] * particleRVectors[0][0] + signs[corner][1] * particleRVectors[1][0] + signs[corner][2] * particleRVectors[2][0];
    y = particlePosition[1] + signs[corner][0] * particleRVectors[0][1] + signs[corner][1] * particleRVectors[1][1] + signs[corner][2] * particleRVectors[2][1];
    z = particlePosition[2] + signs[corner][0] * particleRVectors[0][2] + signs[corner][1] * particleRVectors[1][2] + signs[corner][2] * particleRVectors[2][2];

    real64 xRel = (x - cornerNodePosition[0]) / hEl[0];
    real64 yRel = (y - cornerNodePosition[1]) / hEl[1];
    real64 zRel = (z - cornerNodePosition[2]) / hEl[2];

    for( int i=0; i<2; ++i )
    {
      real64 xWeight = i * xRel + (1 - i) * (1.0 - xRel);
      for( int j=0; j<2; ++j )
      {
        real64 yWeight = j * yRel + (1 - j) * (1.0 - yRel);
        for( int k=0; k<2; ++k )
        {
          real64 zWeight = k * zRel + (1 - k) * (1.0 - zRel);
          real64 weight = xWeight * yWeight * zWeight;

          mappedNodes[node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
          shapeFunctionValues[node] = 0.125 * weight;
          shapeFunctionGradientValues[node][0] = alpha[corner][0] * weight;
          shapeFunctionGradientValues[node][1] = alpha[corner][1] * weight;
          shapeFunctionGradientValues[node][2] = alpha[corner][2] * weight;
          ++node;
        }
      }
    }
  }
}

// For on the fly computations
GEOS_DEVICE
GEOS_FORCE_INLINE
void SolidMechanicsMPM::mapNodesAndComputeShapeFunctions( arrayView3d< localIndex const > const ijkMap,
                                                          real64 const (& xLocalMin)[3],
                                                          real64 const (& hEl)[3],
                                                          ParticleType particleType,
                                                          arraySlice1d< real64 const > const particlePosition,
                                                          arraySlice2d< real64 const > const particleRVectors,
                                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                          localIndex * const mappedNodes,
                                                          real64 * const shapeFunctionValues,
                                                          real64 shapeFunctionGradientValues[][3] )
{
    // Populate mapping arrays based on particle type
    switch( particleType )
    {
      case ParticleType::SinglePoint:
        {
          computeSinglePointParticleShapeFunctions( gridPosition,
                                                    particlePosition,
                                                    ijkMap,
                                                    xLocalMin,
                                                    hEl,
                                                    mappedNodes,
                                                    shapeFunctionValues,
                                                    shapeFunctionGradientValues );

          break;
        }
      case ParticleType::CPDI:
        {
          computeCPDIParticleShapeFunctions( gridPosition,
                                             particlePosition,
                                             particleRVectors,
                                             ijkMap,
                                             xLocalMin,
                                             hEl,
                                             mappedNodes,
                                             shapeFunctionValues,
                                             shapeFunctionGradientValues );

          break;
        }

      default:
        {
          GEOS_ERROR( "Particle type \"" << particleType << "\" is not yet supported." );
          break;
        }
    }
}

GEOS_DEVICE
GEOS_FORCE_INLINE
void SolidMechanicsMPM::computeSinglePointParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                                  arraySlice1d< real64 const > const particlePosition,
                                                                  arrayView3d< int const > const ijkMap,
                                                                  real64 const (& xLocalMin)[3],
                                                                  real64 const (& hEl)[3],
                                                                  localIndex * const mappedNodes,
                                                                  real64 * const shapeFunctionValues,
                                                                  real64 shapeFunctionGradientValues[][3] )
{
  // get IJK associated with particle center
  int centerIJK[3];
  for( int i=0; i<3; ++i )
  {
    centerIJK[i] = LvArray::math::floor( ( particlePosition[i] - xLocalMin[i] ) / hEl[i] );
  }

  // get node IDs, weights and grad weights
  int node = 0;
  int corner = ijkMap[centerIJK[0]][centerIJK[1]][centerIJK[2]];
  auto corner_x = gridPosition[corner];

  real64 xRel = (particlePosition[0] - corner_x[0]) / hEl[0];
  real64 yRel = (particlePosition[1] - corner_x[1]) / hEl[1];
  real64 zRel = (particlePosition[2] - corner_x[2]) / hEl[2];

  for( int i=0; i<2; ++i )
  {
    real64 xWeight = i*xRel + (1-i)*(1.0-xRel);
    real64 dxWeight = i/hEl[0] - (1-i)/hEl[0];
    for( int j=0; j<2; ++j )
    {
      real64 yWeight = j*yRel + (1-j)*(1.0-yRel);
      real64 dyWeight = j/hEl[1] - (1-j)/hEl[1];
      for( int k=0; k<2; ++k )
      {
        real64 zWeight = k*zRel + (1-k)*(1.0-zRel);
        real64 dzWeight = k/hEl[2] - (1-k)/hEl[2];
        mappedNodes[node] = ijkMap[centerIJK[0]+i][centerIJK[1]+j][centerIJK[2]+k];
        shapeFunctionValues[node] = xWeight * yWeight * zWeight;
        shapeFunctionGradientValues[node][0] = dxWeight * yWeight * zWeight;
        shapeFunctionGradientValues[node][1] = xWeight * dyWeight * zWeight;
        shapeFunctionGradientValues[node][2] = xWeight * yWeight * dzWeight;
        ++node;
      }
    }
  }
}

GEOS_DEVICE
GEOS_FORCE_INLINE
void SolidMechanicsMPM::computeCPDIParticleShapeFunctions( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition,
                                                           arraySlice1d< real64 const > const particlePosition,
                                                           arraySlice2d< real64 const > const particleRVectors,
                                                           arrayView3d< int const > const ijkMap,
                                                           real64 const (& xLocalMin)[3],
                                                           real64 const (&hEl)[3],
                                                           localIndex * const mappedNodes,
                                                           real64 * const shapeFunctionValues,
                                                           real64 shapeFunctionGradientValues[][3] )
{
  int const signs[8][3] = { { -1, -1, -1 },
                            {  1, -1, -1 },
                            {  1,  1, -1 },
                            { -1,  1, -1 },
                            { -1, -1,  1 },
                            {  1, -1,  1 },
                            {  1,  1,  1 },
                            { -1,  1,  1 } };

  real64 alpha[8][8];
  real64 cpdiVolume, oneOverV;
  real64 p_r1[3], p_r2[3], p_r3[3]; // allowing 1-indexed r-vectors to persist to torture future postdocs >:)

  for( int i=0; i<3; ++i )
  {
    p_r1[i] = particleRVectors[0][i];
    p_r2[i] = particleRVectors[1][i];
    p_r3[i] = particleRVectors[2][i];
  }

  // We need this because of CPDI domain scaling
  cpdiVolume = 8.0*(-(p_r1[2]*p_r2[1]*p_r3[0]) + p_r1[1]*p_r2[2]*p_r3[0] + p_r1[2]*p_r2[0]*p_r3[1] - p_r1[0]*p_r2[2]*p_r3[1] - p_r1[1]*p_r2[0]*p_r3[2] + p_r1[0]*p_r2[1]*p_r3[2]);
  oneOverV = 1.0 / cpdiVolume;

  alpha[0][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[0][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[0][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
  alpha[1][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[1][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[1][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[2][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[2][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[2][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[3][0]=oneOverV*(p_r1[2]*p_r2[1] - p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[3][1]=oneOverV*(-(p_r1[2]*p_r2[0]) + p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[3][2]=oneOverV*(p_r1[1]*p_r2[0] - p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
  alpha[4][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[4][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[4][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);
  alpha[5][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] - p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] + p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[5][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] + p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] - p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[5][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] - p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] + p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[6][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] - p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] + p_r2[1]*p_r3[2]);
  alpha[6][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] + p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] - p_r2[0]*p_r3[2]);
  alpha[6][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] - p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] + p_r2[0]*p_r3[1]);
  alpha[7][0]=oneOverV*(-(p_r1[2]*p_r2[1]) + p_r1[1]*p_r2[2] + p_r1[2]*p_r3[1] + p_r2[2]*p_r3[1] - p_r1[1]*p_r3[2] - p_r2[1]*p_r3[2]);
  alpha[7][1]=oneOverV*(p_r1[2]*p_r2[0] - p_r1[0]*p_r2[2] - p_r1[2]*p_r3[0] - p_r2[2]*p_r3[0] + p_r1[0]*p_r3[2] + p_r2[0]*p_r3[2]);
  alpha[7][2]=oneOverV*(-(p_r1[1]*p_r2[0]) + p_r1[0]*p_r2[1] + p_r1[1]*p_r3[0] + p_r2[1]*p_r3[0] - p_r1[0]*p_r3[1] - p_r2[0]*p_r3[1]);

  // get IJK associated with each corner
  int cornerIJK[8][3]; // CPDI can map to up to 8 cells
  for( int corner=0; corner < 8; ++corner )
  {
    for( int i=0; i<3; ++i )
    {
      real64 cornerPositionComponent = particlePosition[i] +
                                       signs[corner][0] * particleRVectors[0][i] +
                                       signs[corner][1] * particleRVectors[1][i] +
                                       signs[corner][2] * particleRVectors[2][i];

      cornerIJK[corner][i] = LvArray::math::floor( ( cornerPositionComponent - xLocalMin[i] ) / hEl[i] ); // TODO: Temporarily store the CPDI
                                                                                                // corners since they're re-used
                                                                                                // below?
    }
  }

  // get node IDs associated with each corner from IJK map, along with weights and grad weights
  // *** The order in which we access the IJK map must match the order we evaluate the shape functions! ***
  int node = 0;
  for( int corner=0; corner<8; ++corner )
  {
    int cornerNode = ijkMap[cornerIJK[corner][0]][cornerIJK[corner][1]][cornerIJK[corner][2]];
    auto cornerNodePosition = gridPosition[cornerNode];

    real64 x, y, z;
    x = particlePosition[0] + signs[corner][0] * particleRVectors[0][0] + signs[corner][1] * particleRVectors[1][0] + signs[corner][2] * particleRVectors[2][0];
    y = particlePosition[1] + signs[corner][0] * particleRVectors[0][1] + signs[corner][1] * particleRVectors[1][1] + signs[corner][2] * particleRVectors[2][1];
    z = particlePosition[2] + signs[corner][0] * particleRVectors[0][2] + signs[corner][1] * particleRVectors[1][2] + signs[corner][2] * particleRVectors[2][2];

    real64 xRel = (x - cornerNodePosition[0]) / hEl[0];
    real64 yRel = (y - cornerNodePosition[1]) / hEl[1];
    real64 zRel = (z - cornerNodePosition[2]) / hEl[2];

    for( int i=0; i<2; ++i )
    {
      real64 xWeight = i * xRel + (1 - i) * (1.0 - xRel);
      for( int j=0; j<2; ++j )
      {
        real64 yWeight = j * yRel + (1 - j) * (1.0 - yRel);
        for( int k=0; k<2; ++k )
        {
          real64 zWeight = k * zRel + (1 - k) * (1.0 - zRel);
          real64 weight = xWeight * yWeight * zWeight;

          mappedNodes[node] = ijkMap[cornerIJK[corner][0]+i][cornerIJK[corner][1]+j][cornerIJK[corner][2]+k];
          shapeFunctionValues[node] = 0.125 * weight;
          shapeFunctionGradientValues[node][0] = alpha[corner][0] * weight;
          shapeFunctionGradientValues[node][1] = alpha[corner][1] * weight;
          shapeFunctionGradientValues[node][2] = alpha[corner][2] * weight;
          ++node;
        }
      }
    }
  }
}

//CC: Either need to return or pass body force variables by reference
inline void SolidMechanicsMPM::computeGeneralizedVortexMMSBodyForce( real64 const time_n,
                                                                                 ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // Method of Manufactured Solutions - generalized vortex.
  // As Described in K. Kamojjala, R Brannon, A Sadeghirad, J. Guilkey "Verification Tests
  // in Solid Mechanics", Engineering with computers, (2015).

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    // ParticleType particleType = subRegion.getParticleType();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();

    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    // CC: this feels wrong
    // Constitutive model for body force calculation
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );

    array1d< real64 > shearModulus;
    string constitutiveModelName = constitutiveModel.getCatalogName();
    if( constitutiveModelName == "Hyperelastic" ){
      Hyperelastic & hyperelastic = dynamic_cast< Hyperelastic & >( constitutiveModel );
      shearModulus = hyperelastic.shearModulus();
    }

    if( constitutiveModelName == "HyperelasticMMS" ){
      HyperelasticMMS & hyperelasticMMS = dynamic_cast< HyperelasticMMS & >( constitutiveModel );
      shearModulus = hyperelasticMMS.shearModulus();
    }

    if( constitutiveModelName == "ElasticIsotropic" || constitutiveModelName == "CeramicDamage" || constitutiveModelName == "StrainHardeningPolymer" || constitutiveModelName == "VonMisesJ" ){
      ElasticIsotropic & elasticIsotropic = dynamic_cast< ElasticIsotropic & >( constitutiveModel );
      shearModulus = elasticIsotropic.shearModulus();
    }

    if( constitutiveModelName == "ElasticTransverseIsotropic" || constitutiveModelName == "ElasticTransverseIsotropicPressureDependent" ){
      ElasticTransverseIsotropic & elasticTransverseIsotropic = dynamic_cast< ElasticTransverseIsotropic & >( constitutiveModel );
      shearModulus = elasticTransverseIsotropic.effectiveShearModulus();
    }

    if( constitutiveModelName == "Graphite" )
    {
      Graphite & graphite = dynamic_cast< Graphite & >( constitutiveModel );
      shearModulus = graphite.effectiveShearModulus();
    }

    GEOS_ERROR_IF( !constitutiveModel.hasWrapper( constitutive::SolidBase::viewKeyStruct:: defaultDensityString() ) , "Constitutive model must have particle density for the generalized vortex problem!");
    real64 const initialDensity = constitutiveModel.getReference< real64 >( constitutive::SolidBase::viewKeyStruct:: defaultDensityString() );

    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      real64 mu = shearModulus[p],
             rho0 = initialDensity, // CC: WRONG, check if this is supposed to be density
             //mu = 384.615384615384585, // Make consistent with elastic properties
             //lambda = 577,             // Make consistent with elastic properties
             //rho0 = 1000.0,            // Make consistent with initial density.
             A = 1.0, // Amplitude of the deformation.
             ri = 0.75,
             ro = 1.25;

      real64 pi = 3.141592653589793;

      real64 x = particlePosition[p][0],
             y = particlePosition[p][1];

      // Radial and angular coordinates:
      real64 R = sqrt( x * x + y * y );

      // CC: Comment from Mike in old geos?
      // I believe the examples in the paper used the reference radius to compute the body force, which
      // greatly increases accuracy.  I've disabled this option, because it is better test to use
      // the reference in the current configuration.
      //    real64 x0 = p_x0( 0 ),
      //           y0 = p_x0( 1 );
      //    real64 R = sqrt( x0 * x0 + y0 * y0 );

      // In either case use the current angle
      real64 Theta = atan2( y, x );

      if( ( R > ri ) && ( R < ro ) )
      {
        // Evaluate some temporary variables:
        real64 p1 = 4096.0 * R * LvArray::math::pow( 15.0 - 47.0 * R + 48.0 * R * R - 16.0 * R * R * R, 2 ) * mu * LvArray::math::pow(
            sin( pi * time_n ), 4 ) / rho0,
            p2 = pi * pi * R * LvArray::math::pow( 15.0 - 32.0 * R + 16.0 * R * R, 4 ) * pow( sin( 2.0 * pi * time_n ), 2 ),
            p3 = -16.0 * ( -45.0 + 188.0 * R - 240.0 * R * R + 96.0 * R * R * R ),
            p4 = -45.0 + 188.0 * R - 240.0 * R * R + 96.0 * R * R * R,
            p5 = LvArray::math::pow( 15.0 - 32.0 * R + 16.0 * R * R, 2 );

        // Evaluate radial component of the body force:
        real64 br = p1 - p2;

        // Evaluate circumferential component of the body force:
        real64 bt = ( 2.0 * mu * p3 + 2.0 * cos( 2.0 * pi * time_n ) * ( 16.0 * mu * p4 + pi * pi * R * rho0 * p5 ) ) / rho0;

        // Evaluate the rotation angle
        real64 alpha = A * ( 1.0 - cos( 2.0 * pi * time_n ) ) * ( 1.0 - 32.0 * LvArray::math::pow( R - 1.0, 2 ) + 256.0 * LvArray::math::pow(
                        R - 1.0, 4 ) ) / 2.0;

        // Evaluate the deformed angular coordinate
        real64 theta = Theta + alpha;

        // evaluate the cartesian components of the body force:
        particleBodyForce[p][0] += br * cos( theta ) - bt * sin( theta );
        particleBodyForce[p][1] += br * sin( theta ) + bt * cos( theta );
        particleBodyForce[p][2] += 0.0;
      }
      } ); // particle loop
      ++subRegionIndex;
  } ); // subregion loop

}


void SolidMechanicsMPM::computeBodyForce( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Particle fields
    arrayView2d< real64 > const particleBodyForce = subRegion.getField< fields::mpm::particleBodyForce >();

    int const numDims = m_numDims;
    real64 bodyForce[3] = {};
    LvArray::tensorOps::copy< 3 >( bodyForce, m_bodyForce );
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      for( int i=0; i < numDims; ++i)
      {
        particleBodyForce[p][i] = bodyForce[i];
      }
    } ); // particle loop

    ++subRegionIndex;
  } ); // subregion loop
}

void SolidMechanicsMPM::computeArtificialViscosity( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  real64 length = m_planeStrain == 1 ? LvArray::math::min( m_hEl[0], m_hEl[1] ) : LvArray::math::min( m_hEl[0], LvArray::math::min( m_hEl[1], m_hEl[2] ) );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 const > const particleWavespeed = subRegion.getField< fields::mpm::particleWavespeed >();
    arrayView1d< real64 const > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
    arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();


    real64 const artificialViscosityQ0 = m_artificialViscosityQ0;
    real64 const artificialViscosityQ1 = m_artificialViscosityQ1;
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // This is computed for du/dx = tr(D) = tr(L);
      real64 trD = LvArray::tensorOps::trace< 3 >( particleVelocityGradient[p] ); // characteristic increment = dx*(du/dx), where du/dx is the characteristic velocity gradient.
      if( trD < 0.0 )
      {
        particleArtificialViscosity[p] = particleDensity[p] * ( artificialViscosityQ1 * ( trD * length ) * ( trD * length )
                                         + artificialViscosityQ0 * fabs( particleWavespeed[p] * trD * length ) );
      }
      else
      {
        particleArtificialViscosity[p] = 0.0;
      }
    } );
  } );
}


void SolidMechanicsMPM::computeInternalEnergyAndTemperature( const real64 dt,
                                                             ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int numDims = m_numDims;
  int shockHeating = m_shockHeating;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get constitutive model reference
    // string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    // ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );
    // if( constitutiveModel.hasWrapper( "strengthScale" ) )
    // {
    // arrayView1d< real64 > const constitutiveStrengthScale = constitutiveModel.getReference< array1d< real64 > >( "strengthScale" );

      arrayView1d< real64 const > const particleDensity = subRegion.getField< fields::mpm::particleDensity >();
      arrayView1d< real64 const > const particleHeatCapacity = subRegion.getField< fields::mpm::particleHeatCapacity >();
      arrayView1d< real64 const > const particleArtificialViscosity = subRegion.getField< fields::mpm::particleArtificialViscosity >();
      arrayView2d< real64 const > const particleStress = subRegion.getField< fields::mpm::particleStress >();
      arrayView3d< real64 const > const particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();

      arrayView1d< real64 > const particleTemperature = subRegion.getParticleTemperature();
      arrayView1d< real64 > const particleInternalEnergy = subRegion.getField< fields::mpm::particleInternalEnergy >();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        //  ========================================================================
        //  Update internal energy with stress power, where Fdot and PK1 stress are a work-
        //  Compute Stress and strain measure The rate of internal energy: udot = (1/rho0)*PK1:Fdot
        //  use L=Fdot*Finv, rho*J=rho0, and PK1=J*sigma(i,k)*Finv(j,k)
        //          udot = (1/rho0)*PK1:Fdot
        //          udot = (1/rho0)*(J*sigma(i,k)*Finv(j,k))*Fdot(i,j)
        //          udot = (1/rho0)*J*sigma(i,k)*L(i,k)
        //          udot = (1/rho)*sigma(i,k)*L(i,k)
        //
        // To compute the temperature it will be necessary to partition the energy into a thermal part
        // and to use some model for the heat capacity.

        real64 energyIncrement = 0.0;

        // 1st half-step internal energy update
        localIndex voigtMap[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };
        for( int i = 0 ; i < numDims ; ++i )
        {
          for( int j = 0 ; j < numDims ; ++j )
          {
            energyIncrement += 0.5 * dt * particleStress[p][voigtMap[i][j]] * particleVelocityGradient[p][i][j] / particleDensity[p]; //  General case for increment in internal energy
          }
          energyIncrement -= shockHeating * 0.5 * dt * particleArtificialViscosity[p] * particleVelocityGradient[p]( i, i ) / particleDensity[p]; // Artificial viscosity contribution
        }

        particleTemperature[p] += energyIncrement / particleHeatCapacity[p];
        particleInternalEnergy[p] += energyIncrement;

        //      // 1st half-step internal energy update
        //      for( int i = 0 ; i < m_ndim ; ++i )
        //      {
        //        for( int j = 0 ; j < m_ndim ; ++j )
        //        {
        //          p_internalEnergy[pp] += 0.5 * dt * p_stress[pp]( i, j ) * p_L[pp]( i, j ) / p_rho[pp]; //  General case for increment in internal energy
        //        }
        //        p_internalEnergy[pp] -= m_shock_heating * 0.5 * dt * p_q[pp] * p_L[pp]( i, i ) / p_rho[pp]; // Artificial viscosity contribution
        //      }

        // Prevent non-physical negative internal energy from integration error.
        particleInternalEnergy[p] = ( particleInternalEnergy[p] < 0 ) ? 0 : particleInternalEnergy[p];

        // ***************************************************************************
        // TODO: This is the internal energy calculation from Uintah:
        // We need to consider whether to include the artificial viscosity heating.
        //
        //    R2Tensor avgStress = (pStress_new[idx] + pStress[idx])*0.5;
        //    double avgVolume  = (pVolume_deformed[idx]+pVolume[idx])*0.5;
        //
        //    double pSpecificStrainEnergy = (tensorD(0,0)*avgStress(0,0) +
        //                                    tensorD(1,1)*avgStress(1,1) +
        //                                    tensorD(2,2)*avgStress(2,2) +
        //                               2.0*(tensorD(0,1)*avgStress(0,1) +
        //                                    tensorD(0,2)*avgStress(0,2) +
        //                                    tensorD(1,2)*avgStress(1,2)))*
        //                                    avgVolume*delT/pMass[idx];
        //
        //    // Compute rate of change of specific volume
        //    double Vdot = (pVolume_deformed[idx] - pVolume[idx])/(pMass[idx]*delT);
        //
        //    pEnergy_new[idx] = pEnergy[idx] + pSpecificStrainEnergy
        //                                    - p_q[idx]*Vdot*delT*include_AV_heating;
        //
        //    totalStrainEnergy += pSpecificStrainEnergy*pMass[idx];
        // ***************************************************************************

      } );
    // }
  } );
}

void SolidMechanicsMPM::computeSPHJacobian( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  int const planeStrain = m_planeStrain;
  int const jacobianType = m_computeSPHJacobian;
  real64 hEl[3];
  LvArray::tensorOps::copy< 3 >( hEl, m_hEl );

  // Compute local density based on SPH method (actually compute the local jacobian to account
  // for the possibility of different densities between materials.
  ParticleManager::ParticleViewAccessor< arrayView2d< real64 const > > particlePositionAccessor = particleManager.constructArrayViewAccessor< real64, 2 >( "particleCenter" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleReferenceVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleReferenceVolume" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleVolumeAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleVolume" );
  ParticleManager::ParticleViewAccessor< arrayView1d< real64 const > > particleMassAccessor = particleManager.constructArrayViewAccessor< real64, 1 >( "particleMass" );

  ParticleManager::ParticleViewConst< arrayView2d< real64 const > > particlePositionView = particlePositionAccessor.toNestedViewConst();
  // ParticleManager::ParticleViewConst< arrayView1d< real64 const > > particleReferenceVolumeView = particleReferenceVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > particleVolumeView = particleVolumeAccessor.toNestedViewConst();
  ParticleManager::ParticleViewConst< arrayView1d< real64 const > > particleMassView = particleMassAccessor.toNestedViewConst();

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Get neighbor list
    OrderedVariableToManyParticleRelation & neighborList = subRegion.neighborList();
    arrayView1d< localIndex const > const numNeighborsAll = neighborList.m_numParticles.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborRegions = neighborList.m_toParticleRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborSubRegions = neighborList.m_toParticleSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const neighborIndices = neighborList.m_toParticleIndex.toViewConst();

    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView1d< real64 const > const particleReferenceVolume = subRegion.getField< fields::mpm::particleReferenceVolume >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleSPHJacobian = subRegion.getField< fields::mpm::particleSPHJacobian >();

    ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
    localIndex regionIndexOfSubRegion = region.getIndexInParent();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      localIndex numNeighbors = numNeighborsAll[p];
      arraySlice1d< localIndex const > const regionIndices = neighborRegions[pp];
      arraySlice1d< localIndex const > const subRegionIndices = neighborSubRegions[pp];
      arraySlice1d< localIndex const > const particleIndices = neighborIndices[pp];

      // Compute the kernel scalar field at a point, for a given list of neighbor particles.
      // The lists xp, fp, and the length np could refer to all the particles in the patch,
      // but generally this function will be evaluated with x equal to some particle center,
      // and xp, fp, will be lists for just the neighbors of the particle.

      // Initialize
      real64 relativePosition[3];
      real64 f = 0.0,
             k = 0.0;

      if( jacobianType == 1 || jacobianType == 2 )
      {
        for(localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
        {
          localIndex regionIndex = regionIndices[neighborIndex];
          localIndex subRegionIndex = subRegionIndices[neighborIndex];
          localIndex particleIndex = particleIndices[neighborIndex];

          real64 neighborVolume = 0.0;
          real64 fp = 0.0;
          switch( jacobianType )
          {
            case 1:
              fp = ( planeStrain == 1 ? hEl[2] : 1.0 ) / particleVolumeView[regionIndex][subRegionIndex][particleIndex];
              neighborVolume = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
              break;
            case 2:
              if( regionIndex == regionIndexOfSubRegion )
              {
                fp = particleMassView[regionIndex][subRegionIndex][particleIndex] / particleVolumeView[regionIndex][subRegionIndex][particleIndex];
                neighborVolume = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
              }
              break;
            default:
              //Throw error
              break;
          }

          LvArray::tensorOps::copy< 3 >( relativePosition, particlePosition[p]);
          LvArray::tensorOps::subtract< 3 >( relativePosition, particlePositionView[regionIndex][subRegionIndex][particleIndex] );
          real64 kernelVal = kernel( LvArray::tensorOps::l2Norm< 3 >( relativePosition ) );
          k += neighborVolume * kernelVal;
          f += neighborVolume * fp * kernelVal;
        }

        // Return the normalized kernel field (which eliminates edge effects)
        real64 sphDensity = k > 0.0 ? f / k : 0.0;

        // type 1
        real64 initialParticleNumberDensity = 1.0 / particleReferenceVolume[p];
        if( planeStrain == 1 && jacobianType == 1 )
        {
          // This is actually the particles per unit area.  This assumes that there is only one particles in
          // the z-direction, which will only be the case if the particleFileWriter script is called
          // with the python flag planeStrain=1 (in addition to setting that GEOS flag)
          initialParticleNumberDensity *= hEl[2];
        }

        //type 2
        if( jacobianType == 2 )
        {
          initialParticleNumberDensity *= particleMass[p];
        }

        particleSPHJacobian[p] = initialParticleNumberDensity / fmax( 1.0e-9, sphDensity );
    }
    else
    {
      particleSPHJacobian[p] = particleVolume[p]/particleReferenceVolume[p];
    }
      // // Populate neighbor positions
      // for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
      // {
      //   localIndex regionIndex = regionIndices[neighborIndex];
      //   localIndex subRegionIndex = subRegionIndices[neighborIndex];
      //   localIndex particleIndex = particleIndices[neighborIndex];
      //   neighborPositions[neighborIndex][0] = particlePositionView[regionIndex][subRegionIndex][particleIndex][0];
      //   neighborPositions[neighborIndex][1] = particlePositionView[regionIndex][subRegionIndex][particleIndex][1];
      //   neighborPositions[neighborIndex][2] = particlePositionView[regionIndex][subRegionIndex][particleIndex][2];
      // }

      // // Compare SPH number density to initial particle-wise number density. This gives poor results for
      // // non-unform particle spacing and/or local particle refinement.
      // if( m_computeSPHJacobian == 1 )
      // {
      //   for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
      //   {
      //     localIndex regionIndex = regionIndices[neighborIndex];
      //     localIndex subRegionIndex = subRegionIndices[neighborIndex];
      //     localIndex particleIndex = particleIndices[neighborIndex];
      //     neighborNumberDensities[neighborIndex] = ( planeStrain ? hEl[2] : 1.0 ) / particleVolumeView[regionIndex][subRegionIndex][particleIndex];
      //     neighborVolumes[neighborIndex] = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
      //   }

      //   // TODO modify this function to mirror points near symmetry planes (within neighbor radius of the boundary)
      //   // to obtain good behavior at sym boundaries.

      //   // Evaluate kernel field
      //   real64 sphParticleNumberDensity = computeKernelField( particlePosition[p], // query point
      //                                                         numNeighbors,
      //                                                         neighborPositions, // List of neighbor particle locations
      //                                                         neighborVolumes,
      //                                                         neighborNumberDensities );

      //   // Determine J
      //   real64 initialParticleNumberDensity;
      //   if( planeStrain )
      //   {
      //     // This is actually the particles per unit area.  This assumes that there is only one particles in
      //     // the z-direction, which will only be the case if the particleFileWriter script is called
      //     // with the python flag planeStrain=1 (in addition to setting that GEOS flag)
      //     initialParticleNumberDensity = m_hEl[2] / particleReferenceVolume[p];
      //   }
      //   else
      //   {
      //     initialParticleNumberDensity = 1.0 / particleReferenceVolume[p];
      //   }
      //   particleSPHJacobian[p] = initialParticleNumberDensity / fmax( 1.0e-9, sphParticleNumberDensity );
      // }
      // // Use SPH mass density to determine J. Only neighbors that are the same material as the ppth particle are considered.
      // // This approach is more robust against non-uniform spacing/local particle refinement.
      // else if ( m_computeSPHJacobian == 2 )
      // {
      //   for( localIndex neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++ )
      //   {
      //     localIndex regionIndex = regionIndices[neighborIndex];
      //     if( regionIndex == regionIndexOfSubRegion ) // same material
      //     {
      //       localIndex subRegionIndex = subRegionIndices[neighborIndex];
      //       localIndex particleIndex = particleIndices[neighborIndex];
      //       neighborMassDensities[neighborIndex] = particleMassView[regionIndex][subRegionIndex][particleIndex] / particleVolumeView[regionIndex][subRegionIndex][particleIndex]; // CC: Should we just used particle density directly?
      //       neighborVolumes[neighborIndex] = particleVolumeView[regionIndex][subRegionIndex][particleIndex];
      //     }
      //     else
      //     {
      //       neighborMassDensities[neighborIndex] = 0.0;
      //       neighborVolumes[neighborIndex] = 0.0;
      //     }
      //   }

      //   // Evaluate kernel field
      //   real64 sphMassDensity = computeKernelField( particlePosition[p], // query point
      //                                               numNeighbors,
      //                                               neighborPositions,    // List of neighbor particle locations
      //                                               neighborVolumes,
      //                                               neighborMassDensities );

      //   particleSPHJacobian[p] = ( particleMass[p] / particleReferenceVolume[p] ) / fmax( 1.0e-9, sphMassDensity );
      // }
      // else
      // {
      //   particleSPHJacobian[p] = particleVolume[p]/particleReferenceVolume[p];
      // }
    } );
  } );
}


void SolidMechanicsMPM::sphOverlapCorrection( real64 const dt,
                                              ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;
  int const numDims = m_numDims;
  int const planeStrain = m_planeStrain;
  real64 const overlapThreshold1 = m_overlapThreshold1;
  real64 const overlapThreshold2 = m_overlapThreshold2;

  // Compute the grid field gradient as the max of the gradients of the particles mapping to the node
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView3d< real64 > particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > particleFDot = subRegion.getField< fields::mpm::particleFDot >();
    arrayView3d< real64 > particleVelocityGradient = subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 const > const particleSPHJacobian = subRegion.getField< fields::mpm::particleSPHJacobian >();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      real64 Fold[3][3] = {};
      LvArray::tensorOps::copy< 3, 3 >( Fold, particleDeformationGradient[p] );

      // Update the Jacobian based on particle data.
      real64 J = LvArray::tensorOps::determinant< 3 >( Fold );

      if ( J <= 0.0 )
      {
        return; // CC: For lambda we would return instead of continue?
      }

      // If there is overdensification, the jacobian as computed from the sph kernel will be much less
      // than that from the particle def. grad, so overlap>1 may indicate overdensification. But, since
      // the neighbor radius is small, there is noise in J_sph, so we only want to scale F if the
      // overlap is above a threshold.
      real64 overlap = J / particleSPHJacobian[p];

      // If the densification exceeds 1, slowly transition from the particle J to the non-local SPH J,
      // and scale Fij accordingly.  Perform this so that J=Jsph if J/J_sph exceeds a threshold.
      // Note that at free surfaces, undeformed Jsph will be >1. This is also true at sym boundaries
      // until we can fix the calculation of the kernel field at boundaries.
      if( overlap > overlapThreshold1 )
      {
        real64 alpha = 1.0; // Smoothing parameter: Jnew = (1-alpha)*Jold + alpha*Jsph
        if( overlap <= overlapThreshold2 )
        {
          // smooth step function
          real64 xi = ( overlap - overlapThreshold1 ) / ( overlapThreshold2 - overlapThreshold1 );
          alpha = 3.0 * xi * xi - 2.0 * xi * xi * xi;
        }

        real64 scale;
        if( planeStrain == 1 )
        {
          scale = LvArray::math::sqrt( 1.0 + alpha * ( 1.0 / overlap - 1.0 ) );
        }
        else
        {
          scale = LvArray::math::pow( 1.0 + alpha * ( 1.0 / overlap - 1.0 ), 1.0 / 3.0 );
        }

        // put a limit on scale so it can't change the density too much in a single time-step.
        scale = fmax( 0.99, scale );

        // TODO modify p_L as well, consistent with the change to p_F in case a hypoelastic constitutive model
        // is used, and so the update to the internal energy is consistent with the change in F.
        //
        // The simple attempt, commented out below, produces extremely high internal energy when
        // the overlap correction is used in a DFG brittle damage problem.

        for( int i = 0 ; i < numDims; ++i )
        {
          for( int j = 0 ; j < numDims; ++j )
          {
            particleDeformationGradient[p][i][j] *= scale;
            particleFDot[p][i][j] += ( particleDeformationGradient[p][i][j] - Fold[i][j] ) / dt;

            // Update velocity gradient without taking the inverse of the particle deformation gradient
            particleVelocityGradient[p][i][j] /= scale;
          }

          particleVelocityGradient[p][i][i] += ( scale - 1.0 ) / ( scale * dt );
        }
        // Modify p_L as well, consistent with the change to p_F in case a hypoelastic constitutive model
        // is used, and so the update to the internal energy is consistent with the change in F.

        // real64 invF[3][3] = {};
        // LvArray::tensorOps::invert< 3 >( invF, particleDeformationGradient[p] );
        // LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleVelocityGradient[p], particleFDot[p], invF );
      }
    } );
  } );

}

void SolidMechanicsMPM::resetDeformationGradient( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION
  int const planeStrain = m_planeStrain;
  int const resetDefGradForFullyDamagedParticles = m_resetDefGradForFullyDamagedParticles;
  int const cpdiDomainScaling = m_cpdiDomainScaling;

  // Reset the deformation gradient to be the spherical part.
  // This should only be used for cases where the deviatoric part
  // of the constitutive model is hypoelastic.
  //
  // We also keep to rotation to make the plotting look a bit better.
  //
  // We also use this for the gas constitutive model.
  // TODO: Make this a more general option (per material type)
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    const ContinuumBase & constitutiveModel = getConstitutiveModel< ContinuumBase >( subRegion, solidMaterialName );
    string constitutiveModelName = constitutiveModel.getCatalogName();

    arrayView1d< real64 const > const particleDamage = subRegion.getParticleDamage();

    arrayView1d< int const > particleDomainScaledFlag;
    if( cpdiDomainScaling == 1 )
    {
      particleDomainScaledFlag = subRegion.getField< fields::mpm::particleDomainScaledFlag >();
    }

    if( constitutiveModelName == "Gas" || resetDefGradForFullyDamagedParticles == 1 )
    {
      arrayView3d< real64 > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];
        if( constitutiveModelName == "Gas" || ( particleDamage[p] > 0.9999999 && ( cpdiDomainScaling == 1 && particleDomainScaledFlag[p] == 1 ) ) )
        {
          real64 rotation[3][3] = {};
          real64 deformationGradient[3][3] = {};
          LvArray::tensorOps::copy< 3, 3 >( deformationGradient, particleDeformationGradient[p] );
          LvArray::tensorOps::polarDecomposition< 3 >( rotation, deformationGradient ); // Start-of-step polar decomposition

          real64 J = LvArray::tensorOps::determinant< 3 >( particleDeformationGradient[p]);

          real64 U[3][3];
          if ( planeStrain == 1 )
          {
            real64 JtoOneHalf = LvArray::math::sqrt( J );
            U[0][0] = JtoOneHalf;
            U[1][1] = JtoOneHalf;
            U[2][2] = 1.0;
          }
          else
          {
            real64 JtoOneThird = LvArray::math::pow( J, 1.0 / 3.0 );
            U[0][0] = JtoOneThird;
            U[1][1] = JtoOneThird;
            U[2][2] = JtoOneThird;
          }

          LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( particleDeformationGradient[p], rotation, U);
        }
      } );
    }
  } );
}


void SolidMechanicsMPM::unscaleCPDIVectors( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  // This undoes CPDI domain scaling. Can be used to plot unscaled domains.
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Only perform on CPDI particles
    if( subRegion.getParticleType() == ParticleType::CPDI )
    {
      arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();
      arrayView3d< real64 const > const particleReferenceRVectors = subRegion.getField< fields::mpm::particleReferenceRVectors >();
      arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();

      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
      forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][0], particleDeformationGradient[p], particleReferenceRVectors[p][0] );
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][1], particleDeformationGradient[p], particleReferenceRVectors[p][1] );
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( particleRVectors[p][2], particleDeformationGradient[p], particleReferenceRVectors[p][2] );
      } );
    }
  } );
}

void SolidMechanicsMPM::computeKineticEnergy( ParticleManager & particleManager )
{
  GEOS_MARK_FUNCTION;

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< real64 > const particleKineticEnergy = subRegion.getField< fields::mpm::particleKineticEnergy >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< parallelDevicePolicy<> >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // TODO: add micro kinetic energy term based on L?
      particleKineticEnergy[p] = 0.5 * particleMass[p] * LvArray::tensorOps::l2NormSquared< 3 >( particleVelocity[p] );
    } );
  } );
}

// This function will do sums of nodal data to write a profile to a file of the
// state averaged across the y-z plane at each x-nodal position.  This is very slow
// but is useful for on-the-fly extraction of shock profiles from mesoscale simulations.
// Typically this is not done every time step.
// Sum local arrays to compute x-profile to write to file
void SolidMechanicsMPM::computeXProfile( int const cycleNumber,
                                         real64 const time,
                                         real64 const dt,
                                         NodeManager & nodeManager,
                                         SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;

  int numXNodesLocal = m_nEl[0] + 1;
  // int numYNodesLocal = m_nEl[0] + 1;
  // int numZNodesLocal = m_nEl[0] + 1;

  int numXNodesGlobal = static_cast< int >( round( ( m_xGlobalMax[0] - m_xGlobalMin[0]) / m_hEl[0] ) ) + 1;
  int numYNodesGlobal = static_cast< int >( round( ( m_xGlobalMax[1] - m_xGlobalMin[1]) / m_hEl[1] ) ) + 1;
  int numZNodesGlobal = static_cast< int >( round( ( m_xGlobalMax[2] - m_xGlobalMin[2]) / m_hEl[2] ) ) + 1;

  array1d< real64 > xProfileDensityLocal( numXNodesGlobal ),
                    xProfileDensityGlobal( numXNodesGlobal ),
                    xProfileDamageLocal( numXNodesGlobal ),
                    xProfileDamageGlobal( numXNodesGlobal ),
                    xProfileXStressLocal( numXNodesGlobal ),
                    xProfileXStressGlobal( numXNodesGlobal ),
                    xProfileYStressLocal( numXNodesGlobal ),
                    xProfileYStressGlobal( numXNodesGlobal ),
                    xProfileZStressLocal( numXNodesGlobal ),
                    xProfileZStressGlobal( numXNodesGlobal ),
                    xProfileKineticEnergyLocal( numXNodesGlobal ),
                    xProfileKineticEnergyGlobal( numXNodesGlobal );

  for( int n=0; n < numXNodesGlobal; ++n )
  {
    xProfileDamageLocal[n] = 0.0;
    xProfileDamageGlobal[n] = 0.0;

    xProfileDensityLocal[n] = 0.0;
    xProfileDensityGlobal[n] = 0.0;

    xProfileXStressLocal[n] = 0.0;
    xProfileXStressGlobal[n] = 0.0;

    xProfileYStressLocal[n] = 0.0;
    xProfileYStressGlobal[n] = 0.0;

    xProfileZStressLocal[n] = 0.0;
    xProfileZStressGlobal[n] = 0.0;

    xProfileKineticEnergyLocal[n] = 0.0;
    xProfileKineticEnergyGlobal[n] = 0.0;
  }

  // array will store density at each x-nodal position, averaged over the domain.
  // Compute the damage, which is stored in the mass-weighted variable
  // g_damage = (sum m_p*D_p)/m_i for each field (which also assumes p_flag=2 equates to D=1)
  // dividing by the summed mass should give mass fractino of damaged material.

  // array will store density at each x-nodal position, averaged over the domain.
  // Compute density first by summing mass across all nodes and fields.
  // set profile to 0 before sum (may not be needed).

  // index in the global x-position array of the current node.
  int iGlobal,
      jGlobal,
      kGlobal,
      iMin,
      iMax,
      jMin,
      jMax,
      kMin,
      kMax;

  real64 smallX = 1.0e-8*(m_xLocalMax[0] - m_xLocalMin[0]);
  real64 smallY = 1.0e-8*(m_xLocalMax[1] - m_xLocalMin[1]);
  real64 smallZ = 1.0e-8*(m_xLocalMax[2] - m_xLocalMin[2]);

  // limits are defined so that a node is counted if imin <= i <=imax
  iMin = fabs(m_xLocalMin[0] - m_xGlobalMin[0]) < smallX ? 0 : round( ( m_xLocalMin[0] - m_xGlobalMin[0] ) / m_hEl[0] ) + 2;
  iMax = fabs(m_xLocalMax[0] - m_xGlobalMax[0]) < smallX ? numXNodesGlobal : round( ( m_xLocalMax[0] - m_xGlobalMin[0] ) / m_hEl[0] ) - 1;

  jMin = fabs(m_xLocalMin[1] - m_xGlobalMin[1]) < smallY ? 0 : round( ( m_xLocalMin[1] - m_xGlobalMin[1] ) / m_hEl[1] ) + 2;
  jMax = fabs(m_xLocalMax[1] - m_xGlobalMax[1]) < smallY ? numYNodesGlobal : round( ( m_xLocalMax[1] - m_xGlobalMin[1] ) / m_hEl[1] ) - 1;

  kMin = fabs(m_xLocalMin[2] - m_xGlobalMin[2]) < smallZ ? 0 : round( ( m_xLocalMin[2] - m_xGlobalMin[2] ) / m_hEl[2] ) + 2;
  kMax = fabs(m_xLocalMax[2] - m_xGlobalMax[2]) < smallZ ? numZNodesGlobal : round( ( m_xLocalMax[2] - m_xGlobalMin[2] ) / m_hEl[2] ) - 1;


  arrayView1d< int const > const periodic = partition.getPeriodic();
  // int numNodes = nodeManager.size();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 > const gridMass = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView2d< real64 > const gridMassWeightedDamage = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassWeightedDamageString() );
  arrayView3d< real64 > const gridNormalStress = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridNormalStressString() );
  arrayView3d< real64 > const gridVelocity = nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );

  // sum all values to local array.
  for( localIndex g = 0 ; g < numXNodesLocal; ++g )
  {
    iGlobal = round( ( gridPosition[g][0] - m_xGlobalMin[0] ) / m_hEl[0] ); // index in the global x-position array of the current node.
    jGlobal = round( ( gridPosition[g][1] - m_xGlobalMin[1] ) / m_hEl[1] ); // index in the global y-position array of the current node.
    kGlobal = round( ( gridPosition[g][2] - m_xGlobalMin[2] ) / m_hEl[2] ); // index in the global z-position array of the current node.



    for( localIndex fieldIndex = 0; fieldIndex < m_numVelocityFields; ++fieldIndex )
    {
      // Avoid double counting ghost nodes for mass, which is being taken from a synced field..
      if( iGlobal >= iMin && iGlobal <= iMax &&
          jGlobal >= jMin && jGlobal <= (jMax-periodic[1]) &&
          kGlobal >= kMin && kGlobal <= (kMax-periodic[2]) )
      {
        // CC: Moved this inside fieldIndex loop, NEED TO DOUBLE CHECK!
        // mass weighted damage is computed from master particles and never synchronized between patches
        // so a simple sum is adequate

        xProfileDamageLocal[iGlobal] += gridMassWeightedDamage[g][fieldIndex];

        // Normal stress components are mapped to nodes from master particles and never syncronized between patches,
        // so a simple sum should be adequate
        xProfileXStressLocal[iGlobal] += gridNormalStress[g][fieldIndex][0];
        xProfileYStressLocal[iGlobal] += gridNormalStress[g][fieldIndex][1];
        xProfileZStressLocal[iGlobal] += gridNormalStress[g][fieldIndex][2];

        // Kinetic energy is computed relative to some initial velocity ( v0=[vx0,0,0] ) whre vx0 is input.  This allows
        // for micro kinetic energy profiles to be computed for reverse ballistic simulations.
        xProfileKineticEnergyLocal[iGlobal] += 0.5 * gridMass[g][fieldIndex] * (
            ( gridVelocity[g][fieldIndex][0] - m_xProfileVx0 ) * ( gridVelocity[g][fieldIndex][0] - m_xProfileVx0 ) +
            gridVelocity[g][fieldIndex][1] * gridVelocity[g][fieldIndex][1] +
            gridVelocity[g][fieldIndex][2] * gridVelocity[g][fieldIndex][2] );

        xProfileDensityLocal[iGlobal] += gridMass[g][fieldIndex];
      }
    } // end loop over field index
  } // end loop over grid nodes

  // do additive sync on global mass-eighted damage array (not yet divided by mass)
  MpiWrapper::allReduce( xProfileDamageLocal,
                         xProfileDamageGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // do additive sync on global mass array (e.g. mass)
  MpiWrapper::allReduce( xProfileDensityLocal,
                         xProfileDensityGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // do additive sync on global normalStress arrays (e.g. summed stress*volume)
  MpiWrapper::allReduce( xProfileXStressLocal,
                         xProfileXStressGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( xProfileYStressLocal,
                         xProfileYStressGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( xProfileZStressLocal,
                         xProfileZStressGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // do additive sync on global kinetic energy array (summed 0.5*m*V*V)
  MpiWrapper::allReduce( xProfileKineticEnergyLocal,
                         xProfileKineticEnergyGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // Divide the summed mass-weighted damage by the summed nodal mass (the latter is still stored
  // in the not-yet-normalized-by-volume "density" array.  Divide element-wise by mass, but
  // avoid division by zero for 0-mass nodes.
  //
  // The resulting array is now the mass fraction of damaged material at each x-nodal position.
  for( int i = 0 ; i < numXNodesGlobal; ++i )
  {
    xProfileDamageGlobal[i] = ( xProfileDensityGlobal[i] > m_smallMass ) ? ( xProfileDamageGlobal[i] / xProfileDensityGlobal[i] ) : 0.0;
  }

  // Compute cross section of global domain and divide to convert summed mass into density profile.
  // This is global domain cross section in y-z plane, excluding ghost cells.
  real64 volume = m_hEl[0] * ( m_xGlobalMax[1] - m_xGlobalMin[1] - (2.0 + periodic[1]) * m_hEl[1] ) * ( m_xGlobalMax[2] - m_xGlobalMin[2] - (2.0 + periodic[2]) * m_hEl[2] );

  for( localIndex g = 0 ; g < numXNodesLocal; ++g )
  {
    xProfileDensityGlobal[g] *= 1.0 / volume;
    xProfileXStressGlobal[g] *= 1.0 / volume;
    xProfileYStressGlobal[g] *= 1.0 / volume;
    xProfileZStressGlobal[g] *= 1.0 / volume;
    xProfileKineticEnergyGlobal[g] *= 1.0 / volume;
  }

  // proc 0 writes the profile data.
  int rank = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if( rank == 0 )
  {
    // ---------------------------------------------
    // Mass Profile
    //
    std::ofstream file;
    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_density.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal vaerage for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileDensityGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // Damage Profile
    //
    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_damage.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal vaerage for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileDamageGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // xStress Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_xStress.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileXStressGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // yStress Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_yStress.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < numXNodesGlobal ; ++i )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileYStressGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // zStress Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_zStress.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.

    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < ( numXNodesGlobal ) ; ++i )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileZStressGlobal[i];
    }
    file << std::endl;
    file.close();

    // ---------------------------------------------
    // KINETICENERGY Profile
    //

    //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
    //file.exceptions(file.exceptions() | std::ios::failbit);
    file.open( "x_profile_kineticEnergy.csv", std::ios::out | std::ios::app );
    if( file.fail() )
      throw std::ios_base::failure( std::strerror( errno ) );
    //make sure write fails with exception if something is wrong
    file.exceptions( file.exceptions() | std::ios::failbit | std::ifstream::badbit );

    if( cycleNumber == 0 )
    {
      // Write the header to the file.
      // time, x_0, x_1, ..., x_n

      file << "time";
      for( int i = 0 ; i < numXNodesGlobal; ++i )
      {
        file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << m_xGlobalMin[0] + i * m_hEl[0];
      }
      file << std::endl;
    }

    // Write the current time and the nodal average for each x position.
    file << std::setprecision( std::numeric_limits<long double>::digits10 ) << time + dt;
    for( int i = 0 ; i < numXNodesGlobal; ++i )
    {
      file << std::setprecision( std::numeric_limits<long double>::digits10 ) << "," << xProfileKineticEnergyGlobal[i];
    }
    file << std::endl;
    file.close();
  }
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsMPM, string const &, dataRepository::Group * const )

}
