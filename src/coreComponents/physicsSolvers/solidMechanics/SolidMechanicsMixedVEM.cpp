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
 * @file SolidMechanicsMixedVEM.cpp
 */

#include "physicsSolvers/solidMechanics/SolidMechanicsMixedVEM.hpp"

#include "constitutive/solid/SolidBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mixedVEM/HybridMixedVEM.hpp"
#include "mixedVEM/MixedVEMAssembly.hpp"
#include "mixedVEM/MixedVEMCellOutput.hpp"
#include "mixedVEM/MixedVEMDiscretization.hpp"
#include "mixedVEM/MixedVEMFields.hpp"
#include "mixedVEM/MixedVEMManager.hpp"

#include "common/FieldSpecificationOps.hpp"
#include "common/MpiWrapper.hpp"

#include <algorithm>
#include <numeric>
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace mixedVEM;

namespace
{

/// Face roles. A prescribed traction is essential here, a prescribed displacement is natural.
constexpr integer FACE_INTERIOR = 0;
constexpr integer FACE_BOUNDARY = 1;

/// All three displacement components prescribed.
constexpr integer FULL_DISPLACEMENT_MASK = 7;

/**
 * @brief Per subregion scratch, allocated once because a cell subregion has a fixed face count.
 */
struct ElementScratch
{
  explicit ElementScratch( integer const numFacesPerElement )
    : numFaces( numFacesPerElement ),
    numStressDof( NUM_FACE_DOF * numFacesPerElement ),
    faceGeom( static_cast< std::size_t >( numFacesPerElement ) )
  {
    divergence.resize( NUM_RM_DOF, numStressDof );
    divReconstruction.resize( NUM_RM_DOF, numStressDof );
    projection.resize( NUM_SYM_COMP, numStressDof );
    workspace.resize( NUM_SYM_COMP, numStressDof );
    stiffness.resize( numStressDof, numStressDof );
    factorization.resize( numStressDof, numStressDof );
    couplingTranspose.resize( NUM_RM_DOF, numStressDof );
    schur.resize( numStressDof, numStressDof );

    stressDofIndices.resize( static_cast< std::size_t >( numStressDof ) );
    multiplierDofIndices.resize( static_cast< std::size_t >( numStressDof ) );
    packedColumns.resize( static_cast< std::size_t >( numStressDof ) );
    packedValues.resize( static_cast< std::size_t >( numStressDof ) );
    stressRhs.resize( static_cast< std::size_t >( numStressDof ) );
    multiplier.resize( static_cast< std::size_t >( numStressDof ) );
    stress.resize( static_cast< std::size_t >( numStressDof ) );
  }

  integer numFaces;
  integer numStressDof;

  std::vector< FaceGeometry > faceGeom;
  ElementMoments moments;

  array2d< real64 > divergence;
  array2d< real64 > divReconstruction;
  array2d< real64 > projection;
  array2d< real64 > workspace;
  array2d< real64 > stiffness;

  array2d< real64 > factorization;
  array2d< real64 > couplingTranspose;
  array2d< real64 > schur;
  real64 inverseDivGram[NUM_RM_DOF][NUM_RM_DOF] = { { 0.0 } };

  std::vector< globalIndex > stressDofIndices;
  std::vector< globalIndex > multiplierDofIndices;
  std::vector< globalIndex > packedColumns;
  std::vector< real64 > packedValues;
  std::vector< real64 > stressRhs;
  std::vector< real64 > multiplier;
  std::vector< real64 > stress;
  real64 bodyForce[NUM_RM_DOF] = { 0.0 };
};

} // namespace

SolidMechanicsMixedVEM::SolidMechanicsMixedVEM( string const & name,
                                                Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_useHybridization( false )
{
  getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::REQUIRED );
}

void SolidMechanicsMixedVEM::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

  LinearSolverParameters & params = m_linearSolverParameters.get();
  params.isSymmetric = true;
}

void SolidMechanicsMixedVEM::initializePreSubGroups()
{
  PhysicsSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedVEMManager const & vemManager = numericalMethodManager.getMixedVEMManager();

  MixedVEMDiscretization const & discretization = vemManager.getDiscretization( m_discretizationName );

  m_useHybridization = discretization.useHybridization();

  // the interface operator is positive definite, the saddle point system is not
  LinearSolverParameters & params = m_linearSolverParameters.get();
  params.isSymmetric = true;

  if( m_useHybridization )
  {
    // H carries six unknowns per face, so multigrid is told the block size. The null space
    // is deliberately left at the constant modes: hypre's elasticity interpolation is built
    // for three unknowns per node, and on six it roughly triples the operator complexity and
    // loses the coarse grids under refinement. Its own per-unknown constants already carry a
    // global translation, because the three constant modes are written in the global frame.
    // amgNullSpaceType="rigidBodyModes" turns that interpolation back on. It wins on a
    // bending dominated problem at moderate size and loses the advantage under refinement.
    if( params.amg.numFunctions == 1 )
    {
      params.amg.numFunctions = NUM_FACE_DOF;
    }
  }
  else
  {
    params.direct.checkResidual = 0;

    // the saddle point system is reduced onto the rigid body motions by MGR
    if( params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr )
    {
      if( params.mgr.strategy == LinearSolverParameters::MGR::StrategyType::invalid )
      {
        params.mgr.strategy = LinearSolverParameters::MGR::StrategyType::solidMechanicsMixedVEM;
      }

      // the reduced system is an elasticity operator on the cells, so its coarse levels
      // need the rigid body motions
      if( params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::constantModes )
      {
        params.amg.nullSpaceType = LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes;
      }
    }
  }
}

void SolidMechanicsMixedVEM::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  PhysicsSolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( viewKeyStruct::solidMaterialNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
  solidMaterialName = PhysicsSolverBase::getConstitutiveName< SolidBase >( subRegion );

  GEOS_ERROR_IF( solidMaterialName.empty(),
                 GEOS_FMT( "{}: SolidBase model not found on subregion {}", getDataContext(), subRegion.getName() ) );
}

void SolidMechanicsMixedVEM::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    faceManager.registerField< fields::mixedVEM::faceStress >( getName() ).
      reference().resizeDimension< 1 >( NUM_FACE_DOF );

    faceManager.registerField< fields::mixedVEM::multiplier >( getName() ).
      reference().resizeDimension< 1 >( NUM_FACE_DOF );

    faceManager.registerField< fields::mixedVEM::displacementTrace >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    faceManager.registerField< fields::mixedVEM::traction >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    faceManager.registerField< fields::mixedVEM::boundaryType >( getName() );

    faceManager.registerField< fields::mixedVEM::displacementMask >( getName() );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::mixedVEM::rigidMotion >( getName() ).
        reference().resizeDimension< 1 >( NUM_RM_DOF );

      subRegion.registerField< fields::mixedVEM::displacement >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< fields::mixedVEM::rotation >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< fields::mixedVEM::stress >( getName() ).
        reference().resizeDimension< 1 >( NUM_SYM_COMP );
    } );
  } );
}

void SolidMechanicsMixedVEM::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                        DofManager & dofManager ) const
{
  if( m_useHybridization )
  {
    // the only global unknown is the displacement trace on the skeleton
    dofManager.addField( fields::mixedVEM::multiplier::key(),
                         FieldLocation::Face,
                         NUM_FACE_DOF,
                         getMeshTargets() );

    dofManager.addCoupling( fields::mixedVEM::multiplier::key(),
                            fields::mixedVEM::multiplier::key(),
                            DofManager::Connector::Elem );
  }
  else
  {
    dofManager.addField( fields::mixedVEM::faceStress::key(),
                         FieldLocation::Face,
                         NUM_FACE_DOF,
                         getMeshTargets() );

    dofManager.addCoupling( fields::mixedVEM::faceStress::key(),
                            fields::mixedVEM::faceStress::key(),
                            DofManager::Connector::Elem );

    dofManager.addField( fields::mixedVEM::rigidMotion::key(),
                         FieldLocation::Elem,
                         NUM_RM_DOF,
                         getMeshTargets() );

    dofManager.addCoupling( fields::mixedVEM::faceStress::key(),
                            fields::mixedVEM::rigidMotion::key(),
                            DofManager::Connector::Elem );
  }
}

void SolidMechanicsMixedVEM::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                real64 const & GEOS_UNUSED_PARAM( dt ),
                                                DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{}

void SolidMechanicsMixedVEM::computeNearNullSpace( DomainPartition & domain ) const
{
  LinearSolverParameters const & params = m_linearSolverParameters.get();

  if( params.amg.nullSpaceType != LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes ||
      !m_nearNullSpace.empty() )
  {
    return;
  }

  string const fieldKey = m_useHybridization
    ? fields::mixedVEM::multiplier::key()
    : fields::mixedVEM::rigidMotion::key();

  string const dofKey = m_dofManager.getKey( fieldKey );
  localIndex const numLocalDof = m_dofManager.numLocalDofs( fieldKey );
  // the offset of this field's block on this rank, not rankOffset, which counts only the
  // dofs of the field on previous ranks and would land past the block in a two field system
  globalIndex const dofOffset = m_dofManager.globalOffset( fieldKey );

  integer constexpr numModes = 2 * 3;

  m_nearNullSpace.resize( numModes );

  stdVector< arrayView1d< real64 > > values( numModes );
  for( integer k = 0; k < numModes; ++k )
  {
    m_nearNullSpace[k].create( numLocalDof, MPI_COMM_GEOS );
    values[k] = m_nearNullSpace[k].open();
    values[k].zero();
  }

  // the translations come first and the rotations last, the order hypre expects: an
  // unknown based coarsening already reproduces the constant modes, so only the last
  // three are handed over as interpolation vectors
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {
    if( m_useHybridization )
    {
      NodeManager const & nodeManager = mesh.getNodeManager();
      FaceManager const & faceManager = mesh.getFaceManager();

      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions = nodeManager.referencePosition();
      ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();
      arrayView2d< real64 const > const faceNormals = faceManager.faceNormal();
      arrayView1d< integer const > const ghostRank = faceManager.ghostRank();
      arrayView1d< globalIndex const > const dofNumber =
        faceManager.getReference< array1d< globalIndex > >( dofKey );

      real64 const origin[3] = { 0.0, 0.0, 0.0 };

      for( localIndex f = 0; f < faceManager.size(); ++f )
      {
        if( ghostRank[f] >= 0 )
        {
          continue;
        }

        localIndex const localDof = LvArray::integerConversion< localIndex >( dofNumber[f] - dofOffset );
        if( localDof < 0 || localDof + NUM_FACE_DOF > numLocalDof )
        {
          continue;
        }

        IndexedNodeCoordinates const X { nodePositions, faceToNodes[f] };
        integer const numNodes = LvArray::integerConversion< integer >( faceToNodes.sizeOfArray( f ) );

        real64 const normal[3] = { faceNormals( f, 0 ), faceNormals( f, 1 ), faceNormals( f, 2 ) };

        FaceGeometry geom;
        computeFaceGeometry( X, numNodes, normal, geom );

        real64 N[NUM_FACE_DOF][3][3];
        computeFaceBasisMoments( geom, origin, N );

        for( integer j = 0; j < NUM_FACE_DOF; ++j )
        {
          // translation e_k: lambda_j = e_k . int_f phi_j, which is one on its own mode
          if( j < 3 )
          {
            values[j][localDof + j] = 1.0;
          }

          // rotation e_k ^ x: lambda_j = e_k . int_f x ^ phi_j
          real64 b[3];
          faceBasisRotationalMoment( N[j], b );

          for( integer k = 0; k < 3; ++k )
          {
            values[3 + k][localDof + j] = b[k];
          }
        }
      }
    }
    else
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               CellElementSubRegion const & subRegion )
      {
        arrayView2d< real64 const > const elemCenters = subRegion.getElementCenter();
        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
        arrayView1d< globalIndex const > const dofNumber =
          subRegion.getReference< array1d< globalIndex > >( dofKey );

        for( localIndex k = 0; k < subRegion.size(); ++k )
        {
          if( ghostRank[k] >= 0 )
          {
            continue;
          }

          localIndex const localDof = LvArray::integerConversion< localIndex >( dofNumber[k] - dofOffset );
          if( localDof < 0 || localDof + NUM_RM_DOF > numLocalDof )
          {
            continue;
          }

          real64 const center[3] = { elemCenters( k, 0 ), elemCenters( k, 1 ), elemCenters( k, 2 ) };

          // u = a + omega ^ x restricted to E has coefficients (a + omega ^ x_E, omega)
          for( integer m = 0; m < 3; ++m )
          {
            values[m][localDof + m] = 1.0;

            real64 axis[3] = { 0.0, 0.0, 0.0 };
            axis[m] = 1.0;

            real64 rotation[3];
            LvArray::tensorOps::crossProduct( rotation, axis, center );

            for( integer i = 0; i < 3; ++i )
            {
              values[3 + m][localDof + i] = rotation[i];
            }
            values[3 + m][localDof + 3 + m] = 1.0;
          }
        }
      } );
    }
  } );

  for( integer k = 0; k < numModes; ++k )
  {
    m_nearNullSpace[k].close();
  }

  // orthonormalize so the rotations handed to multigrid are independent of the translations
  for( integer k = 0; k < numModes; ++k )
  {
    for( integer l = 0; l < k; ++l )
    {
      m_nearNullSpace[k].axpy( -m_nearNullSpace[k].dot( m_nearNullSpace[l] ), m_nearNullSpace[l] );
    }

    real64 const norm = m_nearNullSpace[k].norm2();
    GEOS_ERROR_IF( !( norm > 0.0 ),
                   GEOS_FMT( "{}: near null space mode {} is empty, so the degree of freedom "
                             "indexing of field {} is wrong", getDataContext(), k, fieldKey ) );

    m_nearNullSpace[k].scale( 1.0 / norm );
  }
}

std::unique_ptr< PreconditionerBase< LAInterface > >
SolidMechanicsMixedVEM::createPreconditioner( DomainPartition & domain ) const
{
  LinearSolverParameters const & params = m_linearSolverParameters.get();

  // the saddle point form only ever reaches multigrid through MGR; the hybridized form is
  // handed straight to it, and then only if the input asked for the rigid body modes
  bool const wantsNearNullSpace =
    params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes &&
    ( m_useHybridization
      ? params.preconditionerType == LinearSolverParameters::PreconditionerType::amg
      : params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr );

  if( wantsNearNullSpace )
  {
    computeNearNullSpace( domain );
    return LAInterface::createPreconditioner( params, m_nearNullSpace.toViewConst() );
  }

  return PhysicsSolverBase::createPreconditioner( domain );
}

real64 SolidMechanicsMixedVEM::solverStep( real64 const & time_n,
                                           real64 const & dt,
                                           integer const cycleNumber,
                                           DomainPartition & domain )
{
  // linearImplicitStep does not build the system itself, so do it on the first call and
  // again whenever the mesh changes
  Timestamp const meshTimestamp = getMeshModificationTimestamp( domain );
  if( meshTimestamp > getSystemSetupTimestamp() )
  {
    m_dofManager.clear();
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );
    setSystemSetupTimestamp( meshTimestamp );

    m_nearNullSpace.clear();
    m_precond = createPreconditioner( domain );
  }

  // the problem is linear, so one assembly and one solve are enough
  return linearImplicitStep( time_n, dt, cycleNumber, domain );
}

void SolidMechanicsMixedVEM::resetStateToBeginningOfStep( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{}

void SolidMechanicsMixedVEM::updateState( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // the unknowns are the primary fields, there is no secondary state to refresh
}

namespace
{

/**
 * @brief Remove the columns of the constrained degrees of freedom, keeping the matrix symmetric.
 * @param[in,out] localMatrix the local part of the system matrix
 * @param[in,out] localRhs the local residual
 * @param[in] rankOffset the first global row owned by this rank
 * @param[in,out] dofs the constrained degrees of freedom, in any order
 * @param[in,out] deltas the increment each of them is being driven to
 *
 * Zeroing the row of a constrained unknown fixes its value but leaves its column standing, so
 * the matrix stops being symmetric. That is fatal here: the whole point of the hybridized form
 * is a symmetric positive definite interface operator, and a conjugate gradient has no right
 * to work on anything else. Since the increment of a constrained unknown is known, its column
 * is moved to the residual and dropped, which restores symmetry.
 */
void eliminateConstrainedColumns( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs,
                                  globalIndex const rankOffset,
                                  stdVector< globalIndex > & dofs,
                                  stdVector< real64 > & deltas )
{
  if( dofs.empty() )
  {
    return;
  }

  stdVector< localIndex > order( dofs.size() );
  std::iota( order.begin(), order.end(), localIndex( 0 ) );
  std::sort( order.begin(), order.end(),
             [&]( localIndex const a, localIndex const b ){ return dofs[a] < dofs[b]; } );

  stdVector< globalIndex > sortedDof( dofs.size() );
  stdVector< real64 > sortedDelta( dofs.size() );
  for( std::size_t i = 0; i < order.size(); ++i )
  {
    sortedDof[i] = dofs[ static_cast< std::size_t >( order[i] ) ];
    sortedDelta[i] = deltas[ static_cast< std::size_t >( order[i] ) ];
  }

  for( localIndex r = 0; r < localMatrix.numRows(); ++r )
  {
    globalIndex const globalRow = r + rankOffset;

    arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( r );
    arraySlice1d< real64 > const entries = localMatrix.getEntries( r );

    for( localIndex q = 0; q < localMatrix.numNonZeros( r ); ++q )
    {
      if( columns[q] == globalRow || entries[q] == 0.0 )
      {
        continue;
      }

      auto const it = std::lower_bound( sortedDof.begin(), sortedDof.end(), columns[q] );
      if( it != sortedDof.end() && *it == columns[q] )
      {
        localRhs[r] += entries[q] * sortedDelta[ static_cast< std::size_t >( std::distance( sortedDof.begin(), it ) ) ];
        entries[q] = 0.0;
      }
    }
  }
}

/**
 * @brief Build the geometry, the moments and the element operators of one cell.
 */
void buildElementOperators( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePositions,
                            ArrayOfArraysView< localIndex const > const & faceToNodes,
                            arrayView2d< real64 const > const & faceNormals,
                            arraySlice1d< localIndex const > const & elemToFaces,
                            arraySlice1d< localIndex const > const & elemToNodes,
                            real64 const (&elemCenter)[3],
                            real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP],
                            ElementScratch & scratch )
{
  buildElementGeometry( nodePositions, faceToNodes, faceNormals, elemToFaces,
                        scratch.numFaces, elemCenter, scratch.faceGeom.data(), scratch.moments );

  IndexedNodeCoordinates const X { nodePositions, elemToNodes };
  real64 const diameter = computeElementDiameter( X, LvArray::integerConversion< integer >( elemToNodes.size() ) );

  computeElementOperators( scratch.faceGeom.data(),
                           scratch.numFaces,
                           elemCenter,
                           diameter,
                           scratch.moments,
                           compliance,
                           scratch.divergence.toSlice(),
                           scratch.divReconstruction.toSlice(),
                           scratch.projection.toSlice(),
                           scratch.workspace.toSlice(),
                           scratch.stiffness.toSlice() );
}

/**
 * @brief Load g_E and f_E of one cell: the Dirichlet moments and the gravity load.
 */
void buildElementLoads( arrayView2d< real64 const > const & displacementTrace,
                        arrayView2d< real64 const > const & tractionField,
                        arrayView1d< integer const > const & boundaryType,
                        arrayView1d< integer const > const & displacementMask,
                        arraySlice1d< localIndex const > const & elemToFaces,
                        real64 const density,
                        real64 const (&gravityVector)[3],
                        ElementScratch & scratch )
{
  std::fill( scratch.stressRhs.begin(), scratch.stressRhs.end(), 0.0 );

  for( integer lf = 0; lf < scratch.numFaces; ++lf )
  {
    localIndex const faceIndex = elemToFaces[lf];
    if( boundaryType[faceIndex] == FACE_INTERIOR )
    {
      continue;
    }

    FaceGeometry const & geom = scratch.faceGeom[ static_cast< std::size_t >( lf ) ];

    real64 const g[3] = { displacementTrace( faceIndex, 0 ),
                          displacementTrace( faceIndex, 1 ),
                          displacementTrace( faceIndex, 2 ) };
    real64 const t[3] = { tractionField( faceIndex, 0 ),
                          tractionField( faceIndex, 1 ),
                          tractionField( faceIndex, 2 ) };

    FaceConstraint constraint;
    computeFaceConstraint( geom, displacementMask[faceIndex], g, t, constraint );

    // only the free modes load the right hand side, the essential ones are fixed later
    std::size_t const offset = static_cast< std::size_t >( NUM_FACE_DOF * lf );

    for( integer i = 0; i < NUM_FACE_DOF; ++i )
    {
      if( !constraint.essential[i] )
      {
        scratch.stressRhs[offset + static_cast< std::size_t >( i )] = constraint.load[i];
      }
    }
  }

  // (f_E)_i = int_E rho g . R_i, whose rotational part follows the first moment
  real64 rotational[3];
  LvArray::tensorOps::crossProduct( rotational, scratch.moments.firstMoment, gravityVector );

  for( integer i = 0; i < 3; ++i )
  {
    scratch.bodyForce[i] = density * scratch.moments.volume * gravityVector[i];
    scratch.bodyForce[3 + i] = density * rotational[i];
  }
}

} // namespace

void SolidMechanicsMixedVEM::classifyFaces( real64 const time, MeshLevel & mesh ) const
{
  FaceManager & faceManager = mesh.getFaceManager();

  arrayView1d< integer const > const domainBoundary = faceManager.getDomainBoundaryIndicator();

  arrayView1d< integer > const boundaryType = faceManager.getField< fields::mixedVEM::boundaryType >();
  arrayView1d< integer > const displacementMask = faceManager.getField< fields::mixedVEM::displacementMask >();
  arrayView2d< real64 > const displacementTrace = faceManager.getField< fields::mixedVEM::displacementTrace >();
  arrayView2d< real64 > const traction = faceManager.getField< fields::mixedVEM::traction >();

  // a free surface is traction free, which is the essential condition of the mixed form
  forAll< serialPolicy >( faceManager.size(), [=]( localIndex const f )
  {
    // a boundary face carries a traction free condition until something says otherwise
    boundaryType[f] = ( domainBoundary[f] == 1 ) ? FACE_BOUNDARY : FACE_INTERIOR;
    displacementMask[f] = 0;

    for( integer i = 0; i < 3; ++i )
    {
      displacementTrace( f, i ) = 0.0;
      traction( f, i ) = 0.0;
    }
  } );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

  // only the component this specification names becomes a prescribed displacement, the rest
  // of the face keeps whatever traction it was given
  fsManager.applyFieldValue< serialPolicy >( time, mesh, fields::mixedVEM::displacementTrace::key(),
                                             [&]( FieldSpecification const & fs,
                                                  SortedArrayView< localIndex const > const & targetSet )
  {
    integer const component = fs.getComponent();
    integer const bit = ( component < 0 ) ? FULL_DISPLACEMENT_MASK : ( 1 << component );

    forAll< serialPolicy >( targetSet.size(), [=]( localIndex const i )
    {
      displacementMask[targetSet[i]] |= bit;
    } );
  } );

  // a prescribed traction only supplies data, the component is essential either way
  fsManager.applyFieldValue< serialPolicy >( time, mesh, fields::mixedVEM::traction::key(),
                                             [&]( FieldSpecification const &,
                                                  SortedArrayView< localIndex const > const & )
  {} );
}

void SolidMechanicsMixedVEM::assembleSystem( real64 const time,
                                             real64 const dt,
                                             DomainPartition & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  globalIndex const rankOffset = dofManager.rankOffset();

  real64 const gravity[3] = { gravityVector()[0], gravityVector()[1], gravityVector()[2] };

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    // boundary data is evaluated at the end of the step, as in the other solvers
    classifyFaces( time + dt, mesh );

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions = nodeManager.referencePosition();
    ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const > const faceNormals = faceManager.faceNormal();

    arrayView1d< integer const > const boundaryType = faceManager.getField< fields::mixedVEM::boundaryType >();
    arrayView1d< integer const > const displacementMask = faceManager.getField< fields::mixedVEM::displacementMask >();
    arrayView2d< real64 const > const displacementTrace = faceManager.getField< fields::mixedVEM::displacementTrace >();
    arrayView2d< real64 const > const tractionField = faceManager.getField< fields::mixedVEM::traction >();
    arrayView2d< real64 const > const faceStress = faceManager.getField< fields::mixedVEM::faceStress >();
    arrayView2d< real64 const > const multiplierField = faceManager.getField< fields::mixedVEM::multiplier >();

    string const faceDofKey = m_useHybridization
      ? dofManager.getKey( fields::mixedVEM::multiplier::key() )
      : dofManager.getKey( fields::mixedVEM::faceStress::key() );

    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );

    string const elemDofKey = m_useHybridization
      ? string()
      : dofManager.getKey( fields::mixedVEM::rigidMotion::key() );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion const & subRegion )
    {
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      SolidBase const & solid = getConstitutiveModel< SolidBase >( subRegion, solidName );

      arrayView1d< real64 const > const bulkModulus = solid.getBulkModulus();
      arrayView1d< real64 const > const shearModulus = solid.getShearModulus();
      arrayView2d< real64 const > const density = solid.getDensity();

      arrayView2d< localIndex const > const elemToFaces = subRegion.faceList().toViewConst();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = subRegion.nodeList().toViewConst();
      arrayView2d< real64 const > const elemCenters = subRegion.getElementCenter();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const rigidMotion = subRegion.getField< fields::mixedVEM::rigidMotion >();

      integer const numFacesPerElement = LvArray::integerConversion< integer >( elemToFaces.size( 1 ) );
      integer const numStressDof = NUM_FACE_DOF * numFacesPerElement;

      ElementScratch scratch( numFacesPerElement );

      arrayView1d< globalIndex const > elemDofNumber;
      if( !m_useHybridization )
      {
        elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      }

      // ghost elements are visited too: a face row owned here may have its second
      // neighbour on another rank, and that contribution would otherwise be lost. Rows
      // this rank does not own are filtered out by the scatter itself.
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        GEOS_UNUSED_VAR( ghostRank );

        real64 const elemCenter[3] = { elemCenters( k, 0 ), elemCenters( k, 1 ), elemCenters( k, 2 ) };

        real64 const mu = shearModulus[k];
        real64 const lambda = bulkModulus[k] - 2.0 * mu / 3.0;

        real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
        makeIsotropicCompliance( lambda, mu, compliance );

        buildElementOperators( nodePositions, faceToNodes, faceNormals,
                               elemToFaces[k], elemToNodes[k], elemCenter, compliance, scratch );

        buildElementLoads( displacementTrace, tractionField, boundaryType, displacementMask,
                           elemToFaces[k], density( k, 0 ), gravity, scratch );

        if( !m_useHybridization )
        {
          gatherStressDofIndices( elemToFaces[k], numFacesPerElement, faceDofNumber,
                                  scratch.stressDofIndices.data() );

          globalIndex dispDofIndices[NUM_RM_DOF];
          for( integer m = 0; m < NUM_RM_DOF; ++m )
          {
            dispDofIndices[m] = elemDofNumber[k] + m;
          }

          addElementToMatrix< serialAtomic >( localMatrix,
                                              rankOffset,
                                              scratch.stressDofIndices.data(),
                                              numStressDof,
                                              dispDofIndices,
                                              scratch.stiffness.toSliceConst(),
                                              scratch.divergence.toSliceConst(),
                                              scratch.packedValues.data() );

          // residual R = A x - b, evaluated at the current state
          for( integer i = 0; i < numStressDof; ++i )
          {
            localIndex const localRow =
              LvArray::integerConversion< localIndex >( scratch.stressDofIndices[ static_cast< std::size_t >( i ) ] - rankOffset );
            if( localRow < 0 || localRow >= localRhs.size() )
            {
              continue;
            }

            real64 value = -scratch.stressRhs[ static_cast< std::size_t >( i ) ];

            for( integer j = 0; j < numStressDof; ++j )
            {
              localIndex const face = elemToFaces( k, j / NUM_FACE_DOF );
              value += scratch.stiffness( i, j ) * faceStress( face, j % NUM_FACE_DOF );
            }
            for( integer m = 0; m < NUM_RM_DOF; ++m )
            {
              value += scratch.divergence( m, i ) * rigidMotion( k, m );
            }

            RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow], value );
          }

          for( integer m = 0; m < NUM_RM_DOF; ++m )
          {
            localIndex const localRow =
              LvArray::integerConversion< localIndex >( dispDofIndices[m] - rankOffset );
            if( localRow < 0 || localRow >= localRhs.size() )
            {
              continue;
            }

            // the row is -B sigma - f, matching the sign used in the matrix
            real64 value = -scratch.bodyForce[m];
            for( integer j = 0; j < numStressDof; ++j )
            {
              localIndex const face = elemToFaces( k, j / NUM_FACE_DOF );
              value -= scratch.divergence( m, j ) * faceStress( face, j % NUM_FACE_DOF );
            }

            RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow], value );
          }
        }
        else
        {
          bool const condensed = computeLocalCondensation( scratch.stiffness.toSliceConst(),
                                                           scratch.divergence.toSliceConst(),
                                                           numFacesPerElement,
                                                           scratch.factorization.toSlice(),
                                                           scratch.couplingTranspose.toSlice(),
                                                           scratch.inverseDivGram,
                                                           scratch.schur.toSlice() );

          GEOS_ERROR_IF( !condensed,
                         GEOS_FMT( "{}: static condensation failed on element {} of {}",
                                   getDataContext(), k, subRegion.getName() ) );

          // sigma_E at lambda = 0, whose jump loads the interface problem
          std::fill( scratch.multiplier.begin(), scratch.multiplier.end(), 0.0 );

          real64 displacement[NUM_RM_DOF];
          recoverElementSolution( scratch.schur.toSliceConst(),
                                  scratch.couplingTranspose.toSliceConst(),
                                  scratch.inverseDivGram,
                                  scratch.faceGeom.data(),
                                  numFacesPerElement,
                                  scratch.multiplier.data(),
                                  scratch.stressRhs.data(),
                                  scratch.bodyForce,
                                  scratch.stress.data(),
                                  displacement );

          gatherStressDofIndices( elemToFaces[k], numFacesPerElement, faceDofNumber,
                                  scratch.multiplierDofIndices.data() );

          applyContinuityOperator( scratch.faceGeom.data(), numFacesPerElement, scratch.schur.toSlice() );

          addElementToInterfaceMatrix< serialAtomic >( localMatrix,
                                                       rankOffset,
                                                       scratch.multiplierDofIndices.data(),
                                                       numStressDof,
                                                       scratch.schur.toSliceConst(),
                                                       scratch.packedColumns.data(),
                                                       scratch.packedValues.data() );

          // R = H lambda + C_E sigma_E^0 - q, with q the prescribed traction of a free face
          for( integer i = 0; i < numStressDof; ++i )
          {
            localIndex const localRow =
              LvArray::integerConversion< localIndex >( scratch.multiplierDofIndices[ static_cast< std::size_t >( i ) ] - rankOffset );
            if( localRow < 0 || localRow >= localRhs.size() )
            {
              continue;
            }

            integer const lf = i / NUM_FACE_DOF;
            integer const mode = i % NUM_FACE_DOF;
            localIndex const face = elemToFaces( k, lf );

            FaceGeometry const & geom = scratch.faceGeom[ static_cast< std::size_t >( lf ) ];

            real64 value = geom.outwardSign * scratch.stress[ static_cast< std::size_t >( i ) ];

            for( integer j = 0; j < numStressDof; ++j )
            {
              localIndex const columnFace = elemToFaces( k, j / NUM_FACE_DOF );
              value += scratch.schur( i, j ) * multiplierField( columnFace, j % NUM_FACE_DOF );
            }

            if( boundaryType[face] != FACE_INTERIOR )
            {
              real64 const g[3] = { displacementTrace( face, 0 ),
                                    displacementTrace( face, 1 ),
                                    displacementTrace( face, 2 ) };
              real64 const t[3] = { tractionField( face, 0 ),
                                    tractionField( face, 1 ),
                                    tractionField( face, 2 ) };

              FaceConstraint constraint;
              computeFaceConstraint( geom, displacementMask[face], g, t, constraint );

              if( constraint.essential[mode] )
              {
                value -= geom.outwardSign * constraint.value[mode];
              }
            }

            RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow], value );
          }
        }
      }
    } );
  } );
}


void SolidMechanicsMixedVEM::applyBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time ),
                                                      real64 const GEOS_UNUSED_PARAM( dt ),
                                                      DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  globalIndex const rankOffset = dofManager.rankOffset();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions = nodeManager.referencePosition();
    ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const > const faceNormals = faceManager.faceNormal();

    arrayView2d< localIndex const > const faceToElementRegion = faceManager.elementRegionList();
    arrayView2d< localIndex const > const faceToElementSubRegion = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const faceToElement = faceManager.elementList();

    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< integer const > const boundaryType = faceManager.getField< fields::mixedVEM::boundaryType >();
    arrayView1d< integer const > const displacementMask = faceManager.getField< fields::mixedVEM::displacementMask >();
    arrayView2d< real64 const > const displacementTrace = faceManager.getField< fields::mixedVEM::displacementTrace >();
    arrayView2d< real64 const > const tractionField = faceManager.getField< fields::mixedVEM::traction >();
    arrayView2d< real64 const > const faceStress = faceManager.getField< fields::mixedVEM::faceStress >();
    arrayView2d< real64 const > const multiplierField = faceManager.getField< fields::mixedVEM::multiplier >();

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenters =
      elemManager.constructArrayViewAccessor< real64, 2 >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

    string const faceDofKey = m_useHybridization
      ? dofManager.getKey( fields::mixedVEM::multiplier::key() )
      : dofManager.getKey( fields::mixedVEM::faceStress::key() );

    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );

    stdVector< globalIndex > constrainedDof;
    stdVector< real64 > constrainedDelta;

    for( localIndex f = 0; f < faceManager.size(); ++f )
    {
      if( boundaryType[f] == FACE_INTERIOR )
      {
        continue;
      }

      integer const side = ( faceToElement( f, 0 ) >= 0 ) ? 0 : 1;
      localIndex const er = faceToElementRegion( f, side );
      localIndex const esr = faceToElementSubRegion( f, side );
      localIndex const ei = faceToElement( f, side );

      real64 const elemCenter[3] = { elemCenters[er][esr]( ei, 0 ),
                                     elemCenters[er][esr]( ei, 1 ),
                                     elemCenters[er][esr]( ei, 2 ) };

      IndexedNodeCoordinates const X { nodePositions, faceToNodes[f] };
      integer const numNodes = LvArray::integerConversion< integer >( faceToNodes.sizeOfArray( f ) );

      real64 const normal[3] = { faceNormals( f, 0 ), faceNormals( f, 1 ), faceNormals( f, 2 ) };

      FaceGeometry geom;
      computeFaceGeometry( X, numNodes, normal, geom );
      orientFace( elemCenter, geom );

      real64 const g[3] = { displacementTrace( f, 0 ), displacementTrace( f, 1 ), displacementTrace( f, 2 ) };
      real64 const t[3] = { tractionField( f, 0 ), tractionField( f, 1 ), tractionField( f, 2 ) };

      FaceConstraint constraint;
      computeFaceConstraint( geom, displacementMask[f], g, t, constraint );

      for( integer j = 0; j < NUM_FACE_DOF; ++j )
      {
        // the saddle point form fixes the stress of an essential mode; the hybridized form
        // instead freezes the multiplier of a free mode, whose data already sits in g_E
        bool const constrained = m_useHybridization ? !constraint.essential[j] : constraint.essential[j];
        if( !constrained )
        {
          continue;
        }

        globalIndex const dof = faceDofNumber[f] + j;
        localIndex const localRow = LvArray::integerConversion< localIndex >( dof - rankOffset );

        real64 const current = m_useHybridization ? multiplierField( f, j ) : faceStress( f, j );
        real64 const target = m_useHybridization ? 0.0 : constraint.value[j];

        // recorded whether or not this rank owns the row: a column of a locally owned row may
        // belong to a face this rank only ghosts, and it has to be eliminated as well
        constrainedDof.emplace_back( dof );
        constrainedDelta.emplace_back( target - current );

        if( faceGhostRank[f] >= 0 || localRow < 0 || localRow >= localRhs.size() )
        {
          continue;
        }

        real64 rhsValue = 0.0;
        FieldSpecificationEqual::SpecifyFieldValue( dof, rankOffset, localMatrix, rhsValue,
                                                    target, current );
        localRhs[localRow] = rhsValue;
      }
    }

    eliminateConstrainedColumns( localMatrix, localRhs, rankOffset,
                                 constrainedDof, constrainedDelta );

    // A constrained row is now an identity row, so a rigid body motion that is still nonzero
    // there is no longer a near null vector of the matrix multigrid is handed. The modes are
    // trimmed to the free unknowns instead. The preconditioner holds these vectors by pointer
    // and reads them when it sets up, which happens after this, so updating them in place is
    // enough.
    if( !m_nearNullSpace.empty() )
    {
      for( integer k = 0; k < m_nearNullSpace.size(); ++k )
      {
        arrayView1d< real64 > const values = m_nearNullSpace[k].open();

        for( globalIndex const dof : constrainedDof )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( dof - rankOffset );
          if( localRow >= 0 && localRow < values.size() )
          {
            values[localRow] = 0.0;
          }
        }
        m_nearNullSpace[k].close();
      }

      for( integer k = 0; k < m_nearNullSpace.size(); ++k )
      {
        real64 const norm = m_nearNullSpace[k].norm2();
        if( norm > 0.0 )
        {
          m_nearNullSpace[k].scale( 1.0 / norm );
        }
      }
    }

  } );
}

real64 SolidMechanicsMixedVEM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time ),
                                                      real64 const & GEOS_UNUSED_PARAM( dt ),
                                                      DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                      DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                      arrayView1d< real64 const > const & localRhs )
{
  real64 localSum = 0.0;
  for( localIndex i = 0; i < localRhs.size(); ++i )
  {
    localSum += localRhs[i] * localRhs[i];
  }

  real64 const residual = std::sqrt( MpiWrapper::sum( localSum, MPI_COMM_GEOS ) );

  GEOS_LOG_RANK_0_IF( getLogLevel() >= 1, GEOS_FMT( "        ( Rmixedvem ) = ( {:4.2e} )", residual ) );

  return residual;
}

void SolidMechanicsMixedVEM::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  real64 const GEOS_UNUSED_PARAM( dt ),
                                                  DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  DofManager::CompMask const mask( NUM_FACE_DOF, 0, NUM_FACE_DOF );

  FieldIdentifiers fieldsToBeSync;

  if( m_useHybridization )
  {
    dofManager.addVectorToField( localSolution,
                                 fields::mixedVEM::multiplier::key(),
                                 fields::mixedVEM::multiplier::key(),
                                 scalingFactor,
                                 mask );

    fieldsToBeSync.addFields( FieldLocation::Face, { fields::mixedVEM::multiplier::key() } );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 fields::mixedVEM::faceStress::key(),
                                 fields::mixedVEM::faceStress::key(),
                                 scalingFactor,
                                 mask );

    dofManager.addVectorToField( localSolution,
                                 fields::mixedVEM::rigidMotion::key(),
                                 fields::mixedVEM::rigidMotion::key(),
                                 scalingFactor,
                                 DofManager::CompMask( NUM_RM_DOF, 0, NUM_RM_DOF ) );

    fieldsToBeSync.addFields( FieldLocation::Face, { fields::mixedVEM::faceStress::key() } );
    fieldsToBeSync.addElementFields( { fields::mixedVEM::rigidMotion::key() }, getMeshTargets().begin()->second );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void SolidMechanicsMixedVEM::implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time ),
                                                   real64 const & GEOS_UNUSED_PARAM( dt ),
                                                   DomainPartition & domain )
{
  computeCellFields( domain );
}

void SolidMechanicsMixedVEM::computeCellFields( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  real64 const gravity[3] = { gravityVector()[0], gravityVector()[1], gravityVector()[2] };

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions = nodeManager.referencePosition();
    ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const > const faceNormals = faceManager.faceNormal();

    arrayView1d< integer const > const boundaryType = faceManager.getField< fields::mixedVEM::boundaryType >();
    arrayView1d< integer const > const displacementMask = faceManager.getField< fields::mixedVEM::displacementMask >();
    arrayView2d< real64 const > const displacementTrace = faceManager.getField< fields::mixedVEM::displacementTrace >();
    arrayView2d< real64 const > const tractionField = faceManager.getField< fields::mixedVEM::traction >();
    arrayView2d< real64 const > const multiplierField = faceManager.getField< fields::mixedVEM::multiplier >();
    arrayView2d< real64 > const faceStress = faceManager.getField< fields::mixedVEM::faceStress >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion & subRegion )
    {
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      SolidBase const & solid = getConstitutiveModel< SolidBase >( subRegion, solidName );

      arrayView1d< real64 const > const bulkModulus = solid.getBulkModulus();
      arrayView1d< real64 const > const shearModulus = solid.getShearModulus();
      arrayView2d< real64 const > const density = solid.getDensity();

      arrayView2d< localIndex const > const elemToFaces = subRegion.faceList().toViewConst();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = subRegion.nodeList().toViewConst();
      arrayView2d< real64 const > const elemCenters = subRegion.getElementCenter();

      arrayView2d< real64 > const rigidMotion = subRegion.getField< fields::mixedVEM::rigidMotion >();
      arrayView2d< real64 > const cellDisplacement = subRegion.getField< fields::mixedVEM::displacement >();
      arrayView2d< real64 > const cellRotation = subRegion.getField< fields::mixedVEM::rotation >();
      arrayView2d< real64 > const cellStress = subRegion.getField< fields::mixedVEM::stress >();

      integer const numFacesPerElement = LvArray::integerConversion< integer >( elemToFaces.size( 1 ) );
      integer const numStressDof = NUM_FACE_DOF * numFacesPerElement;

      ElementScratch scratch( numFacesPerElement );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        real64 const elemCenter[3] = { elemCenters( k, 0 ), elemCenters( k, 1 ), elemCenters( k, 2 ) };

        real64 const mu = shearModulus[k];
        real64 const lambda = bulkModulus[k] - 2.0 * mu / 3.0;

        real64 compliance[NUM_SYM_COMP][NUM_SYM_COMP];
        makeIsotropicCompliance( lambda, mu, compliance );

        buildElementOperators( nodePositions, faceToNodes, faceNormals,
                               elemToFaces[k], elemToNodes[k], elemCenter, compliance, scratch );

        if( m_useHybridization )
        {
          buildElementLoads( displacementTrace, tractionField, boundaryType, displacementMask,
                             elemToFaces[k], density( k, 0 ), gravity, scratch );

          bool const condensed = computeLocalCondensation( scratch.stiffness.toSliceConst(),
                                                           scratch.divergence.toSliceConst(),
                                                           numFacesPerElement,
                                                           scratch.factorization.toSlice(),
                                                           scratch.couplingTranspose.toSlice(),
                                                           scratch.inverseDivGram,
                                                           scratch.schur.toSlice() );
          GEOS_ERROR_IF( !condensed,
                         GEOS_FMT( "{}: static condensation failed during recovery", getDataContext() ) );

          for( integer i = 0; i < numStressDof; ++i )
          {
            localIndex const face = elemToFaces( k, i / NUM_FACE_DOF );
            scratch.multiplier[ static_cast< std::size_t >( i ) ] = multiplierField( face, i % NUM_FACE_DOF );
          }

          real64 displacement[NUM_RM_DOF];
          recoverElementSolution( scratch.schur.toSliceConst(),
                                  scratch.couplingTranspose.toSliceConst(),
                                  scratch.inverseDivGram,
                                  scratch.faceGeom.data(),
                                  numFacesPerElement,
                                  scratch.multiplier.data(),
                                  scratch.stressRhs.data(),
                                  scratch.bodyForce,
                                  scratch.stress.data(),
                                  displacement );

          for( integer m = 0; m < NUM_RM_DOF; ++m )
          {
            rigidMotion( k, m ) = displacement[m];
          }
          for( integer i = 0; i < numStressDof; ++i )
          {
            localIndex const face = elemToFaces( k, i / NUM_FACE_DOF );
            faceStress( face, i % NUM_FACE_DOF ) = scratch.stress[ static_cast< std::size_t >( i ) ];
          }
        }
        else
        {
          for( integer i = 0; i < numStressDof; ++i )
          {
            localIndex const face = elemToFaces( k, i / NUM_FACE_DOF );
            scratch.stress[ static_cast< std::size_t >( i ) ] = faceStress( face, i % NUM_FACE_DOF );
          }
        }

        real64 displacementValue[3], rotationValue[3];
        computeCellDisplacement( rigidMotion[k].dataIfContiguous(), displacementValue, rotationValue );

        real64 stressValue[NUM_SYM_COMP];
        computeCellStress( scratch.projection.toSliceConst(), numFacesPerElement,
                           scratch.stress.data(), stressValue );

        for( integer i = 0; i < 3; ++i )
        {
          cellDisplacement( k, i ) = displacementValue[i];
          cellRotation( k, i ) = rotationValue[i];
        }
        for( integer a = 0; a < NUM_SYM_COMP; ++a )
        {
          cellStress( k, a ) = stressValue[a];
        }
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsMixedVEM, string const &, Group * const )


} // namespace geos
