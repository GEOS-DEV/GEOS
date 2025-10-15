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
 * @file SinglePhaseFVM.cpp
 */

#include "SinglePhaseStressPathDrivenFVM.hpp"

#include "physicsSolvers/LogLevelsInfo.hpp"
#include "constitutive/contact/HydraulicApertureRelationSelector.hpp"
#include "mesh/DomainPartition.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/solidMechanics/contact/FractureState.hpp"

/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace constitutive;

SinglePhaseStressPathDrivenFVM::SinglePhaseStressPathDrivenFVM( const string & name,
                                                                  dataRepository::Group * const parent ):
  SinglePhaseFVM< SinglePhaseBase >( name, parent ),
  m_oedometricStressPath( groupKeyStruct::oedometricStressPathString(), this )
{
  SinglePhaseFVM< SinglePhaseBase >::template addLogLevel< logInfo::Convergence >();

  registerGroup( groupKeyStruct::oedometricStressPathString(), &m_oedometricStressPath );
}

void SinglePhaseStressPathDrivenFVM::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::setConstitutiveNamesCallSuper( subRegion );
  
  // To make Barton Bandis constitutive law mandatory
  if( dynamic_cast< SurfaceElementSubRegion * >( &subRegion ) )
  {
    this->template setConstitutiveName< constitutive::BartonBandis >( subRegion,
                                                                      viewKeyStruct::hydraulicApertureRelationNameString(), "hydraulic aperture" );
  }
}

void SinglePhaseStressPathDrivenFVM::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  // As in PoromechanicsConformingFractures.hpp

  // remove the contribution of the hydraulic aperture from the stencil weights
  prepareStencilWeights( domain );

  updateHydraulicAperture( domain );

  // update the stencil weights using the updated hydraulic aperture
  updateStencilWeights( domain );

  SinglePhaseBase::updateState( domain ); // updates the permeability of the cell and the fracture 
}

void SinglePhaseStressPathDrivenFVM::updateHydraulicAperture( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames, [&]( localIndex const, 
                                                                                             auto & subRegion )
    {
      updateFractureAperture( subRegion );
    } );
  } );
}

void SinglePhaseStressPathDrivenFVM::updateFractureAperture( SurfaceElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;
  
  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 > const & newHydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
  arrayView2d< real64 const > const & normalVector = subRegion.getField< fields::normalVector >(); // mesh/MeshFields.hpp
  
  string const & hydraulicApertureRelationName = 
    subRegion.template getReference< string >( viewKeyStruct::hydraulicApertureRelationNameString()  );
  BartonBandis const & hydraulicApertureModel = 
    this->template getConstitutiveModel< BartonBandis >( subRegion, hydraulicApertureRelationName );

  BartonBandisUpdates hydraulicApertureWrapper = hydraulicApertureModel.createKernelWrapper(); 

  real64 sumAperture = 0.0;
  // not used
  real64 dHydraulicAperture_aperture = 0.0;
  real64 dHydraulicAperture_dNormalTraction = 0.0;
  forAll< parallelDevicePolicy<> >( subRegion.size(), [&] GEOS_DEVICE ( localIndex const k )
  {
    R1Tensor const normal = { normalVector[k][0], normalVector[k][1], normalVector[k][2] };
    
    real64 const sigmaN_N = m_oedometricStressPath.computeFractureStress( pressure[k], normal);
    // The traction sigmaN_N must be negative
    newHydraulicAperture[k] = hydraulicApertureWrapper.computeHydraulicAperture( 0.0,
                                                    -sigmaN_N,
                                                    fields::contact::FractureState::Slip, // not open
                                                    dHydraulicAperture_aperture,
                                                    dHydraulicAperture_dNormalTraction );

    sumAperture += newHydraulicAperture[k];
  } );

  real64 const averageAperture = sumAperture / subRegion.size();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [&newHydraulicAperture, averageAperture] GEOS_DEVICE ( localIndex const k )
  {
    newHydraulicAperture[k] = averageAperture;
  } );
}


REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseStressPathDrivenFVM, string const &, dataRepository::Group * const )

} /* namespace geos */
