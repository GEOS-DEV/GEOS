/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file QuasiDynamicEQ.cpp
 */

#include "QuasiDynamicEQ.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "kernels/RateAndStateKernels.hpp"
#include "rateAndStateFields.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

QuasiDynamicEQ::QuasiDynamicEQ( const string & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_stressSolver( nullptr )
{
  this->registerWrapper( viewKeyStruct::maxNumberOfNewtonIterationsString(), &m_maxNewtonIterations ).
    setInputFlag( InputFlags::REQUIRED ).
    setApplyDefaultValue( 5 ).
    setDescription( "Maximum number of Newton iterations string." );

  this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of solver for computing stress" );  

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );  
}

void QuasiDynamicEQ::postInputInitialization()
{
  
  // Initialize member stress solver as specified in XML input
  if( !m_stressSolverName.empty() )
  {
    m_stressSolver = &this->getParent().getGroup< SolverBase >( m_stressSolverName );
  }

  SolverBase::postInputInitialization();
}

QuasiDynamicEQ::~QuasiDynamicEQ()
{
  // TODO Auto-generated destructor stub
}

void QuasiDynamicEQ::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< rateAndState::stateVariable >( getName() );
      subRegion.registerField< rateAndState::slipRate >( getName() );
      subRegion.registerField< rateAndState::stateVariable_n >( getName() );
      subRegion.registerField< rateAndState::slipRate_n >( getName() );
    } );
  } );
}

real64 QuasiDynamicEQ::solverStep( real64 const & time_n,
                                   real64 const & dt,
                                   const int cycleNumber,
                                   DomainPartition & domain )
{
  real64 const dtStress = updateStresses( time_n, dt, cycleNumber, domain );

  // Loop over subRegions to solve for seismicity rate
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     SurfaceElementSubRegion & subRegion )
    {
      // solve rate and state equations.
      rateAndStateKernels::createAndLaunch( subRegion, viewKeyStruct::frictionLawNameString(), m_maxNewtonIterations, time_n, dt );
      // save old state
      saveOldState( subRegion );
    } );
  } );

  // return time step size achieved by stress solver
  return dtStress;
}

real64 QuasiDynamicEQ::updateStresses( real64 const & time_n,
                                       real64 const & dt,
                                       const int cycleNumber,
                                       DomainPartition & domain ) const
{
  // Call member variable stress solver to update the stress state
  if( m_stressSolver )
  {

    // 1. Solve the momentum balance
    real64 const dtStress =  m_stressSolver->solverStep( time_n, dt, cycleNumber, domain );

    return dtStress;
  }
  else
  {
    // Spring-slider version
  }
  return dt;
}

void QuasiDynamicEQ::saveOldState( ElementSubRegionBase & subRegion ) const
{
  arrayView1d< real64 > const stateVariable   = subRegion.getField< rateAndState::stateVariable >();
  arrayView1d< real64 > const stateVariable_n = subRegion.getField< rateAndState::stateVariable_n >();
  arrayView1d< real64 > const slipRate        = subRegion.getField< rateAndState::slipRate >();
  arrayView1d< real64 > const slipRate_n      = subRegion.getField< rateAndState::slipRate_n >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    slipRate_n[k]      = slipRate[k];
    stateVariable_n[k] = stateVariable[k];
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, QuasiDynamicEQ, string const &, dataRepository::Group * const )
} // namespace geos
