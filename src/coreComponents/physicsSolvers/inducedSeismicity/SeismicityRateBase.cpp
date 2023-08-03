/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SeismicityRateBase.cpp
 */

#include "SeismicityRateBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;

//START_SPHINX_INCLUDE_CONSTRUCTOR
SeismicityRateBase::SeismicityRateBase( const string & name,
                              Group * const parent ):
  SolverBase( name, parent ),
  m_stressSolver( nullptr )
  {
    this->registerWrapper( viewKeyStruct::stressSolverNameString(), &m_stressSolverName ).
          setInputFlag( InputFlags::REQUIRED ).
          setDescription( "Name of solver for computing stress" );
    this->registerWrapper( viewKeyStruct::faultNormalString(), &m_faultNormal ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Fault normal direction" );
    this->registerWrapper( viewKeyStruct::faultShearString(), &m_faultShear ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Fault shear direction" );
    this->registerWrapper( viewKeyStruct::initialSigmaString(), &m_initialSigma ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Initial normal stress" );
    this->registerWrapper( viewKeyStruct::initialTauString(), &m_initialTau ).
          setInputFlag( InputFlags::OPTIONAL ).
          setDescription( "Initial shear stress" );
  }
//END_SPHINX_INCLUDE_CONSTRUCTOR

//START_SPHINX_INCLUDE_REGISTERDATAONMESH
void SeismicityRateBase::registerDataOnMesh( Group & meshBodies )
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
      subRegion.registerField< geos::fields::inducedSeismicity::meanStress >( getName() ).
        reference().resizeDimension< 1 >( 6 );

      subRegion.registerField< geos::fields::inducedSeismicity::initialmeanNormalStress >( getName() );
      subRegion.registerField< geos::fields::inducedSeismicity::initialmeanShearStress >( getName() );

      subRegion.registerField< geos::fields::inducedSeismicity::meanNormalStress >( getName() );
      subRegion.registerField< geos::fields::inducedSeismicity::meanNormalStress_n >( getName() );
      subRegion.registerField< geos::fields::inducedSeismicity::meanShearStress >( getName() );
      subRegion.registerField< geos::fields::inducedSeismicity::meanShearStress_n >( getName() );

      subRegion.registerField< geos::fields::inducedSeismicity::seismicityRate >( getName() );
    } );
   } );
}

void SeismicityRateBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // 1. Validate various models against each other (must have same phases and components)
  // validateConstitutiveModels( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // Initialize initial stresses as specified by user (normal and shear need to be specified if only calling flow solver)
      if ( m_initialSigma < 0 )
      {
        arrayView1d< real64 > const tempSigIni = subRegion.getField< geos::fields::inducedSeismicity::initialmeanNormalStress >();
        tempSigIni.setValues< parallelHostPolicy >( m_initialSigma );
        arrayView1d< real64 > const tempSig = subRegion.getField< geos::fields::inducedSeismicity::meanNormalStress >();
        tempSig.setValues< parallelHostPolicy >( m_initialSigma );
        arrayView1d< real64 > const tempSig_n = subRegion.getField< geos::fields::inducedSeismicity::meanNormalStress_n >();
        tempSig_n.setValues< parallelHostPolicy >( m_initialSigma );
      }

      if ( m_initialTau < 0 )
      {
        arrayView1d< real64 > const tempTauIni = subRegion.getField< geos::fields::inducedSeismicity::initialmeanShearStress >();
        tempTauIni.setValues< parallelHostPolicy >( m_initialTau );
        arrayView1d< real64 > const tempTau = subRegion.getField< geos::fields::inducedSeismicity::meanShearStress >();
        tempTau.setValues< parallelHostPolicy >( m_initialTau );
        arrayView1d< real64 > const tempTau_n = subRegion.getField< geos::fields::inducedSeismicity::meanShearStress_n >();
        tempTau_n.setValues< parallelHostPolicy >( m_initialTau );
      }
                 
      arrayView1d< real64 > const tempR = subRegion.getField< geos::fields::inducedSeismicity::seismicityRate >();
      tempR.setValues< parallelHostPolicy >( 1.0 );
    } );
  } );
}

void SeismicityRateBase::postProcessInput()
{
  initializeFaultOrientation();
  m_stressSolver = &this->getParent().getGroup< SolverBase >( m_stressSolverName );
  SolverBase::postProcessInput();
}

void SeismicityRateBase::initializeFaultOrientation()
{
  m_faultNormalVoigt[0] = m_faultNormal[0]*m_faultNormal[0];
  m_faultNormalVoigt[1] = m_faultNormal[1]*m_faultNormal[1];
  m_faultNormalVoigt[2] = m_faultNormal[2]*m_faultNormal[2];
  m_faultNormalVoigt[3] = 2*m_faultNormal[1]*m_faultNormal[2];
  m_faultNormalVoigt[4] = 2*m_faultNormal[0]*m_faultNormal[2];
  m_faultNormalVoigt[5] = 2*m_faultNormal[0]*m_faultNormal[1];

  m_faultShearVoigt[0] = m_faultShear[0]*m_faultNormal[0];
  m_faultShearVoigt[1] = m_faultShear[1]*m_faultNormal[1];
  m_faultShearVoigt[2] = m_faultShear[2]*m_faultNormal[2];
  m_faultShearVoigt[3] = m_faultShear[1]*m_faultNormal[2] + m_faultShear[2]*m_faultNormal[1];
  m_faultShearVoigt[4] = m_faultShear[0]*m_faultNormal[2] + m_faultShear[2]*m_faultNormal[0];
  m_faultShearVoigt[5] = m_faultShear[0]*m_faultNormal[1] + m_faultShear[1]*m_faultNormal[0];
}

SeismicityRateBase::~SeismicityRateBase()
{
  // TODO Auto-generated destructor stub
}

} // namespace geos
