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
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

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
      subRegion.registerField< inducedSeismicity::meanStress >( getName() ).
        reference().resizeDimension< 1 >( 6 );

      subRegion.registerField< inducedSeismicity::initialMeanNormalStress >( getName() );
      subRegion.registerField< inducedSeismicity::initialMeanShearStress >( getName() );

      subRegion.registerField< inducedSeismicity::meanNormalStress >( getName() );
      subRegion.registerField< inducedSeismicity::meanNormalStress_n >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStress >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStress_n >( getName() );

      subRegion.registerField< inducedSeismicity::seismicityRate >( getName() );
    } );
   } );
}

void SeismicityRateBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const tempR = subRegion.getField< inducedSeismicity::seismicityRate >();
      tempR.setValues< parallelHostPolicy >( 1.0 );
    } );
  } );
}

void SeismicityRateBase::updateMeanSolidStress( ElementSubRegionBase & subRegion )
{
  arrayView2d< real64 > const meanStress = subRegion.getField< inducedSeismicity::meanStress >();
  meanStress.setValues< parallelHostPolicy >( 0.0 );

  arrayView1d< real64 > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
  arrayView1d< real64 > const sig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
  
  arrayView1d< real64 > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
  arrayView1d< real64 > const tau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();

  // get total stresses
  string const & solidModelName = subRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
  SolidBase const & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidModelName );
  arrayView3d< real64 const, solid::STRESS_USD > const stress = solidModel.getStress();
  
  // loop over elements 
  forAll< parallelDevicePolicy<> >(  sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    // update previous solutions
    sig_n[k] = sig[k];
    tau_n[k] = tau[k];

    // update mean stress at element
    for( int q = 0; q < stress.size(1); q++ ) 
    {
      LvArray::tensorOps::add< 6 >( meanStress[k], stress[k][q] );
    }

    // average
    LvArray::tensorOps::scale< 6 >(meanStress[k], 1./stress.size(1));

    // Compute fault normal and shear stresses
    sig[k] = LvArray::tensorOps::AiBi< 6 >( meanStress[k], m_faultNormalVoigt);
    tau[k] = LvArray::tensorOps::AiBi< 6 >( meanStress[k], m_faultShearVoigt);
  } );
}

void SeismicityRateBase::initializeMeanSolidStress(integer const cycleNumber, DomainPartition & domain)
{
  if ( cycleNumber == 0 ) {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )

    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
        arrayView1d< real64 > const sig_i = subRegion.getField< inducedSeismicity::initialMeanNormalStress >();
  
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
        arrayView1d< real64 > const tau_i = subRegion.getField< inducedSeismicity::initialMeanShearStress >();

        // update mean solid stress at each element
        if ( subRegion.hasWrapper( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString() ) )
        {
          initializeFaultOrientation();

          updateMeanSolidStress( subRegion );
          forAll< parallelDevicePolicy<> >(  sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            // update previous solutions
            sig_i[k] = sig[k];
            tau_i[k] = tau[k];
          } );

        }
        else
        {
          arrayView1d< real64 > const tempSigIni = subRegion.getField< inducedSeismicity::initialMeanNormalStress >();
          tempSigIni.setValues< parallelHostPolicy >( m_initialSigma );
          arrayView1d< real64 > const tempSig = subRegion.getField< inducedSeismicity::meanNormalStress >();
          tempSig.setValues< parallelHostPolicy >( m_initialSigma );
          arrayView1d< real64 > const tempSig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
          tempSig_n.setValues< parallelHostPolicy >( m_initialSigma );

          arrayView1d< real64 > const tempTauIni = subRegion.getField< inducedSeismicity::initialMeanShearStress >();
          tempTauIni.setValues< parallelHostPolicy >( m_initialTau );
          arrayView1d< real64 > const tempTau = subRegion.getField< inducedSeismicity::meanShearStress >();
          tempTau.setValues< parallelHostPolicy >( m_initialTau );
          arrayView1d< real64 > const tempTau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();
          tempTau_n.setValues< parallelHostPolicy >( m_initialTau );
        }
      } );
    } );
  }
}

void SeismicityRateBase::postProcessInput()
{
  m_stressSolver = &this->getParent().getGroup< SolverBase >( m_stressSolverName );
  SolverBase::postProcessInput();
}

void SeismicityRateBase::initializeFaultOrientation()
{
  if ( m_faultNormal.size()==1 )
  {
    // THROW ERROR
  }
  else
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
}

SeismicityRateBase::~SeismicityRateBase()
{
  // TODO Auto-generated destructor stub
}

} // namespace geos
