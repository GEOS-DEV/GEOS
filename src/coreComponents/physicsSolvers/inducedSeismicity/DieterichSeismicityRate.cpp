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
 * @file DieterichSeismicityRate.cpp
 */

// Source includes
#include "DieterichSeismicityRate.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

/*----------------------------------------------------------------------------------
 * DieterichSeismicityRate: Solving the ODE for seismicity rate from Dieterich, 1994
 * ---------------------------------------------------------------------------------
 *
 * What does this solver do?
 * --------------------------
 *
 * This solver finds a solution R(x, t) - the seismicity rate - to the ODE formulated 
 * by Dieterich, 1994 given a certain stressing history. The stressing history can consist
 * of mechanical stresses and pore pressure. The solver class includes a member variable 
 * pointing to the stress solver that is specified in the XML file. SolverStep for the 
 * stress solver is then called in the SolverStep function for the seismicity rate, to take
 * the updated stress history as the input.
 * 
 * Solving the ODE is currently implemented in two different ways: 1. By backward-Euler 
 * discretization of the ODE and 2. Computing the closed-form solution to the ODE which 
 * involves numerical calculation of an integral of a stress functional. Both methods solve
 * for the log of the seismicity rate in order to avoid overflow that typically accompanies
 * the exponential of the stress history in the original formulations. 
 * 
 * 
 * Where can I find an example of what it does?
 * --------------------------------------------
 *
 * TODO
 *
 * ---------------------------------------------------------------------------------
 */


/* CONSTRUCTOR */

//START_SPHINX_INCLUDE_CONSTRUCTOR
DieterichSeismicityRate::DieterichSeismicityRate( const string & name,
                                                  Group * const parent ):
  SeismicityRateBase( name, parent ) 
  {
    this->registerWrapper( viewKeyStruct::directEffectString(), &m_directEffect ).
          setInputFlag( InputFlags::REQUIRED ).
          setDescription( "Rate-and-state friction direct effect parameter" );
    this->registerWrapper( viewKeyStruct::backgroundStressingRateString(), &m_backgroundStressingRate ).
          setInputFlag( InputFlags::REQUIRED ).
          setDescription( "Background stressing rate" );
  }
//END_SPHINX_INCLUDE_CONSTRUCTOR

DieterichSeismicityRate::~DieterichSeismicityRate()
{
  // TODO Auto-generated destructor stub
}

//START_SPHINX_INCLUDE_REGISTERDATAONMESH
void DieterichSeismicityRate::registerDataOnMesh( Group & meshBodies )
{
  SeismicityRateBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< inducedSeismicity::directEffect >( getName() );
      subRegion.registerField< inducedSeismicity::backgroundStressingRate >( getName() );

      subRegion.registerField< inducedSeismicity::logSeismicityRate >( getName() );
      subRegion.registerField< inducedSeismicity::logSeismicityRate_n >( getName() );
      subRegion.registerField< inducedSeismicity::logDenom >( getName() );
      subRegion.registerField< inducedSeismicity::logDenom_n >( getName() );
    } );
   } );
}
//END_SPHINX_INCLUDE_REGISTERDATAONMESH

real64 DieterichSeismicityRate::solverStep( real64 const & time_n,
                                  real64 const & dt,
                                  const int cycleNumber,
                                  DomainPartition & domain )
{
  // Call stress solver
  real64 dtStress = m_stressSolver->solverStep(time_n, dt, cycleNumber, domain );

  // update mean solid stress at each element
  updateMeanSolidStress( domain );

  // Solve for seismicity rate
  integralSolverStep( time_n, dt, cycleNumber, domain ); 

  return dtStress;
}

// Solve backward-Euler discretization of ODE by Newton-Raphson
void DieterichSeismicityRate::odeSolverStep( real64 const & time_n,
                    real64 const & dt,
                    const int cycleNumber,
                    DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        arrayView1d< real64 > const R = subRegion.getField< inducedSeismicity::seismicityRate >();
        arrayView1d< real64 > const h = subRegion.getField< inducedSeismicity::logSeismicityRate >();
        arrayView1d< real64 > const h_n = subRegion.getField< inducedSeismicity::logSeismicityRate_n >();
  
        arrayView1d< real64 const > const directEffect = subRegion.getField< inducedSeismicity::directEffect >();
        arrayView1d< real64 const > const backgroundStressingRate = subRegion.getField< inducedSeismicity::backgroundStressingRate >();
        
        arrayView1d< real64 const > const p_i = subRegion.getField< flow::initialPressure >();
        arrayView1d< real64 const > const p = subRegion.getField< flow::pressure >();
        arrayView1d< real64 const > const pDot = subRegion.getField< inducedSeismicity::pressureRate >();
  
        arrayView1d< real64 const > const sig_i = subRegion.getField< inducedSeismicity::initialmeanNormalStress >();
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
        arrayView1d< real64 const > const sigDot = subRegion.getField< inducedSeismicity::meanNormalStressRate >();
        
        arrayView1d< real64 const > const tau_i = subRegion.getField< inducedSeismicity::initialmeanShearStress >();
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
        arrayView1d< real64 const > const tauDot = subRegion.getField< inducedSeismicity::meanShearStressRate >();
  
        // solve for logarithm of seismicity rate
        real64 tol = 1e-8;

        real64 c = 1e-3;
        forAll< parallelDevicePolicy<> >(  h.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          real64 error = 100;
          real64 h_prev = h[k];

          while (error > tol) 
          {
            // Hard code log meanShear stress history
            real64 curTau = directEffect[k]*sig_i[k]*std::log(c*(time_n+dt)+1) + tau_i[k];
            real64 curTauDot = c*directEffect[k]*sig_i[k]/(c*(time_n+dt)+1);

            // real64 gdot = ((tauDot[k]+aSig[k]/t_a[k])*(sig[k]-p[k]) + (tau[k]+aSig[k]/t_a[k]*(time_n+dt))*pDot[k]) / 
            //                   (aSig[k]/sig_i[k]*std::pow((sig[k]-p[k]), 2));
            real64 gdot = ((curTauDot+backgroundStressingRate[k])*(sig[k]-p[k]) + 
                                          (curTau+backgroundStressingRate[k]*(time_n+dt))*pDot[k]) / 
                            (directEffect[k]*std::pow((sig[k]-p[k]), 2));

            real64 f = h[k] + dt/(directEffect[k]*sig_i[k]/backgroundStressingRate[k])*LvArray::math::exp(h[k]) - (h_n[k]  + dt*gdot);
            real64 dfdh = 1 + dt/(directEffect[k]*sig_i[k]/backgroundStressingRate[k])*LvArray::math::exp(h[k]);
  
            h[k] = h_prev - f/dfdh;
  
            error = std::abs(h[k]-h_prev);
            h_prev=h[k];
          }
  
          h_n[k] = h[k];
          R[k]=LvArray::math::exp( h[k] );
        } );
      } );
    } );
}


// Solve integral solution to ODE
void DieterichSeismicityRate::integralSolverStep( real64 const & time_n,
                    real64 const & dt,
                    const int cycleNumber,
                    DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        arrayView1d< real64 > const R = subRegion.getField< inducedSeismicity::seismicityRate >();
        arrayView1d< real64 > const h = subRegion.getField< inducedSeismicity::logSeismicityRate >();
        arrayView1d< real64 > const logDenom = subRegion.getField< inducedSeismicity::logDenom >();
        arrayView1d< real64 > const logDenom_n = subRegion.getField< inducedSeismicity::logDenom_n >();

        arrayView1d< real64 const > const directEffect = subRegion.getField< inducedSeismicity::directEffect >();
        arrayView1d< real64 const > const backgroundStressingRate = subRegion.getField< inducedSeismicity::backgroundStressingRate >();
        
        arrayView1d< real64 const > const p_i = subRegion.getField< flow::initialPressure >();
        arrayView1d< real64 const > const p = subRegion.getField< flow::pressure >();
        arrayView1d< real64 const > const p_n = subRegion.getField< flow::pressure_n >();
        
        arrayView1d< real64 const > const sig_i = subRegion.getField< inducedSeismicity::initialmeanNormalStress >();
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
        arrayView1d< real64 const > const sig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
  
        arrayView1d< real64 const > const tau_i = subRegion.getField< inducedSeismicity::initialmeanShearStress >();
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
        arrayView1d< real64 const > const tau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();
                  
        // solve for logarithm of seismicity rate
        real64 c = 1e-3;

        forAll< parallelDevicePolicy<> >(  h.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          // Hard code log meanShear stress history
          // real64 curTau = directEffect[k]*sig_i[k]*std::log(c*(time_n+dt)+1) + tau_i[k];
          // real64 curTau_n = directEffect[k]*sig_i[k]*std::log(c*time_n+1) + tau_i[k];

          // real64 g = (curTau + backgroundStressingRate[k]*(time_n+dt))/(directEffect[k]*(sig[k]-p[k])) 
          //                     - tau_i[k]/(directEffect[k]*(sig_i[k]-p_i[k]));
          // real64 g_n = (curTau_n + backgroundStressingRate[k]*time_n)/(directEffect[k]*(sig_n[k]-p_n[k])) 
          //                     - tau_i[k]/(directEffect[k]*(sig_i[k]-p_i[k]));

          real64 g = (tau[k] + backgroundStressingRate[k]*(time_n+dt))/(directEffect[k]*(-sig[k]-p[k])) 
                              - tau_i[k]/(directEffect[k]*(-sig_i[k]-p_i[k]));
          real64 g_n = (tau_n[k] + backgroundStressingRate[k]*time_n)/(directEffect[k]*(-sig_n[k]-p_n[k])) 
                              - tau_i[k]/(directEffect[k]*(-sig_i[k]-p_i[k]));

          real64 deltaLogDenom = std::log(1 + dt/(2*(directEffect[k]*-sig_i[k]/backgroundStressingRate[k]))
                                              *(std::exp(g - logDenom_n[k]) + std::exp(g_n - logDenom_n[k]) ));
          logDenom[k] = logDenom_n[k] + deltaLogDenom;
  
          h[k] = g - logDenom[k];
          R[k] = LvArray::math::exp( h[k] );
          
          logDenom_n[k] = logDenom[k];
        } );
      } );
    } );
}

void DieterichSeismicityRate::updateMeanSolidStress( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        arrayView2d< real64 > const meanStress = subRegion.getField< inducedSeismicity::meanStress >();
        meanStress.setValues< parallelHostPolicy >( 0.0 );

        arrayView1d< real64 > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
        arrayView1d< real64 > const sig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
  
        arrayView1d< real64 > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
        arrayView1d< real64 > const tau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();

        // get solid stresses
        // string const & solidModelName = subRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
        // SolidBase const & solidModel = getConstitutiveModel< SolidBase >( subRegion, solidModelName );
        // arrayView3d< real64 const, solid::STRESS_USD > stress = solidModel.getStress();
        
        string const porousSolidName = subRegion.template getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
        PorousSolidBase & porousSolid = subRegion.template getConstitutiveModel< PorousSolidBase >( porousSolidName );

        constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
        {
          typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousMaterialWrapper = castedPorousSolid.createKernelUpdates();
          // loop over elements 
          forAll< parallelDevicePolicy<> >(  sig.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
          {
            // update previous solutions
            sig_n[k] = sig[k];
            tau_n[k] = tau[k];

            // loop over quadrature points and sum
            // use a normal loop

            for(integer q = 0; q < numQuadPoints; q++ ) 
            {
              LvArray::tensorOps::add< 6 >( meanStress[k], stress[k][q] );
            }

            // average
            LvArray::tensorOps::scale< 6 >(meanStress[k], 1./stress.size(1));

            // Compute fault normal and shear stresses
            sig[k] = LvArray::tensorOps::AiBi< 6 >( meanStress[k], m_faultNormalVoigt);
            tau[k] = LvArray::tensorOps::AiBi< 6 >( meanStress[k], m_faultShearVoigt);
          } );
        } );
      } );
    } );
}


void DieterichSeismicityRate::initializePreSubGroups()
{
  SeismicityRateBase::initializePreSubGroups();

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
      arrayView1d< real64 > const tempH = subRegion.getField< inducedSeismicity::logSeismicityRate >();
      tempH.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempH_n = subRegion.getField< inducedSeismicity::logSeismicityRate_n >();
      tempH_n.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempLogDenom = subRegion.getField< inducedSeismicity::logDenom >();
      tempLogDenom.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempLogDenom_n = subRegion.getField< inducedSeismicity::logDenom_n >();
      tempLogDenom_n.setValues< parallelHostPolicy >( 0.0 );

      arrayView1d< real64 > const tempA = subRegion.getField< inducedSeismicity::directEffect >();
      tempA.setValues< parallelHostPolicy >( m_directEffect );
      arrayView1d< real64 > const tempTaur = subRegion.getField< inducedSeismicity::backgroundStressingRate >();
      tempTaur.setValues< parallelHostPolicy >( m_backgroundStressingRate );
    } );
  } );
}

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, DieterichSeismicityRate, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geos */
