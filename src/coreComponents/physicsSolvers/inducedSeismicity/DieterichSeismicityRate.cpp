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

namespace geos
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;
using namespace fields;

/*----------------------------------------------------------------------------------
 * LaplaceFEM: Solving Laplace's partial differential equation with finite elements
 * ---------------------------------------------------------------------------------
 *
 * What does this solver do?
 * --------------------------
 *
 * This solver finds a solution f(x,y,z) to the Laplace equation: div ( grad ( f )) = 0
 * This common elliptic PDE represents the solution of a steady-state heat transfer, for instance.
 *
 * Where can I find an example of what it does?
 * --------------------------------------------
 *
 * Integrated tests associated to this solver are found in the ./integratedTests/ folder
 * These tests consist of computing the steady-state temperature profile in a simple cube-shaped domain
 * with fixed temperatures applied on two opposite cube faces ("Dirichlet" boundary conditions: imposing a value).
 * Feel free to run these tests cases, check out the XML input files, and inspect the output.
 *
 * Implementation: before we start:
 * ---------------------------------
 * In this implementation, the solution function (called above f) is called m_fieldName.
 * The variable m_fieldName is a string that points to a data container (an array) that
 * holds the numerical values of the PDE solution for each location at which f is evaluated.
 *
 * Let's take a look at the implementation step by step.
 *
 * ---------------------------------------------------------------------------------
 */


/* CONSTRUCTOR
   First, let us inspect the constructor of a "LaplaceFEM" object.
   This constructor does three important things:
   1 - It constructs an instance of the LaplaceFEM class (here: using the SolverBase constructor and passing through the arguments).
   2 - It sets some default values for the LaplaceFEM-specific private variables (here: m_fieldName and m_timeIntegrationOption).
   3 - It creates and activates a "registerWrapper" for each private variable.
   This is where the private variables are declared either as REQUIRED or OPTIONAL.
   An error is thrown if a REQUIRED variable is not specified in the XML file,
   along with the description of this variable and possible enum values if relevant.
   The description that is set is used in auto-generated documentation and console error messages.
 */

//START_SPHINX_INCLUDE_CONSTRUCTOR
DieterichSeismicityRate::DieterichSeismicityRate( const string & name,
                                                  Group * const parent ):
  SeismicityRateBase( name, parent ) {}
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
      subRegion.registerField< inducedSeismicity::t_a >( getName() );
      subRegion.registerField< inducedSeismicity::aSigma >( getName() );

      subRegion.registerField< inducedSeismicity::initialmeanNormalStress >( getName() );
      subRegion.registerField< inducedSeismicity::initialmeanShearStress >( getName() );

      subRegion.registerField< inducedSeismicity::pressureRate >( getName() );
      subRegion.registerField< inducedSeismicity::meanNormalStress >( getName() );
      subRegion.registerField< inducedSeismicity::meanNormalStress_n >( getName() );
      subRegion.registerField< inducedSeismicity::meanNormalStressRate >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStress >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStress_n >( getName() );
      subRegion.registerField< inducedSeismicity::meanShearStressRate >( getName() );

      subRegion.registerField< inducedSeismicity::seismicityRate >( getName() );
      subRegion.registerField< inducedSeismicity::h >( getName() );
      subRegion.registerField< inducedSeismicity::h_n >( getName() );
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
  // odeSolverStep( time_n, dt, cycleNumber, domain ); 
  real64 dtStress = m_stressSolver->solverStep(time_n, dt, cycleNumber, domain );
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
        arrayView1d< real64 > const h = subRegion.getField< inducedSeismicity::h >();
        arrayView1d< real64 > const h_n = subRegion.getField< inducedSeismicity::h_n >();
  
        arrayView1d< real64 const > const t_a = subRegion.getField< inducedSeismicity::t_a >();
        arrayView1d< real64 const > const aSig = subRegion.getField< inducedSeismicity::aSigma >();
        
        arrayView1d< real64 const > const p_i = subRegion.getField< fields::flow::initialPressure >();
        arrayView1d< real64 const > const p = subRegion.getField< fields::flow::pressure >();
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
            real64 curTau = aSig[k]*std::log(c*(time_n+dt)+1) + tau_i[k];
            real64 curTauDot = c*aSig[k]/(c*(time_n+dt)+1);

            // real64 gdot = ((tauDot[k]+aSig[k]/t_a[k])*(sig[k]-p[k]) + (tau[k]+aSig[k]/t_a[k]*(time_n+dt))*pDot[k]) / 
            //                   (aSig[k]/sig_i[k]*std::pow((sig[k]-p[k]), 2));
            real64 gdot = ((curTauDot+aSig[k]/t_a[k])*(sig[k]-p[k]) + (curTau+aSig[k]/t_a[k]*(time_n+dt))*pDot[k]) / 
                            (aSig[k]/sig_i[k]*std::pow((sig[k]-p[k]), 2));

            real64 f = h[k] + dt/t_a[k]*LvArray::math::exp(h[k]) - (h_n[k]  + dt*gdot);
            real64 dfdh = 1 + dt/t_a[k]*LvArray::math::exp(h[k]);
  
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
        arrayView1d< real64 > const h = subRegion.getField< inducedSeismicity::h >();
        arrayView1d< real64 > const logDenom = subRegion.getField< inducedSeismicity::logDenom >();
        arrayView1d< real64 > const logDenom_n = subRegion.getField< inducedSeismicity::logDenom_n >();

        arrayView1d< real64 const > const t_a = subRegion.getField< inducedSeismicity::t_a >();
        arrayView1d< real64 const > const aSig = subRegion.getField< inducedSeismicity::aSigma >();
        
        arrayView1d< real64 const > const p_i = subRegion.getField< fields::flow::initialPressure >();
        arrayView1d< real64 const > const p = subRegion.getField< fields::flow::pressure >();
        arrayView1d< real64 const > const p_n = subRegion.getField< fields::flow::pressure_n >();
        
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
          real64 curTau = aSig[k]*std::log(c*(time_n+dt)+1) + tau_i[k];
          real64 curTau_n = aSig[k]*std::log(c*time_n+1) + tau_i[k];

          real64 g = (tau[k] + aSig[k]/t_a[k]*(time_n+dt))/(aSig[k]/sig_i[k]*(sig[k]-p[k])) 
                              - tau_i[k]/(aSig[k]/sig_i[k]*(sig_i[k]-p_i[k]));
          real64 g_n = (tau_n[k] + aSig[k]/t_a[k]*time_n)/(aSig[k]/sig_i[k]*(sig_n[k]-p_n[k])) 
                              - tau_i[k]/(aSig[k]/sig_i[k]*(sig_i[k]-p_i[k]));

          real64 deltaLogDenom = std::log(1 + dt/(2*t_a[k])*(std::exp(g - logDenom_n[k]) + std::exp(g_n - logDenom_n[k]) ));
          logDenom[k] = logDenom_n[k] + deltaLogDenom;
  
          h[k] = g - logDenom[k];
          R[k] = LvArray::math::exp( h[k] );
          
          logDenom_n[k] = logDenom[k];
        } );
      } );
    } );
}

void DieterichSeismicityRate::initializePreSubGroups()
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
      // Hard coded stressing histories for now
      arrayView1d< real64 > const tempPDot = subRegion.getField< inducedSeismicity::pressureRate >();
      tempPDot.setValues< parallelHostPolicy >( 0.0 );

      arrayView1d< real64 > const tempSigIni = subRegion.getField< inducedSeismicity::initialmeanNormalStress >();
      tempSigIni.setValues< parallelHostPolicy >( 100e6 );
      arrayView1d< real64 > const tempSig = subRegion.getField< inducedSeismicity::meanNormalStress >();
      tempSig.setValues< parallelHostPolicy >( 100e6 );
      arrayView1d< real64 > const tempSig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
      tempSig_n.setValues< parallelHostPolicy >( 100e6 );
      arrayView1d< real64 > const tempSigDot = subRegion.getField< inducedSeismicity::meanNormalStressRate >();
      tempSigDot.setValues< parallelHostPolicy >( 0.0 );
      
      arrayView1d< real64 > const tempTauIni = subRegion.getField< inducedSeismicity::initialmeanShearStress >();
      tempTauIni.setValues< parallelHostPolicy >( 60e6 );
      arrayView1d< real64 > const tempTau = subRegion.getField< inducedSeismicity::meanShearStress >();
      tempTau.setValues< parallelHostPolicy >( 60e6 );
      arrayView1d< real64 > const tempTau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();
      tempTau_n.setValues< parallelHostPolicy >( 60e6 );
      arrayView1d< real64 > const tempTauDot = subRegion.getField< inducedSeismicity::meanShearStressRate >();
      tempTauDot.setValues< parallelHostPolicy >( 0.0 );
      
      arrayView1d< real64 > const tempR = subRegion.getField< inducedSeismicity::seismicityRate >();
      tempR.setValues< parallelHostPolicy >( 1.0 );
      arrayView1d< real64 > const tempH = subRegion.getField< inducedSeismicity::h >();
      tempH.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempH_n = subRegion.getField< inducedSeismicity::h_n >();
      tempH_n.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempLogDenom = subRegion.getField< inducedSeismicity::logDenom >();
      tempLogDenom.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempLogDenom_n = subRegion.getField< inducedSeismicity::logDenom_n >();
      tempLogDenom_n.setValues< parallelHostPolicy >( 0.0 );

      arrayView1d< real64 > const tempTa = subRegion.getField< inducedSeismicity::t_a >();
      tempTa.setValues< parallelHostPolicy >( 0.01*100e6/3.171e-5 );
      arrayView1d< real64 > const tempASig = subRegion.getField< inducedSeismicity::aSigma >();
      tempASig.setValues< parallelHostPolicy >( 0.01*100e6 );

    } );
  } );
}

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, DieterichSeismicityRate, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geos */
