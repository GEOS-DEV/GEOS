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
 * Solving the ODE is currently implemented by computing the closed-form interal solution 
 * to the ODE which involves numerical calculation of an integral of a stress functional. 
 * We initially solve for the log of the seismicity rate in order to avoid overflow that 
 * typically occurs in the exponential of the stress history. 
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
  // Save initial stress state on pre-defined fault orienations to field variables
  m_stressSolver->solverStep(time_n, 0.0, cycleNumber, domain );
  initializeMeanSolidStress(cycleNumber, domain); 

  // Call member variable stress solver to update the stress state
  real64 dtStress = m_stressSolver->solverStep(time_n, dt, cycleNumber, domain );

  // Loop over subRegions to update stress on faults and solver for seismicity rate
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // project new stress state to update stress on fault
      if ( subRegion.hasWrapper( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString() ) )
      {
        updateMeanSolidStress( subRegion );
      }

      // solve for the seismcity rate given new stresses on faults
      integralSolverStep( time_n, dtStress, subRegion ); 
    });
  });

  // return time step size achieved by stress solver
  return dtStress;
}

// Solve integral solution to ODE
void DieterichSeismicityRate::integralSolverStep( real64 const & time_n,
                    real64 const & dt,
                    ElementSubRegionBase & subRegion )
{
  // Retrieve field variables
  arrayView1d< real64 > const R = subRegion.getField< inducedSeismicity::seismicityRate >();
  arrayView1d< real64 > const logDenom = subRegion.getField< inducedSeismicity::logDenom >();
  arrayView1d< real64 > const logDenom_n = subRegion.getField< inducedSeismicity::logDenom_n >();

  arrayView1d< real64 const > const directEffect = subRegion.getField< inducedSeismicity::directEffect >();
  arrayView1d< real64 const > const backgroundStressingRate = subRegion.getField< inducedSeismicity::backgroundStressingRate >();
  
  arrayView1d< real64 const > const sig_i = subRegion.getField< inducedSeismicity::initialMeanNormalStress >();
  arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::meanNormalStress >();
  arrayView1d< real64 const > const sig_n = subRegion.getField< inducedSeismicity::meanNormalStress_n >();
  
  arrayView1d< real64 const > const tau_i = subRegion.getField< inducedSeismicity::initialMeanShearStress >();
  arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::meanShearStress >();
  arrayView1d< real64 const > const tau_n = subRegion.getField< inducedSeismicity::meanShearStress_n >();

  // Retrieve pore pressure variables if flow solver is being used. If not, initialize them as zero arrays.
  array1d< real64 > temp( R.size() );
  arrayView1d< real64 const > p_i = temp.toView();
  arrayView1d< real64 const > p = p_i;
  arrayView1d< real64 const > p_n = p_i;
  if ( subRegion.hasWrapper( FlowSolverBase::viewKeyStruct::fluidNamesString() ) )
  {
    p_i = subRegion.getField< flow::initialPressure >();
    p = subRegion.getField< flow::pressure >();
    p_n = subRegion.getField< flow::pressure_n >();
  }

  forAll< parallelDevicePolicy<> >( R.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    // arguments of stress exponential at current and previous time step
    real64 g = (tau[k] + backgroundStressingRate[k]*(time_n+dt))/(directEffect[k]*(-sig[k]-p[k])) 
                        - tau_i[k]/(directEffect[k]*(-sig_i[k]-p_i[k]));
    real64 g_n = (tau_n[k] + backgroundStressingRate[k]*time_n)/(directEffect[k]*(-sig_n[k]-p_n[k])) 
                        - tau_i[k]/(directEffect[k]*(-sig_i[k]-p_i[k]));

    // checkExpArgument();

    // Compute the difference of the log of the denominator of closed for integral solution.
    // This avoids directly computing the exponential of the current stress state which is more prone to overflow.
    logDenom[k] = logDenom_n[k] + std::log(1 + dt/(2*(directEffect[k]*-sig_i[k]/backgroundStressingRate[k]))
                                        *(std::exp(g - logDenom_n[k]) + std::exp(g_n - logDenom_n[k]) ));
  
    // Convert log seismicity rate to raw value
    R[k] = LvArray::math::exp( g - logDenom[k] );
    
    // Update log of the denominator for next time step
    logDenom_n[k] = logDenom[k];
  } );
}

void DieterichSeismicityRate::initializePreSubGroups()
{
  SeismicityRateBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
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

void checkExpArgument( real64 arg ){
  // TODO:
  // 1. CHECK IF CLOSE TO LITHOSTATIC PRESSURE
  // 2. CHECK IF STRESSING RATE IS TOO LARGE
  // 3. CHECK IF a IS TOO SMALL
}

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, DieterichSeismicityRate, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geos */
