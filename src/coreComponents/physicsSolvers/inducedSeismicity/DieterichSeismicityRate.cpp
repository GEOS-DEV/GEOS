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
#include "DieterichSeismicityRateKernels.hpp"

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
  SeismicityRateBase( name, parent ), 
  m_directEffect(0.01),
  m_bStressRate(3.171e-5),
  m_initialSigma(100e6),
  m_initialTau(60e6)
  {
  this->registerWrapper( viewKeyStruct::directEffect(), &m_directEffect ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Rate-and-state direct effect parameter" );
  this->registerWrapper( viewKeyStruct::bStressRate(), &m_bStressRate ).
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
      subRegion.registerField< inducedSeismicity::t_a >( getName() );
      subRegion.registerField< inducedSeismicity::aSigma >( getName() );

      subRegion.registerField< inducedSeismicity::pressure >( getName() );
      subRegion.registerField< inducedSeismicity::pressure_n >( getName() );
      subRegion.registerField< inducedSeismicity::pressureRate >( getName() );
      subRegion.registerField< inducedSeismicity::normalStress >( getName() );
      subRegion.registerField< inducedSeismicity::normalStress_n >( getName() );
      subRegion.registerField< inducedSeismicity::normalStressRate >( getName() );
      subRegion.registerField< inducedSeismicity::shearStress >( getName() );
      subRegion.registerField< inducedSeismicity::shearStress_n >( getName() );
      subRegion.registerField< inducedSeismicity::shearStressRate >( getName() );

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
  // Solve backward-Euler discretization of ODE by Newton-Raphson
  // odeSolverStep( time_n, dt, cycleNumber, domain ); 
  integralSolverStep( time_n, dt, cycleNumber, domain ); 

  // return this->linearImplicitStep( time_n, dt, cycleNumber, domain );
  return dt;
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
      
        arrayView1d< real64 const > const p = subRegion.getField< inducedSeismicity::pressure >();
        arrayView1d< real64 const > const pDot = subRegion.getField< inducedSeismicity::pressureRate >();
  
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::normalStress >();
        arrayView1d< real64 const > const sigDot = subRegion.getField< inducedSeismicity::normalStressRate >();
  
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::shearStress >();
        arrayView1d< real64 const > const tauDot = subRegion.getField< inducedSeismicity::shearStressRate >();
  
        // solve for logarithm of seismicity rate
        real64 tol = 1e-8;

        real64 deltaTau = m_directEffect * m_initialSigma;
        real64 c = 1e-3;
        forAll< parallelDevicePolicy<> >(  h.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          real64 error = 100;
          real64 h_prev = h[k];
  
          while (error > tol) 
          {
            // Hard code log shear stress history
            real64 curTau = deltaTau*std::log(c*(time_n+dt)+1) + m_initialTau;
            real64 curTauDot = c*deltaTau/(c*(time_n+dt)+1);

            // real64 gdot = ((tauDot[k]+m_bStressRate)*(sig[k]-p[k]) + (tau[k]+m_bStressRate*(time_n+dt))*pDot[k]) / 
            //                   (m_directEffect*std::pow((sig[k]-p[k]), 2));
            real64 gdot = ((curTauDot+m_bStressRate)*(sig[k]-p[k]) + (curTau+m_bStressRate*(time_n+dt))*pDot[k]) / 
                            (m_directEffect*std::pow((sig[k]-p[k]), 2));

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

// Solve backward-Euler discretization of ODE by Newton-Raphson
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
      
        arrayView1d< real64 const > const p = subRegion.getField< inducedSeismicity::pressure >();
        arrayView1d< real64 const > const p_n = subRegion.getField< inducedSeismicity::pressure_n >();
  
        arrayView1d< real64 const > const sig = subRegion.getField< inducedSeismicity::normalStress >();
        arrayView1d< real64 const > const sig_n = subRegion.getField< inducedSeismicity::normalStress_n >();
  
        arrayView1d< real64 const > const tau = subRegion.getField< inducedSeismicity::shearStress >();
        arrayView1d< real64 const > const tau_n = subRegion.getField< inducedSeismicity::shearStress_n >();

          
        // solve for logarithm of seismicity rate
        real64 deltaTau = m_directEffect * m_initialSigma;
        real64 c = 1e-3;

        forAll< parallelDevicePolicy<> >(  h.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
        {
          // Hard code log shear stress history
          real64 curTau = deltaTau*std::log(c*(time_n+dt)+1) + m_initialTau;
          real64 curTau_n = deltaTau*std::log(c*time_n+1) + m_initialTau;

          real64 g = (curTau + m_bStressRate*(time_n+dt))/(m_directEffect*(sig[k]-p[k])) - m_initialTau/(m_directEffect*m_initialSigma);
          real64 g_n = (curTau_n + m_bStressRate*time_n)/(m_directEffect*(sig_n[k]-p_n[k])) - m_initialTau/(m_directEffect*m_initialSigma);

          real64 deltaLogDenom = std::log(1 + dt/(2*t_a[k])*(std::exp(g - logDenom_n[k]) + std::exp(g_n - logDenom_n[k]) ));
          logDenom[k] = logDenom_n[k] + deltaLogDenom;
  
          h[k] = g - logDenom[k];
          R[k] = LvArray::math::exp( h[k] );
          
          logDenom_n[k] = logDenom[k];
        } );
      } );
    } );
}

/* SETUP SYSTEM
   Setting up the system using the base class method
 */
void DieterichSeismicityRate::setupSystem( DomainPartition & domain,
                                           DofManager & dofManager,
                                           CRSMatrix< real64, globalIndex > & localMatrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution,
                                           bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    string const dofKey = dofManager.getKey( m_fieldName );
    arrayView1d< globalIndex const > const &
    dofIndex = nodeManager.getReference< globalIndex_array >( dofKey );

    SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                    dofManager.numGlobalDofs(),
                                                    8*8*3 );

    finiteElement::fillSparsity< CellElementSubRegion,
                                 DieterichSeismicityRateKernel >( mesh,
                                                     regionNames,
                                                     this->getDiscretizationName(),
                                                     dofIndex,
                                                     dofManager.rankOffset(),
                                                     sparsityPattern );

    sparsityPattern.compress();
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );
  } );
}


/*
   ASSEMBLE SYSTEM
   This is the most important method to assemble the matrices needed before sending them to our solver.
   For a system A.x = B (with x the unknown), here, we use:
   - A : "localMatrix" this represents a Compressed Row Storage (optimized for sparse) matrix of real64 values associated with their index,
   - B : "localRhs" this represents a vector (1d array) of real64 numbers specified at the equation's right-hand side.
   The "local" prefix indicates that we are working on a local problem here, and the parallelization is performed at a higher level.
   This assembly step collects all the information needed to create the matrices localMatrix and localRhs, and the computation of values
   is done in a specific Laplace kernel optimized for parallel performance. Here we:
   1 - identify and point to the mesh of this domain,
   2 - find the node manager of this mesh,
   3 - extract the indices of the nodes that will be solved for (ie. the degrees of freedom or "dof")
   4 - pass all this information to a Laplace-specific finite element computation kernel.
   The call to the kernel is a templated call designed for performance (we will not explain the kernel here).
   See the implementation in LaplaceFEMKernel.cpp.
 */
//START_SPHINX_INCLUDE_ASSEMBLY
void DieterichSeismicityRate::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                 real64 const GEOS_UNUSED_PARAM( dt ),
                                 DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    string const dofKey = dofManager.getKey( m_fieldName );
    arrayView1d< globalIndex const > const &
    dofIndex =  nodeManager.getReference< array1d< globalIndex > >( dofKey );

    DieterichSeismicityRateKernelFactory kernelFactory( dofIndex, dofManager.rankOffset(), localMatrix, localRhs, m_fieldName );

    string const dummyString = "dummy";
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< >,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            this->getDiscretizationName(),
                                                            dummyString,
                                                            kernelFactory );

  } );

}
//END_SPHINX_INCLUDE_ASSEMBLY

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
      arrayView1d< real64 > const tempTa = subRegion.getField< inducedSeismicity::t_a >();
      tempTa.setValues< parallelHostPolicy >( m_directEffect*m_initialSigma/m_bStressRate );
      arrayView1d< real64 > const tempASig = subRegion.getField< inducedSeismicity::aSigma >();
      tempASig.setValues< parallelHostPolicy >( m_directEffect*m_initialSigma );
    
      // Hard coded stressing histories for now
      arrayView1d< real64 > const tempP = subRegion.getField< inducedSeismicity::pressure >();
      tempP.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempP_n = subRegion.getField< inducedSeismicity::pressure_n >();
      tempP_n.setValues< parallelHostPolicy >( 0.0 );
      arrayView1d< real64 > const tempPDot = subRegion.getField< inducedSeismicity::pressureRate >();
      tempPDot.setValues< parallelHostPolicy >( 0.0 );

      arrayView1d< real64 > const tempSig = subRegion.getField< inducedSeismicity::normalStress >();
      tempSig.setValues< parallelHostPolicy >( m_initialSigma );
      arrayView1d< real64 > const tempSig_n = subRegion.getField< inducedSeismicity::normalStress_n >();
      tempSig_n.setValues< parallelHostPolicy >( m_initialSigma );
      arrayView1d< real64 > const tempSigDot = subRegion.getField< inducedSeismicity::normalStressRate >();
      tempSigDot.setValues< parallelHostPolicy >( 0.0 );

      arrayView1d< real64 > const tempTau = subRegion.getField< inducedSeismicity::shearStress >();
      tempTau.setValues< parallelHostPolicy >( m_initialTau );
      arrayView1d< real64 > const tempTau_n = subRegion.getField< inducedSeismicity::shearStress_n >();
      tempTau_n.setValues< parallelHostPolicy >( m_initialTau );
      arrayView1d< real64 > const tempTauDot = subRegion.getField< inducedSeismicity::shearStressRate >();
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

    } );
  } );
}

//START_SPHINX_INCLUDE_REGISTER
REGISTER_CATALOG_ENTRY( SolverBase, DieterichSeismicityRate, string const &, Group * const )
//END_SPHINX_INCLUDE_REGISTER
} /* namespace geos */
