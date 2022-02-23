/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

//TJ
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "managers/ProblemManager.hpp"
#include "dataRepository/Group.hpp"
#include "managers/FieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "physicsSolvers/multiphysics/HydrofractureSolver.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geosx
{

namespace SinglePhaseFVMKernels
{

/******************************** FluxKernel ********************************/

/**
 * @struct Structure to contain helper functions for the FluxKernel struct.
 */
struct FluxKernelHelper
{

  /**
   * @tparam INTEGRATION_OPTION This specifies the choice of integration rule for the aperture term
   *         in a lubrication permeability.
   * @param[in] aper0 The beginning of step aperture
   * @param[in] aper The current approximation to the end of step aperture
   * @param[out] aperTerm The resulting
   * @param[out] dAperTerm_dAper
   *
   * Typically in lubrication theory, the permeabilty involves a \f$ aperture^3 \f$ term. The
   * flow residual equation assumes a constant value for all parameters over \f$ dt \f$, which
   * may introduce significant errors given the highly nonlinear nature of the cubic aperture term.
   * The template parameter provides options:
   *  - (0) Forward Euler. This results in no non-linearity since the beginning of step aperture
   *    does not change.
   *  - (1) Exact/Simpson's Rule. This is the result of taking
   *    \f$ \int_{0}^{1} (aperture0 + (aperture-aperture0)x)x^3 dx \f$. This results
   *    in a cubic non-linearity in the resulting set of equations.
   *  .
   *  @note The use of option (1) does not imply that the time integration of the residual
   *        equation is exact, or applying Simpson's Rule. It only means that the integral of
   *        the cubic aperture term in the permeablity is exact. All other components of the
   *        residual equation are assumed constant over the integral, or use a backward
   *        Euler as the case may be. Also, we omit the use of a backward Euler option as
   *        it offers no benefit over the exact integration.
   */
  template< int INTEGRATION_OPTION >
  void static apertureForPermeablityCalculation( real64 const aper0,
                                                 real64 const aper,
                                                 real64 & aperTerm,
                                                 real64 & dAperTerm_dAper );


};

template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 0 >( real64 const aper0,
                                                          real64 const,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = aper0*aper0*aper0;
  dAperTerm_dAper = 0.0;
}

template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 1 >( real64 const aper0,
                                                          real64 const aper,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = 0.25 * ( aper0*aper0*aper0 +
                      aper0*aper0*aper +
                      aper0*aper*aper +
                      aper*aper*aper );

  dAperTerm_dAper = 0.25 * ( aper0*aper0 +
                             2*aper0*aper +
                             3*aper*aper );


  //printf( "aper0, aper, Kf = %4.2e, %4.2e, %4.2e \n", aper0, aper, aperTerm );
}


template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 2 >( real64 const,
                                                          real64 const aper,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = aper*aper*aper;
  dAperTerm_dAper = 3.0*aper*aper;
}


struct FluxKernel
{
  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  /**
   * @brief The type for element-based constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::MaterialViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using MaterialView = typename ElementRegionManager::MaterialViewAccessor< VIEWTYPE >::ViewTypeConst;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] dPres The change in pressure for each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[out] jacobian The linear system matrix
   * @param[out] residual The linear system residual
   */
  template< typename STENCIL_TYPE >
  static void
    Launch( STENCIL_TYPE const & stencil,
            real64 const dt,
            ElementView< arrayView1d< globalIndex const > > const & dofNumber,
            ElementView< arrayView1d< real64 const > > const & pres,
            ElementView< arrayView1d< real64 const > > const & dPres,
            ElementView< arrayView1d< real64 const > > const & gravCoef,
            ElementView< arrayView2d< real64 const > > const & dens,
            ElementView< arrayView2d< real64 const > > const & dDens_dPres,
            ElementView< arrayView1d< real64 const > > const & mob,
            ElementView< arrayView1d< real64 const > > const & dMob_dPres,
            ElementView< arrayView1d< real64 const > > const & aperture0,
            ElementView< arrayView1d< real64 const > > const & aperture,
            ElementView< arrayView1d< R1Tensor const > > const & transTMultiplier,
            R1Tensor const gravityVector,
            real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
            ElementView< arrayView1d< real64 const > > const & s,
            ElementView< arrayView1d< real64 const > > const & dSdAper,
#endif
            ParallelMatrix * const jacobian,
            ParallelVector * const residual,
            CRSMatrixView< real64, localIndex > const & dR_dAper,
	    DomainPartition const * const domain);


  /**
   * @brief Compute flux and its derivatives for a given connection
   *
   * This is a general version that assumes different element regions.
   * See below for a specialized version for fluxes within a region.
   */
  inline static void
  Compute( localIndex const stencilSize,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           ElementView< arrayView1d< real64 const > > const & pres,
           ElementView< arrayView1d< real64 const > > const & dPres,
           ElementView< arrayView1d< real64 const > > const & gravCoef,
           ElementView< arrayView2d< real64 const > > const & dens,
           ElementView< arrayView2d< real64 const > > const & dDens_dPres,
           ElementView< arrayView1d< real64 const > > const & mob,
           ElementView< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {
    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d< real64, numElems >   densWeight( numElems );
    stackArray1d< real64, maxStencil > dDensMean_dP( stencilSize );
    stackArray1d< real64, maxStencil > dFlux_dP( stencilSize );

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 0.5;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for( localIndex ke = 0; ke < numElems; ++ke )
    {
      // density
      real64 const density = dens[seri[ke]][sesri[ke]][sei[ke]][0];
      real64 const dDens_dP = dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0];

      // average density
      densMean        += densWeight[ke] * density;
      dDensMean_dP[ke] = densWeight[ke] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    real64 sumWeightGrav = 0.0;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const er  = seri[ke];
      localIndex const esr = sesri[ke];
      localIndex const ei  = sei[ke];

      real64 weight = stencilWeights[ke];

      real64 const gravD = gravCoef[er][esr][ei];
      real64 const gravTerm = densMean * gravD;
      sumWeightGrav += weight * gravD;

      potDif += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex er_up  = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up  = sei[k_up];

    real64 const mobility     = mob[er_up][esr_up][ei_up];
    real64 const dMobility_dP = dMob_dPres[er_up][esr_up][ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      real64 const weight = stencilWeights[ke];
      dFlux_dP[ke] = mobility * ( weight - dDensMean_dP[ke] * sumWeightGrav);
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    flux[0] = dt * fluxVal;
    flux[1] = -flux[0];

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      fluxJacobian[0][ke] = dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -fluxJacobian[0][ke];
    }
  }


  /**
   * @brief Compute flux and its derivatives for a given connection
   *.
   * This is a specialized version for fluxes within the same region.
   * See above for a general version.
   */
  inline static void
  Compute( localIndex const stencilSize,
           arraySlice1d< localIndex const > const &,
           arraySlice1d< localIndex const > const &,
           arraySlice1d< localIndex const > const & stencilElementIndices,
           arraySlice1d< real64 const > const & stencilWeights,
           arrayView1d< real64 const > const & pres,
           arrayView1d< real64 const > const & dPres,
           arrayView1d< real64 const > const & gravCoef,
           arrayView2d< real64 const > const & dens,
           arrayView2d< real64 const > const & dDens_dPres,
           arrayView1d< real64 const > const & mob,
           arrayView1d< real64 const > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {
    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d< real64, numElems > densWeight( numElems );
    stackArray1d< real64, maxStencil > dDensMean_dP( stencilSize );
    stackArray1d< real64, maxStencil > dFlux_dP( stencilSize );

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 1.0 / numElems;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for( localIndex i = 0; i < numElems; ++i )
    {
      // density
      real64 const density = dens[stencilElementIndices[i]][0];
      real64 const dDens_dP = dDens_dPres[stencilElementIndices[i]][0];

      // average density
      densMean += densWeight[i] * density;
      dDensMean_dP[i] = densWeight[i] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    real64 sumWeightGrav = 0.0;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const ei = stencilElementIndices[ke];
      real64 const weight = stencilWeights[ke];

      real64 const gravD = gravCoef[ei];
      real64 const gravTerm = densMean * gravD;
      sumWeightGrav += weight * gravD;
      potDif += weight * (pres[ei] + dPres[ei] - gravTerm);
    }



    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex ei_up  = stencilElementIndices[k_up];

    real64 const mobility     = mob[ei_up];
    real64 const dMobility_dP = dMob_dPres[ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      real64 const weight = stencilWeights[ke];
      dFlux_dP[ke] = mobility * ( weight - dDensMean_dP[ke] * sumWeightGrav);
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    flux[0] = dt * fluxVal;
    flux[1] = -flux[0];

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      fluxJacobian[0][ke] = dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -fluxJacobian[0][ke];
    }
  }



  /**
   * @brief Compute flux and its derivatives for a given multi-element connector.
   *
   * This is a specialized version that flux in a single region, and uses
   * element pairing instead of a proper junction.
   */
  inline static void
  ComputeJunction( localIndex const numFluxElems,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
                   arraySlice1d< real64 const > const & stencilWeights,
                   arrayView1d< real64 const > const & pres,
                   arrayView1d< real64 const > const & dPres,
                   arrayView1d< real64 const > const & gravCoef,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDens_dPres,
                   arrayView1d< real64 const > const & mob,
                   arrayView1d< real64 const > const & dMob_dPres,
                   arrayView1d< real64 const > const & aperture0,
                   arrayView1d< real64 const > const & aperture,
                   real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                   arrayView1d< real64 const > const &,//s,
                   arrayView1d< real64 const > const &,//dSdAper,
#endif
                   real64 const dt,
                   arraySlice1d< real64 > const & flux,
                   arraySlice2d< real64 > const & fluxJacobian,
                   arraySlice2d< real64 > const & dFlux_dAperture,
		   DomainPartition const * const domain,
		   localIndex const iconn)
  {
    real64 sumOfWeights = 0;
    real64 aperTerm[10];
    real64 dAperTerm_dAper[10];

    std::cout << iconn << std::endl;

    MeshLevel const * mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
    NodeManager const * myNodeManager = mesh->getNodeManager();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp = myNodeManager->totalDisplacement();

    dataRepository::Group const * myProblemManager = domain->getParent();
    PhysicsSolverManager const * myPhysicsSolverManager = myProblemManager->GetGroup< PhysicsSolverManager >("Solvers");
    SurfaceGenerator const * mySurface = myPhysicsSolverManager
				       ->GetGroup< SurfaceGenerator >( "SurfaceGen" );

    HydrofractureSolver const * myHydroSolver = myPhysicsSolverManager
	                               ->GetGroup< HydrofractureSolver >( "hydrofracture" );
    real64 const meshSize = myHydroSolver->getMeshSize();


    SortedArray< localIndex > const trailingFaces = mySurface->getTrailingFaces();
    SortedArray< localIndex > const tipNodes = mySurface->getTipNodes();


    dataRepository::Group const * elementSubRegions = domain->GetGroup("MeshBodies")
					->GetGroup<MeshBody>("mesh1")
					->GetGroup<MeshLevel>("Level0")
					->GetGroup<ElementRegionManager>("ElementRegions")
					->GetRegion< FaceElementRegion >( "Fracture" )
					->GetGroup("elementSubRegions");

    FaceElementSubRegion const * subRegion = elementSubRegions->GetGroup< FaceElementSubRegion >( "default" );
    FaceElementSubRegion::FaceMapType const & faceMap = subRegion->faceList();

    FaceManager const * myFaceManager = mesh->getFaceManager();
    arrayView1d< R1Tensor const > const & faceNormal = myFaceManager->faceNormal();
    arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = myFaceManager->nodeList().toViewConst();


    real64 const shearModulus = domain->GetGroup("Constitutive")
			                  ->GetGroup("rock")
			                  ->getReference<real64>("defaultShearModulus");
    real64 const bulkModulus = domain->GetGroup("Constitutive")
		                         ->GetGroup("rock")
				         ->getReference<real64>("defaultBulkModulus");
    real64 const toughness = mySurface->getReference<real64>("rockToughness");
    real64 const viscosity = domain->GetGroup("Constitutive")
                                   ->GetGroup("water")
				       ->getReference<real64>("defaultViscosity");


    // The unit of injectionRate is kg per second
    real64 const injectionRate = domain->getParent()
                                       ->GetGroup<FieldSpecificationManager>("FieldSpecifications")
                                       ->GetGroup<SourceFluxBoundaryCondition>("sourceTerm")
					   ->getReference<real64>("scale");

    // The injectionRate is only for half domain of the KGD problem,
    // to retrieve the full injection rate, we need to multiply it by 2.0
    real64 const q0 = 2.0 * std::abs(injectionRate) /1.0e3;
    real64 const totalTime = myHydroSolver->getTotalTime();

    real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
    real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
    real64 const Eprime = E/(1.0-nu*nu);
    real64 const PI = 2 * acos(0.0);
    real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;
    real64 const mup = 12.0 * viscosity;

    real64 const dimensionlessK = Kprime/pow(pow(Eprime, 3.0)*mup*q0, 0.25);


    SortedArray< localIndex > tipElements;
    for(auto const & trailingFace : trailingFaces)
    {
      bool found = false;
      // loop over all the face element
      for(localIndex i=0; i<faceMap.size(0); i++)
      {
	// loop over all the (TWO) faces in a face element
	for(localIndex j=0; j<faceMap.size(1); j++)
	{
	  // if the trailingFace is one of the two faces in a face element,
	  // we find it
	  if (faceMap[i][j] == trailingFace)
	  {
	    tipElements.insert(i);
	    found = true;
	    break;
	  }
	} // for localIndex j
	if (found)
	  break;
      } // for localIndex i
    }  // for(auto const & trailingFace : trailingFaces)


    real64 tipLength = -1.0;
    real64 lEr = -1.0;
    localIndex tipElmt = -1;

    for( localIndex k=0; k<numFluxElems; ++k )
    {
      localIndex elmtIndex = stencilElementIndices[k];

      if( std::find( tipElements.begin(), tipElements.end(), elmtIndex ) != tipElements.end() )
      {
        tipElmt = elmtIndex;
      }

      if (tipElmt > -1)
      {
	localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[tipElmt][0] );

	R1Tensor NbarTip = faceNormal[elemsToFaces[tipElmt][0]];
	NbarTip -= faceNormal[elemsToFaces[tipElmt][1]];
	NbarTip.Normalize();

	real64 averageGap = 0.0;

	for (localIndex kf = 0; kf < 2; ++kf)
	{
	  for( localIndex a = 0; a < numNodesPerFace; ++a )
	  {
	    localIndex node = faceToNodeMap( elemsToFaces[tipElmt][kf], a );
	    if ( std::find( tipNodes.begin(), tipNodes.end(), node ) == tipNodes.end() )
	    {
	      R1Tensor temp = disp[node];
	      averageGap += (-pow(-1,kf)) * Dot( temp, NbarTip)/2 ;
	    }
	  }
	}

	if (dimensionlessK >= 4.0) // Toughness-dominated case
	{
	  //TJ: the tip asymptote w = Kprime / Eprime * x^(1/2)
	  tipLength = (averageGap * Eprime / Kprime) * (averageGap * Eprime / Kprime);
	}
	else if (dimensionlessK <= 1.0)// Viscosity-dominated case
	{
	  real64 Lm = pow( Eprime*pow(q0,3.0)*pow(totalTime,4.0)/mup, 1.0/6.0 );
	  real64 gamma_m0 = 0.616;
	  real64 velocity = 2.0/3.0 * Lm * gamma_m0 / totalTime;
	  real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);
	  tipLength = sqrt( Eprime/(mup*velocity) * pow( averageGap/Betam, 3.0) );
	}
	else
	{
          real64 const alpha = (dimensionlessK - 1.0)/3.0;

          real64 Lm = pow( Eprime*pow(q0,3.0)*pow(totalTime,4.0)/mup, 1.0/6.0 );
          real64 gamma_m0 = 0.616;
          real64 velocity = 2.0/3.0 * Lm * gamma_m0 / totalTime;
          real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);

          real64 const A = (1-alpha)*Betam*pow(mup*velocity/Eprime ,1.0/3.0);
          real64 const B = alpha*Kprime/Eprime;

          real64 l0 = 1.0e-9*meshSize;
          real64 dwdl;
          real64 lold = l0;
          real64 deltal;
          real64 wold;
          for (int i = 0; i < 41; i++)
          {
            dwdl = 2.0/3.0*A*pow(lold, -1.0/3.0) + 1.0/2.0*B*pow(lold, -0.5);
            wold = A*pow(lold, 2.0/3.0) + B*pow(lold, 1.0/2.0);
            if (std::abs(wold - averageGap)< 1.0e-6*averageGap)
            {
              break;
            }
            if (i == 40)
            {
              GEOSX_ERROR( "Newton loop cannot converge to solve l from w." );
            }
            deltal = (averageGap - wold)/dwdl;
            lold = lold + deltal;
          }
          l0 = lold;
          tipLength = l0;
	}
      } // if (tipElmt > -1), we have a tip elmt here


#define PERM_CALC 1
//      real64 const aperAdd = aperture0[stencilElementIndices[k]] < 0.09e-3 ? ( 0.09e-3 -
// aperture0[stencilElementIndices[k]] ) : 0.0;
#if PERM_CALC==1
      FluxKernelHelper::
        apertureForPermeablityCalculation< 2 >( aperture0[stencilElementIndices[k]],
                                                aperture[stencilElementIndices[k]],
                                                aperTerm[k],
                                                dAperTerm_dAper[k] );
#elif PERM_CALC==2

      if( s[k] >= 1.0 )
      {
        aperTerm[k] = aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]];
        dAperTerm_dAper[k] = 3*aperture[stencilElementIndices[k]]*aperture[stencilElementIndices[k]];
      }
      else
      {
        aperTerm[k] = aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]]/s[k];
        dAperTerm_dAper[k] = 3*aperture[stencilElementIndices[k]]*aperture[stencilElementIndices[k]]/s[k]
                             - aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]]/(s[k]*s[k]) *
                             dSdAper[k];
      }
#endif
//      aperTerm[k] += aperAdd*aperAdd*aperAdd;

      if (tipElmt > -1) // k is a tip elmt
      {
	//TJ: the real length for pressure gradient
	lEr = 0.5 * tipLength;
	GEOSX_ASSERT_MSG( lEr > 0.0, "tipLength is wrong!");
        sumOfWeights += aperTerm[k] * stencilWeights[k] * 0.5 * meshSize/lEr;
      }
      else  // k is not a tip elmt
      {
        sumOfWeights += aperTerm[k] * stencilWeights[k];
      }
    }

    localIndex k[2];
    for( k[0]=0; k[0]<numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
      {
        real64 dFlux_dP[2] = {0, 0};

        localIndex const ei[2] = { stencilElementIndices[k[0]],
                                   stencilElementIndices[k[1]] };
#if 0
        real64 const weight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
                              ( stencilWeights[k[1]]*aperTerm[k[1]] ) / sumOfWeights;

        real64 const
        dWeight_dAper[2] =
        { ( 1 / aperTerm[k[0]]  - stencilWeights[k[0]] / sumOfWeights ) * weight * dAperTerm_dAper[k[0]],
          ( 1 / aperTerm[k[1]]  - stencilWeights[k[1]] / sumOfWeights ) * weight * dAperTerm_dAper[k[1]]};
#else
        real64 const c = meanPermCoeff;

        real64 harmonicWeight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
                                ( stencilWeights[k[1]]*aperTerm[k[1]] ) / sumOfWeights;

        //TJ: we have a tip elmt in the pair, use the real distance for pressure gradient
        if (tipElmt > -1)
        {
          harmonicWeight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
                           ( stencilWeights[k[1]]*aperTerm[k[1]] ) *
			   0.5*meshSize/lEr
			   / sumOfWeights;
        }

        real64 weight = c * harmonicWeight
                     + (1.0 - c) * 0.25 * ( stencilWeights[k[0]]*aperTerm[k[0]] + stencilWeights[k[1]]*aperTerm[k[1]] );

        if (tipElmt > -1)
        {
          if (ei[0] == tipElmt)
          {
            weight = c * harmonicWeight
                  + (1.0 - c) * 0.25 * ( stencilWeights[k[0]]*aperTerm[k[0]] * 0.5*meshSize/lEr
                                       + stencilWeights[k[1]]*aperTerm[k[1]] );
          }
          else
          {
            weight = c * harmonicWeight
                  + (1.0 - c) * 0.25 * ( stencilWeights[k[0]]*aperTerm[k[0]]
				       + stencilWeights[k[1]]*aperTerm[k[1]] * 0.5*meshSize/lEr );
          }
        }

        real64
        dHarmonicWeight_dAper[2] =
        { ( 1 / aperTerm[k[0]]  - stencilWeights[k[0]] / sumOfWeights ) * harmonicWeight * dAperTerm_dAper[k[0]],
          ( 1 / aperTerm[k[1]]  - stencilWeights[k[1]] / sumOfWeights ) * harmonicWeight * dAperTerm_dAper[k[1]]};

        if (tipElmt > -1)
        {
          real64 dExtraTerm;
	  if (dimensionlessK >= 4.0) // Toughness-dominated case
	  {
	    dExtraTerm = -2.0 * 0.5*meshSize / lEr / aperture[tipElmt];
	  }
	  else if (dimensionlessK <= 1.0) // Viscosity-dominated case
	  {
	    dExtraTerm = -3.0/2.0 * 0.5*meshSize / lEr / aperture[tipElmt];
	  }
	  else
	  {
            real64 const alpha = (dimensionlessK - 1.0)/3.0;

            real64 Lm = pow( Eprime*pow(q0,3.0)*pow(totalTime,4.0)/mup, 1.0/6.0 );
            real64 gamma_m0 = 0.616;
            real64 velocity = 2.0/3.0 * Lm * gamma_m0 / totalTime;
            real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);

            real64 const C = (1-alpha)*3.0/5.0*Betam*pow( mup*velocity/Eprime, 1.0/3.0);
            real64 const D = alpha*2.0/3.0*Kprime/Eprime;

            real64 l0 = 1.0e-9*meshSize;
            real64 dwdl;
            real64 lold = l0;
            real64 deltal;
            real64 wold;
            for (int i = 0; i < 41; i++)
            {
              dwdl = 2.0/3.0*C*pow(lold, -1.0/3.0) + 1.0/2.0*D*pow(lold, -0.5);
              wold = C*pow(lold, 2.0/3.0) + D*pow(lold, 1.0/2.0);
              if (std::abs(wold - aperture[tipElmt])< 1.0e-6*aperture[tipElmt])
              {
                break;
              }
              if (i == 40)
              {
                GEOSX_ERROR( "Newton loop cannot converge to solve l from w_tip." );
              }
              deltal = (aperture[tipElmt] - wold)/dwdl;
              lold = lold + deltal;
            }
            l0 = lold;
            real64 coeff = 0.5 / (2.0/3.0*C*pow(l0, -1.0/3.0) + 1.0/2.0*D*pow(l0, -0.5));
            dExtraTerm = -0.5*meshSize/pow(0.5*l0, 2.0)*coeff;
	  }

          if (ei[0] == tipElmt)
          {
            real64 term1;
            real64 term2;
            term1 = ( stencilWeights[k[0]] * dAperTerm_dAper[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * aperTerm[k[1]]
                  +   stencilWeights[k[0]] * aperTerm[k[0]] * dExtraTerm * stencilWeights[k[1]] * aperTerm[k[1]] )
		  / sumOfWeights;

            term2 = -( stencilWeights[k[0]] * dAperTerm_dAper[k[0]] * 0.5*meshSize/lEr
        	  +    stencilWeights[k[0]] * aperTerm[k[0]] * dExtraTerm ) / pow(sumOfWeights, 2.0)
		  * stencilWeights[k[0]] * aperTerm[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * aperTerm[k[1]];

            dHarmonicWeight_dAper[0] = term1 + term2;

            term1 = stencilWeights[k[0]] * aperTerm[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * dAperTerm_dAper[k[1]]
		  / sumOfWeights;

            term2 = -stencilWeights[k[1]] * dAperTerm_dAper[k[1]] / pow(sumOfWeights, 2.0)
		  * stencilWeights[k[0]] * aperTerm[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * aperTerm[k[1]];

            dHarmonicWeight_dAper[1] = term1 + term2;
          }
          else
          {
            real64 term1;
            real64 term2;

            term1 = stencilWeights[k[0]] * dAperTerm_dAper[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * aperTerm[k[1]]
                  / sumOfWeights;

            term2 = -stencilWeights[k[0]] * dAperTerm_dAper[k[0]] / pow(sumOfWeights, 2.0)
            	  * stencilWeights[k[0]] * aperTerm[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * aperTerm[k[1]];

            dHarmonicWeight_dAper[0] = term1 + term2;

            term1 = ( stencilWeights[k[0]] * aperTerm[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * dAperTerm_dAper[k[1]]
                  +   stencilWeights[k[0]] * aperTerm[k[0]] * dExtraTerm * stencilWeights[k[1]] * aperTerm[k[1]] )
		  / sumOfWeights;

            term2 = -( stencilWeights[k[1]] * dAperTerm_dAper[k[1]] * 0.5*meshSize/lEr
                  +    stencilWeights[k[1]] * aperTerm[k[1]] * dExtraTerm ) / pow(sumOfWeights, 2.0)
             	  * stencilWeights[k[0]] * aperTerm[k[0]] * 0.5*meshSize/lEr * stencilWeights[k[1]] * aperTerm[k[1]];

            dHarmonicWeight_dAper[1] = term1 + term2;
          }
        } //  if (tipElmt > -1)

        real64
        dWeight_dAper[2] =
        { c * dHarmonicWeight_dAper[0] + 0.25 * ( 1.0 - c )*stencilWeights[k[0]]*dAperTerm_dAper[k[0]],
          c * dHarmonicWeight_dAper[1] + 0.25 * ( 1.0 - c )*stencilWeights[k[1]]*dAperTerm_dAper[k[1]] };

        if (tipElmt > -1)
        {
          real64 dExtraTerm;
	  if (dimensionlessK >= 4.0) // Toughness-dominated case
	  {
	    dExtraTerm = -2.0 * 0.5*meshSize / lEr / aperture[tipElmt];
	  }
	  else if (dimensionlessK <= 1.0)// Viscosity-dominated case
	  {
	    dExtraTerm = -3.0/2.0 * 0.5*meshSize / lEr / aperture[tipElmt];
	  }
	  else
	  {
            real64 const alpha = (dimensionlessK - 1.0)/3.0;

            real64 Lm = pow( Eprime*pow(q0,3.0)*pow(totalTime,4.0)/mup, 1.0/6.0 );
            real64 gamma_m0 = 0.616;
            real64 velocity = 2.0/3.0 * Lm * gamma_m0 / totalTime;
            real64 Betam = pow(2.0, 1.0/3.0) * pow(3.0, 5.0/6.0);

            real64 const C = (1-alpha)*3.0/5.0*Betam*pow( mup*velocity/Eprime, 1.0/3.0);
            real64 const D = alpha*2.0/3.0*Kprime/Eprime;

            real64 l0 = 1.0e-9*meshSize;
            real64 dwdl;
            real64 lold = l0;
            real64 deltal;
            real64 wold;
            for (int i = 0; i < 41; i++)
            {
              dwdl = 2.0/3.0*C*pow(lold, -1.0/3.0) + 1.0/2.0*D*pow(lold, -0.5);
              wold = C*pow(lold, 2.0/3.0) + D*pow(lold, 1.0/2.0);
              if (std::abs(wold - aperture[tipElmt])< 1.0e-6*aperture[tipElmt])
              {
                break;
              }
              if (i == 40)
              {
                GEOSX_ERROR( "Newton loop cannot converge to solve l from w_tip." );
              }
              deltal = (aperture[tipElmt] - wold)/dwdl;
              lold = lold + deltal;
            }
            l0 = lold;
            real64 coeff = 0.5 / (2.0/3.0*C*pow(l0, -1.0/3.0) + 1.0/2.0*D*pow(l0, -0.5));
            dExtraTerm = -0.5*meshSize/pow(0.5*l0, 2.0)*coeff;
	  }

          if (ei[0] == tipElmt)
          {
            dWeight_dAper[0] = c * dHarmonicWeight_dAper[0]
			     + 0.25 * ( 1.0 - c )*stencilWeights[k[0]]*dAperTerm_dAper[k[0]]*0.5*meshSize/lEr
			     + 0.25 * ( 1.0 - c )*stencilWeights[k[0]]*aperTerm[k[0]]*dExtraTerm;

            dWeight_dAper[1] = c * dHarmonicWeight_dAper[1]
			     + 0.25 * ( 1.0 - c )*stencilWeights[k[1]]*dAperTerm_dAper[k[1]];
          }
          else
          {
            dWeight_dAper[0] = c * dHarmonicWeight_dAper[0]
			     + 0.25 * ( 1.0 - c )*stencilWeights[k[0]]*dAperTerm_dAper[k[0]];

            dWeight_dAper[1] = c * dHarmonicWeight_dAper[1]
            		     + 0.25 * ( 1.0 - c )*stencilWeights[k[1]]*dAperTerm_dAper[k[1]]*0.5*meshSize/lEr
            		     + 0.25 * ( 1.0 - c )*stencilWeights[k[1]]*aperTerm[k[1]]*dExtraTerm;
          }
        }


#endif
        // average density
        real64 const densMean = 0.5 * ( dens[ei[0]][0] + dens[ei[1]][0] );

        real64 const dDensMean_dP[2] = { 0.5 * dDens_dPres[ei[0]][0],
                                         0.5 * dDens_dPres[ei[1]][0] };

        real64 const potDif =  ( ( pres[ei[0]] + dPres[ei[0]] ) - ( pres[ei[1]] + dPres[ei[1]] ) -
                                 densMean * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) );

        // upwinding of fluid properties (make this an option?)
        localIndex const k_up = (potDif >= 0) ? 0 : 1;

        localIndex ei_up  = stencilElementIndices[k[k_up]];

        real64 const mobility     = mob[ei_up];
        real64 const dMobility_dP = dMob_dPres[ei_up];

        // Compute flux and fill flux rval
        real64 const fluxVal = mobility * weight * potDif * dt;
        flux[k[0]] += fluxVal;
        flux[k[1]] -= fluxVal;

        // compute and fill dFlux_dP
        dFlux_dP[0] = mobility * weight * (  1 - dDensMean_dP[0] * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) ) * dt;
        dFlux_dP[1] = mobility * weight * ( -1 - dDensMean_dP[1] * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) ) * dt;
        dFlux_dP[k_up] += dMobility_dP * weight * potDif * dt;

        fluxJacobian[k[0]][k[0]] += dFlux_dP[0];
        fluxJacobian[k[0]][k[1]] += dFlux_dP[1];
        fluxJacobian[k[1]][k[0]] -= dFlux_dP[0];
        fluxJacobian[k[1]][k[1]] -= dFlux_dP[1];

        real64 const dFlux_dAper[2] = { mobility * dWeight_dAper[0] * potDif * dt,
                                        mobility * dWeight_dAper[1] * potDif * dt };
        dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
        dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
        dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
        dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];
      }
    }
  }

};


} // namespace SinglePhaseFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
