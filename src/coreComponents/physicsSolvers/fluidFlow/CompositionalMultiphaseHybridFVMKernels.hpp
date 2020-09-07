/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "finiteVolume/HybridFVMInnerProduct.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{

/******************************** UpwindingHelper ********************************/

struct UpwindingHelper
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

  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    UpwindViscousTerm( localIndex const er, localIndex const esr, localIndex const ei,
                       localIndex const ifaceLoc,
                       ElementView< arrayView3d< real64 const > > const & phaseDens,
                       ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                       ElementView< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                       ElementView< arrayView2d< real64 const > > const & phaseMob,
                       ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                       ElementView< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                       ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                       ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
                       ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                       ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                       ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
                       real64 ( &upwPhaseViscCoef )[ NF ][ NP ][ NC ],
                       real64 ( &dUpwPhaseViscCoef_dPres )[ NF ][ NP ][ NC ],
                       real64 ( &dUpwPhaseViscCoef_dCompDens )[ NF ][ NP ][ NC ][ NC ],
                       globalIndex ( &upwDofNumber )[ NF ][ NP ] );

};

struct AssemblerKernelHelper
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
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                              arrayView1d< real64 const > const & dFacePres,
                              arrayView1d< real64 const > const & faceGravCoef,
                              arraySlice1d< localIndex const > const & elemToFaces,
                              real64 const & elemPres,
                              real64 const & dElemPres,
                              real64 const & elemGravCoef,
                              arraySlice1d< real64 const > const & elemPhaseDens,
                              arraySlice1d< real64 const > const & dElemPhaseDens_dPres,
                              arraySlice2d< real64 const > const & dElemPhaseDens_dCompFrac,
                              arraySlice2d< real64 const > const & dElemCompFrac_dCompDens,
                              arraySlice2d< real64 const > const & transMatrix,
                              real64 ( &oneSidedVolFlux )[ NF ],
                              real64 ( &dOneSidedVolFlux_dPres )[ NF ],
                              real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                              real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] );

  /**
   * @brief In a given element, collect the upwinded mobilities at this element's faces
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof numbers of all the cell centered pressures (non-local)
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[inout] upwMobility the upwinded mobilities at this element's faces
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or
   * neighbor)
   * @param[inout] upwDofNumber  the dof number of the upwind pressure
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    UpdateUpwindedCoefficients( localIndex const er, localIndex const esr, localIndex const ei,
                                arrayView2d< localIndex const > const & elemRegionList,
                                arrayView2d< localIndex const > const & elemSubRegionList,
                                arrayView2d< localIndex const > const & elemList,
                                SortedArrayView< localIndex const > const & regionFilter,
                                arraySlice1d< localIndex const > const & elemToFaces,
                                ElementView< arrayView3d< real64 const > > const & phaseDens,
                                ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                                ElementView< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                                ElementView< arrayView2d< real64 const > > const & phaseMob,
                                ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                ElementView< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                                ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
                                ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                                ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                                ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
                                real64 const (&oneSidedVolFlux)[ NF ],
                                real64 ( &upwPhaseViscCoef )[ NF ][ NP ][ NC ],
                                real64 ( &dUpwPhaseViscCoef_dPres )[ NF ][ NP ][ NC ],
                                real64 ( &dUpwPhaseViscCoef_dCompDens )[ NF ][ NP ][ NC ][ NC ],
                                globalIndex ( &upwDofNumber )[ NF ][ NP ] );

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] dt the time step size
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwMobility the upwinded mobilities at this element's faces
   * @param[in] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or
   * neighbor)
   * @param[in] upwDofNumber  the dof number of the upwind pressure
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  AssembleOneSidedMassFluxes( arrayView1d< globalIndex const > const & faceDofNumber,
                              arraySlice1d< localIndex const > const & elemToFaces,
                              globalIndex const elemDofNumber,
                              globalIndex const rankOffset,
                              real64 const (&oneSidedVolFlux)[ NF ],
                              real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                              real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                              real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                              real64 const (&upwPhaseViscCoef)[ NF ][ NP ][ NC ],
                              real64 const (&dUpwPhaseViscCoef_dPres)[ NF ][ NP ][ NC ],
                              real64 const (&dUpwPhaseViscCoef_dCompDens)[ NF ][ NP ][ NC ][ NC ],
                              globalIndex const (&upwDofNumber)[ NF ][ NP ],
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                       arrayView1d< integer const > const & faceGhostRank,
                       arraySlice1d< localIndex const > const & elemToFaces,
                       globalIndex const elemDofNumber,
                       globalIndex const rankOffset,
                       real64 const (&oneSidedVolFlux)[ NF ],
                       real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs );

};

/******************************** AssemblerKernel ********************************/

struct AssemblerKernel
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
   * @brief In a given element, assemble the mass conservation equation and the contribution of this element to the face
   * constraints
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank ghost rank of each face
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] eleomGravCoef the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] mobility the mobilities in the domain (non-local)
   * @param[in] dMobility_dPres the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const er, localIndex const esr, localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & dFacePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           ElementView< arrayView3d< real64 const > > const & phaseDens,
           ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementView< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
           ElementView< arrayView2d< real64 const > > const & phaseMob,
           ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementView< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
           ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
           ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs );

};


/******************************** FluxKernel ********************************/

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
   * @brief In a given element, assemble the mass conservation equation and the contribution of this element to the face
   * constraints
   * @param[in] numPhases number of fluid phases
   * @param[in] numComponents number of components
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePotential the phase potentials at the mesh faces at the beginning of the time step
   * @param[in] dFacePotential the accumulated potential updates at the mesh face
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dPres the derivative of the density wrt pressure at this element's center
   * @param[in] dElemDens_dComp the derivative of the density wrt component densities at this element's center
   * @param[in] elemPhaseMobility the mobilities in the domain (non-local)
   * @param[in] dElemPhaseMobility_dPres the derivatives of the mobilities in the domain wrt cell-centered pressure
   *(non-local)
   * @param[in] dElemPhaseMobility_dComp the derivatives of the mobilities in the domain wrt cell-centered component
   * densities (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  static void
  Launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          arrayView1d< real64 const > const & faceGravCoef,
          ElementView< arrayView3d< real64 const > > const & phaseDens,
          ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementView< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
          ElementView< arrayView2d< real64 const > > const & phaseMob,
          ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementView< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
          ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
          ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< typename POLICY, typename REDUCE_POLICY >
  static void
  Launch( arrayView1d< real64 const > const & localResidual,
          globalIndex const rankOffset,
          localIndex const numPhases,
          arrayView1d< globalIndex const > const & facePresDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ElementView< arrayView1d< real64 const > > const & elemVolume,
          ElementView< arrayView1d< real64 const > > const & totalDensOld,
          ElementView< arrayView2d< real64 const > > const & phaseMobOld,
          real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( facePresDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      // if not ghost face and if adjacent to target region
      if( faceGhostRank[iface] < 0 && facePresDofNumber[iface] >= 0 )
      {
        real64 normalizer = 0;
        localIndex elemCounter = 0;
        for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[iface][k];
          localIndex const esr = elemSubRegionList[iface][k];
          localIndex const ei  = elemList[iface][k];

          bool const onBoundary = (er == -1 || esr == -1 || ei == -1);

          // if not on boundary, save the mobility and the upwDofNumber
          if( !onBoundary )
          {
            real64 totalMobOld = 0.0;
            for( localIndex ip = 0; ip < numPhases; ++ip )
            {
              totalMobOld += phaseMobOld[er][esr][ei][ip];
            }
            normalizer += totalDensOld[er][esr][ei] * elemVolume[er][esr][ei]
                          / totalMobOld;
            elemCounter++;
          }
        }
        normalizer /= elemCounter;

        localIndex const lid = facePresDofNumber[iface] - rankOffset;
        real64 const val = localResidual[lid] / normalizer; // to get something dimensionless
        sumScaled += val * val;
      }
    } );

    localResidualNorm = localResidualNorm + sumScaled.get();
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  Launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      if( ghostRank[iface] < 0 )
      {
        localIndex const localRow = dofNumber[iface] - rankOffset;
        {
          real64 const newFacePres = facePres[iface] + dFacePres[iface] + scalingFactor * localSolution[localRow];
          check.min( newFacePres >= 0.0 );
        }
      }
    } );
    return check.get();
  }

};

/******************************** Kernel switches ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void KernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    default: GEOSX_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector( localIndex numFacesInElem, localIndex numComps, localIndex numPhases, ARGS && ... args )
{
  internal::KernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
  {
    if( numPhases == 2 )
    {
      if( numComps == 2 )
      {
        KERNELWRAPPER::template Launch< NF(), 2, 2 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 3 )
      {
        KERNELWRAPPER::template Launch< NF(), 3, 2 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 4 )
      {
        KERNELWRAPPER::template Launch< NF(), 4, 2 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 5 )
      {
        KERNELWRAPPER::template Launch< NF(), 5, 2 >( std::forward< ARGS >( args )... );
      }
      else
      {
        GEOSX_ERROR( "Unsupported number of components: " << numComps );
      }
    }
    else if( numPhases == 3 )
    {
      if( numComps == 2 )
      {
        KERNELWRAPPER::template Launch< NF(), 2, 3 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 3 )
      {
        KERNELWRAPPER::template Launch< NF(), 3, 3 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 4 )
      {
        KERNELWRAPPER::template Launch< NF(), 4, 3 >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 5 )
      {
        KERNELWRAPPER::template Launch< NF(), 5, 3 >( std::forward< ARGS >( args )... );
      }
      else
      {
        GEOSX_ERROR( "Unsupported number of components: " << numComps );
      }
    }
    else
    {
      GEOSX_ERROR( "Unsupported number of phases: " << numPhases );
    }
  } );
}

} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
