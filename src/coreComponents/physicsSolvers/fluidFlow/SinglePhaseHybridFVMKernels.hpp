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
 * @file SinglePhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "finiteVolume/HybridFVMInnerProduct.hpp"

namespace geosx
{

class DofManager;

namespace SinglePhaseHybridFVMKernels
{

/******************************** AssemblerKernelHelper ********************************/

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
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                            arrayView1d< real64 const > const & dFacePres,
                            arrayView1d< real64 const > const & faceGravCoef,
                            arraySlice1d< localIndex const > const elemToFaces,
                            real64 const & elemPres,
                            real64 const & dElemPres,
                            real64 const & elemGravDepth,
                            real64 const & elemDens,
                            real64 const & dElemDens_dp,
                            arraySlice2d< real64 const > const & transMatrix,
                            arraySlice1d< real64 > const & oneSidedVolFlux,
                            arraySlice1d< real64 > const & dOneSidedVolFlux_dp,
                            arraySlice2d< real64 > const & dOneSidedVolFlux_dfp );

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
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  UpdateUpwindedCoefficients( arrayView2d< localIndex const > const & elemRegionList,
                              arrayView2d< localIndex const > const & elemSubRegionList,
                              arrayView2d< localIndex const > const & elemList,
                              SortedArrayView< localIndex const > const & regionFilter,
                              arraySlice1d< localIndex const > const elemToFaces,
                              ElementView< arrayView1d< real64 const > > const & mobility,
                              ElementView< arrayView1d< real64 const > > const & dMobility_dp,
                              ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
                              localIndex const er,
                              localIndex const esr,
                              localIndex const ei,
                              arraySlice1d< real64 const > const & oneSidedVolFlux,
                              arraySlice1d< real64 > const & upwMobility,
                              arraySlice1d< real64 > const & dUpwMobility_dp,
                              arraySlice1d< globalIndex > const & upwDofNumber );

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
  template< localIndex NF >
  //GEOSX_HOST_DEVICE
  static void
  AssembleOneSidedMassFluxes( real64 const & dt,
                              arrayView1d< globalIndex const > const & faceDofNumber,
                              arraySlice1d< localIndex const > const elemToFaces,
                              globalIndex const elemDofNumber,
                              arraySlice1d< real64 const > const & oneSidedVolFlux,
                              arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
                              arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
                              arraySlice1d< real64 const > const & upwMobility,
                              arraySlice1d< real64 const > const & dUpwMobility_dp,
                              arraySlice1d< globalIndex const > const & upwDofNumber,
                              ParallelMatrix * const matrix,
                              ParallelVector * const rhs );


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
  template< localIndex NF >
  //GEOSX_HOST_DEVICE
  static void
  AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                       arraySlice1d< localIndex const > const elemToFaces,
                       globalIndex const elemDofNumber,
                       arraySlice1d< real64 const > const & oneSidedVolFlux,
                       arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
                       arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
                       ParallelMatrix * const matrix,
                       ParallelVector * const rhs );

};

#define INST_AssembleKernelHelper( NF ) \
  extern template \
  void \
  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( arrayView1d< real64 const > const & facePres, \
                                                         arrayView1d< real64 const > const & dFacePres, \
                                                         arrayView1d< real64 const > const & faceGravCoef, \
                                                         arraySlice1d< localIndex const > const elemToFaces, \
                                                         real64 const & elemPres, \
                                                         real64 const & dElemPres, \
                                                         real64 const & elemGravDepth, \
                                                         real64 const & elemDens, \
                                                         real64 const & dElemDens_dp, \
                                                         arraySlice2d< real64 const > const & transMatrix, \
                                                         arraySlice1d< real64 > const & oneSidedVolFlux, \
                                                         arraySlice1d< real64 > const & dOneSidedVolFlux_dp, \
                                                         arraySlice2d< real64 > const & dOneSidedVolFlux_dfp ); \
  extern template \
  void \
  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( real64 const & dt, \
                                                           arrayView1d< globalIndex const > const & faceDofNumber, \
                                                           arraySlice1d< localIndex const > const elemToFaces, \
                                                           globalIndex const elemDofNumber, \
                                                           arraySlice1d< real64 const > const & oneSidedVolFlux, \
                                                           arraySlice1d< real64 const > const & dOneSidedVolFlux_dp, \
                                                           arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp, \
                                                           arraySlice1d< real64 const > const & upwMobility, \
                                                           arraySlice1d< real64 const > const & dUpwMobility_dp, \
                                                           arraySlice1d< globalIndex const > const & upwDofNumber, \
                                                           ParallelMatrix * const matrix, \
                                                           ParallelVector * const rhs ); \
  extern template \
  void \
  AssemblerKernelHelper::AssembleConstraints< NF >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                                    arraySlice1d< localIndex const > const elemToFaces, \
                                                    globalIndex const elemDofNumber, \
                                                    arraySlice1d< real64 const > const & oneSidedVolFlux, \
                                                    arraySlice1d< real64 const > const & dOneSidedVolFlux_dp, \
                                                    arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp, \
                                                    ParallelMatrix * const matrix, \
                                                    ParallelVector * const rhs )

INST_AssembleKernelHelper( 4 );
INST_AssembleKernelHelper( 5 );
INST_AssembleKernelHelper( 6 );

#undef INST_AssembleKernelHelper


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
  template< localIndex NF >
  //GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const er,
           localIndex const esr,
           localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & dFacePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arraySlice1d< localIndex const > const elemToFaces,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           real64 const & elemDens,
           real64 const & dElemDens_dp,
           ElementView< arrayView1d< real64 const > > const & mobility,
           ElementView< arrayView1d< real64 const > > const & dMobility_dp,
           ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
           arraySlice2d< real64 const > const & transMatrix,
           real64 const & dt,
           ParallelMatrix * const matrix,
           ParallelVector * const rhs );

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
   * @brief Assemble the mass conservation equations and face constraints in the cell subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion pointer to the cell element subregion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] mesh the mesh object (single level only)
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemDens the density in the elements of the subregion
   * @param[in] dElemDens_dp the derivative of the density in the elements of the subregion
   * @param[in] mobility the mobilities in the domain (non-local)
   * @param[in] dMobility_dPres the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] lengthTolerance tolerance used in the transmissibility calculations
   * @param[in] dt time step size
   * @param[in] dofManager the dof manager
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF >
  static void
  Launch( localIndex er,
          localIndex esr,
          CellElementSubRegion const & subRegion,
          SortedArrayView< localIndex const > const & regionFilter,
          MeshLevel const & mesh,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView2d< real64 const > const & elemDens,
          arrayView2d< real64 const > const & dElemDens_dp,
          ElementView< arrayView1d< real64 const > > const & mobility,
          ElementView< arrayView1d< real64 const > > const & dMobility_dp,
          real64 const lengthTolerance,
          real64 const dt,
          DofManager const * const dofManager,
          ParallelMatrix * const matrix,
          ParallelVector * const rhs );


};

#define INST_FluxKernel( NF ) \
  extern template \
  void FluxKernel::Launch< NF >( localIndex er, \
                                 localIndex esr, \
                                 CellElementSubRegion const & subRegion, \
                                 SortedArrayView< localIndex const > const & regionFilter, \
                                 MeshLevel const & mesh, \
                                 arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                 arrayView2d< localIndex const > const & elemRegionList, \
                                 arrayView2d< localIndex const > const & elemSubRegionList, \
                                 arrayView2d< localIndex const > const & elemList, \
                                 ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                 arrayView1d< globalIndex const > const & faceDofNumber, \
                                 arrayView1d< real64 const > const & facePres, \
                                 arrayView1d< real64 const > const & dFacePres, \
                                 arrayView1d< real64 const > const & faceGravCoef, \
                                 arrayView2d< real64 const > const & elemDens, \
                                 arrayView2d< real64 const > const & dElemDens_dp, \
                                 ElementView< arrayView1d< real64 const > > const & mobility, \
                                 ElementView< arrayView1d< real64 const > > const & dMobility_dp, \
                                 real64 const lengthTolerance, \
                                 real64 const dt, \
                                 DofManager const * const dofManager, \
                                 ParallelMatrix * const matrix, \
                                 ParallelVector * const rhs )

INST_FluxKernel( 4 );
INST_FluxKernel( 5 );
INST_FluxKernel( 6 );

#undef INST_FluxKernel


namespace helpers
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
    default: GEOSX_ERROR( "Unknown numFaceInElem value: " << value );
  }
}

} // namespace helpers

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector( localIndex numFacesInElem, ARGS && ... args )
{
  helpers::KernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
  {
    KERNELWRAPPER::template Launch< NF() >( std::forward< ARGS >( args )... );
  } );
}


} // namespace SinglePhaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
