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
 * @file SinglePhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "mesh/MeshLevel.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geosx
{
namespace singlePhaseHybridFVMKernels
{

/******************************** Kernel switches ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
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
    case 7:
    { lambda( std::integral_constant< T, 7 >() ); return;}
    case 8:
    { lambda( std::integral_constant< T, 8 >() ); return;}
    case 9:
    { lambda( std::integral_constant< T, 9 >() ); return;}
    case 10:
    { lambda( std::integral_constant< T, 10 >() ); return;}
    case 11:
    { lambda( std::integral_constant< T, 11 >() ); return;}
    case 12:
    { lambda( std::integral_constant< T, 12 >() ); return;}
    case 13:
    { lambda( std::integral_constant< T, 13 >() ); return;}
    default: GEOSX_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_FACE number of faces per element
 * @tparam IP the type of inner product
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_FACE, typename IP >
class ElementBasedAssemblyKernel
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using FlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::mobility,
                      extrinsicMeshData::flow::dMobility_dPressure >;

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] er the index of the region
   * @param[in] esr the index of the subregion
   * @param[in] lengthTolerance tolerance used in the transmissibility computations
   * @param[in] faceDofKey the string key to retrieve the face degrees of freedom numbers
   * @param[in] nodeManager the node manager
   * @param[in] faceManager the face manager
   * @param[in] subRegion the element subregion
   * @param[in] dofNumberAccessor accessors for the element dof numbers
   * @param[in] flowAccessors accessors for the flow variables
   * @param[in] fluid the fluid model
   * @param[in] permeability the permeability model
   * @param[in] regionFilter the region filter
   * @param[in] dt the time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              localIndex const er,
                              localIndex const esr,
                              real64 const & lengthTolerance,
                              string const faceDofKey,
                              NodeManager const & nodeManager,
                              FaceManager const & faceManager,
                              CellElementSubRegion const & subRegion,
                              DofNumberAccessor const & dofNumberAccessor,
                              FlowAccessors const & flowAccessors,
                              constitutive::SingleFluidBase const & fluid,
                              constitutive::PermeabilityBase const & permeability,
                              SortedArrayView< localIndex const > const & regionFilter,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_er( er ),
    m_esr( esr ),
    m_lengthTolerance( lengthTolerance ),
    m_dt( dt ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_elemDofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_faceGhostRank( faceManager.ghostRank() ),
    m_faceDofNumber( faceManager.getReference< array1d< globalIndex > >( faceDofKey ) ),
    m_elemToFaces( subRegion.faceList().toViewConst() ),
    m_elemCenter( subRegion.getElementCenter() ),
    m_elemVolume( subRegion.getElementVolume() ),
    m_elemGravCoef( subRegion.getExtrinsicData< extrinsicMeshData::flow::gravityCoefficient >() ),
    m_faceToNodes( faceManager.nodeList().toViewConst() ),
    m_faceGravCoef( faceManager.getExtrinsicData< extrinsicMeshData::flow::gravityCoefficient >() ),
    m_regionFilter( regionFilter ),
    m_nodePosition( nodeManager.referencePosition() ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_elemPerm( permeability.permeability() ),
    m_transMultiplier( faceManager.getExtrinsicData< extrinsicMeshData::flow::transMultiplier >() ),
    m_elemPres( subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >() ),
    m_facePres( faceManager.getExtrinsicData< extrinsicMeshData::flow::facePressure >() ),
    m_elemDens ( fluid.density() ),
    m_dElemDens_dPres( fluid.dDensity_dPressure() ),
    m_mob( flowAccessors.get( extrinsicMeshData::flow::mobility {} ) ),
    m_dMob_dPres( flowAccessors.get( extrinsicMeshData::flow::dMobility_dPressure {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}


  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {

    GEOSX_HOST_DEVICE
    StackVariables()
      : transMatrix( NUM_FACE, NUM_FACE )
    {}

    stackArray2d< real64, NUM_FACE *NUM_FACE > transMatrix;

    real64 oneSidedVolFlux[NUM_FACE]{};
    real64 dOneSidedVolFlux_dPres[NUM_FACE]{};
    real64 dOneSidedVolFlux_dFacePres[NUM_FACE][NUM_FACE]{};

    real64 divMassFluxes = 0;
    real64 dDivMassFluxes_dElemVars[NUM_FACE+1]{};
    real64 dDivMassFluxes_dFaceVars[NUM_FACE]{};

    localIndex cellCenteredEqnRowIndex = 0;
    localIndex faceCenteredEqnRowIndex[NUM_FACE]{};
    globalIndex elemDofColIndices[NUM_FACE+1]{};
    globalIndex faceDofColIndices[NUM_FACE]{};
  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    stack.cellCenteredEqnRowIndex = m_elemDofNumber[m_er][m_esr][ei] - m_rankOffset;
    stack.elemDofColIndices[0] = m_elemDofNumber[m_er][m_esr][ei];
    for( integer iFaceLoc = 0; iFaceLoc < NUM_FACE; ++iFaceLoc )
    {
      stack.faceCenteredEqnRowIndex[iFaceLoc] = m_faceDofNumber[m_elemToFaces[ei][iFaceLoc]] - m_rankOffset;
      stack.faceDofColIndices[iFaceLoc] = m_faceDofNumber[m_elemToFaces[ei][iFaceLoc]];
    }
  }

  /**
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeGradient( localIndex const ei,
                        StackVariables & stack ) const
  {
    for( integer iFaceLoc = 0; iFaceLoc < NUM_FACE; ++iFaceLoc )
    {
      // now in the following nested loop,
      // we compute the contribution of face jFaceLoc to the one sided total volumetric flux at face iFaceLoc
      for( integer jFaceLoc = 0; jFaceLoc < NUM_FACE; ++jFaceLoc )
      {
        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = m_elemPres[ei];
        real64 const fPres = m_facePres[m_elemToFaces[ei][jFaceLoc]];

        real64 const ccGravCoef = m_elemGravCoef[ei];
        real64 const fGravCoef = m_faceGravCoef[m_elemToFaces[ei][jFaceLoc]];

        real64 const ccDens = m_elemDens[ei][0];
        real64 const dCcDens_dPres = m_dElemDens_dPres[ei][0];
        // no density evaluated at the face center

        // pressure difference
        real64 const presDif = ccPres - fPres;
        real64 const dPresDif_dPres = 1;
        real64 const dPresDif_dFacePres = -1;

        // gravity term
        real64 const gravCoefDif = ccGravCoef - fGravCoef;
        real64 const gravTerm = ccDens * gravCoefDif;
        real64 const dGravTerm_dPres = dCcDens_dPres * gravCoefDif;

        // potential difference
        real64 const potDif = presDif - gravTerm;
        real64 const dPotDif_dPres = dPresDif_dPres - dGravTerm_dPres;
        real64 const dPotDif_dFacePres = dPresDif_dFacePres;

        // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
        stack.oneSidedVolFlux[iFaceLoc] = stack.oneSidedVolFlux[iFaceLoc]
                                          + stack.transMatrix[iFaceLoc][jFaceLoc] * potDif;
        stack.dOneSidedVolFlux_dPres[iFaceLoc] = stack.dOneSidedVolFlux_dPres[iFaceLoc]
                                                 + stack.transMatrix[iFaceLoc][jFaceLoc] * dPotDif_dPres;
        stack.dOneSidedVolFlux_dFacePres[iFaceLoc][jFaceLoc] = stack.dOneSidedVolFlux_dFacePres[iFaceLoc][jFaceLoc]
                                                               + stack.transMatrix[iFaceLoc][jFaceLoc] * dPotDif_dFacePres;
      }
    }
  }

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeFluxDivergence( localIndex const ei,
                              StackVariables & stack ) const
  {
    // upwinded mobility
    real64 upwMobility = 0.0;
    real64 dUpwMobility_dPres = 0.0;
    globalIndex upwDofNumber = 0;

    // for each element, loop over the one-sided faces
    for( integer iFaceLoc = 0; iFaceLoc < NUM_FACE; ++iFaceLoc )
    {

      // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element
      localIndex local[3] = { m_er, m_esr, ei };
      localIndex neighbor[3] = { m_er, m_esr, ei };
      bool const isNeighborFound =
        hybridFVMKernels::CellConnectivity::
          isNeighborFound( local,
                           iFaceLoc,
                           m_elemRegionList,
                           m_elemSubRegionList,
                           m_elemList,
                           m_regionFilter,
                           m_elemToFaces[ei],
                           neighbor );

      // 2) Upwind the mobility at this face
      if( stack.oneSidedVolFlux[iFaceLoc] >= 0 || !isNeighborFound )
      {
        upwMobility = m_mob[m_er][m_esr][ei];
        dUpwMobility_dPres = m_dMob_dPres[m_er][m_esr][ei];
        upwDofNumber = m_elemDofNumber[m_er][m_esr][ei];
      }
      else
      {
        upwMobility = m_mob[neighbor[0]][neighbor[1]][neighbor[2]];
        dUpwMobility_dPres = m_dMob_dPres[neighbor[0]][neighbor[1]][neighbor[2]];
        upwDofNumber = m_elemDofNumber[neighbor[0]][neighbor[1]][neighbor[2]];
      }

      // 3) Add to the flux divergence
      real64 const dt_upwMobility = m_dt * upwMobility;

      // compute the mass flux at the one-sided face plus its derivatives and add the newly computed flux to the sum
      stack.divMassFluxes = stack.divMassFluxes + dt_upwMobility * stack.oneSidedVolFlux[iFaceLoc];
      stack.dDivMassFluxes_dElemVars[0] = stack.dDivMassFluxes_dElemVars[0] + m_dt * upwMobility * stack.dOneSidedVolFlux_dPres[iFaceLoc];
      stack.dDivMassFluxes_dElemVars[iFaceLoc+1] = m_dt * dUpwMobility_dPres * stack.oneSidedVolFlux[iFaceLoc];
      for( integer jFaceLoc = 0; jFaceLoc < NUM_FACE; ++jFaceLoc )
      {
        stack.dDivMassFluxes_dFaceVars[jFaceLoc] = stack.dDivMassFluxes_dFaceVars[jFaceLoc]
                                                   + dt_upwMobility * stack.dOneSidedVolFlux_dFacePres[iFaceLoc][jFaceLoc];
      }

      // collect the relevant dof numbers
      stack.elemDofColIndices[iFaceLoc+1] = upwDofNumber;
    }
  }

  /**
   * @brief Compute the fluxes contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack,
                FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    GEOSX_UNUSED_VAR( ei, stack, kernelOp );

    real64 const perm[ 3 ] = { m_elemPerm[ei][0][0], m_elemPerm[ei][0][1], m_elemPerm[ei][0][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    IP::template compute< NUM_FACE >( m_nodePosition,
                                      m_transMultiplier,
                                      m_faceToNodes,
                                      m_elemToFaces[ei],
                                      m_elemCenter[ei],
                                      m_elemVolume[ei],
                                      perm,
                                      m_lengthTolerance,
                                      stack.transMatrix );

    /*
     * compute auxiliary quantities at the one sided faces of this element:
     * 1) One-sided volumetric fluxes
     * 2) Upwinded mobilities
     */

    // for each one-sided face of the elem,
    // compute the volumetric flux using transMatrix
    computeGradient( ei, stack );

    // at this point, we know the local flow direction in the element
    // so we can upwind the transport coefficients (mobilities) at the one sided faces
    if( m_elemGhostRank[ei] < 0 )
    {

      /*
       * perform assembly in this element in two steps:
       * 1) mass conservation equations
       * 2) face constraints
       */

      // use the computed one sided vol fluxes and the upwinded mobilities
      // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
      computeFluxDivergence( ei, stack );
    }

  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble cell-centered residual and its derivatives

    // we are ready to assemble the local flux and its derivatives
    // no need for atomic adds - each row is assembled by a single thread

    if( m_elemGhostRank[ei] < 0 )
    {

      // residual
      m_localRhs[stack.cellCenteredEqnRowIndex] =
        m_localRhs[stack.cellCenteredEqnRowIndex] + stack.divMassFluxes;

      // jacobian -- derivative wrt elem centered vars
      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( stack.cellCenteredEqnRowIndex,
                                                                  &stack.elemDofColIndices[0],
                                                                  &stack.dDivMassFluxes_dElemVars[0],
                                                                  NUM_FACE+1 );

      // jacobian -- derivatives wrt face centered vars
      m_localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( stack.cellCenteredEqnRowIndex,
                                                                  &stack.faceDofColIndices[0],
                                                                  &stack.dDivMassFluxes_dFaceVars[0],
                                                                  NUM_FACE );

    }

    // Step 2: assemble face-centered residuals and their derivatives

    globalIndex const dofColIndexElemPres = stack.elemDofColIndices[0];

    // for each element, loop over the local (one-sided) faces
    for( integer iFaceLoc = 0; iFaceLoc < NUM_FACE; ++iFaceLoc )
    {
      if( m_faceGhostRank[m_elemToFaces[ei][iFaceLoc]] < 0 )
      {
        // residual
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[stack.faceCenteredEqnRowIndex[iFaceLoc]], stack.oneSidedVolFlux[iFaceLoc] );

        // jacobian -- derivative wrt local cell centered pressure term
        m_localMatrix.addToRow< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[iFaceLoc],
                                                        &dofColIndexElemPres,
                                                        &stack.dOneSidedVolFlux_dPres[iFaceLoc],
                                                        1 );

        // jacobian -- derivatives wrt face pressure terms
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.faceCenteredEqnRowIndex[iFaceLoc],
                                                                            &stack.faceDofColIndices[0],
                                                                            stack.dOneSidedVolFlux_dFacePres[iFaceLoc],
                                                                            NUM_FACE );
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// offset for my MPI rank
  globalIndex const m_rankOffset;

  /// index of the region and sub-region
  localIndex const m_er;
  localIndex const m_esr;

  /// length tolerance
  real64 const m_lengthTolerance;

  /// time step size
  real64 const m_dt;

  /// ghost rank numbers
  arrayView1d< integer const > const m_elemGhostRank;
  ElementViewConst< arrayView1d< globalIndex const > > const m_elemDofNumber;
  arrayView1d< integer const > const m_faceGhostRank;
  arrayView1d< globalIndex const > const m_faceDofNumber;

  /// topological and geometrical data
  arrayView2d< localIndex const > const m_elemToFaces;
  arrayView2d< real64 const > const m_elemCenter;
  arrayView1d< real64 const > const m_elemVolume;
  arrayView1d< real64 const > const m_elemGravCoef;
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  arrayView1d< real64 const > const m_faceGravCoef;

  SortedArrayView< localIndex const > const m_regionFilter;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition;
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;

  /// permeability
  arrayView3d< real64 const > const m_elemPerm;
  arrayView1d< real64 const > const m_transMultiplier;

  /// pressure and fluid data
  arrayView1d< real64 const > const m_elemPres;
  arrayView1d< real64 const > const m_facePres;
  arrayView2d< real64 const > const m_elemDens;
  arrayView2d< real64 const > const m_dElemDens_dPres;
  ElementViewConst< arrayView1d< real64 const > > const m_mob;
  ElementViewConst< arrayView1d< real64 const > > const m_dMob_dPres;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

};

class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] er the index of the region
   * @param[in] esr the index of the subregion
   * @param[in] lengthTolerance tolerance used in the transmissibility computations
   * @param[in] elemDofKey the string key to retrieve the element dof numbers
   * @param[in] faceDofKey the string key to retrieve the face dof numbers
   * @param[in] solverName the name of the solver
   * @param[in] nodeManager the node manager
   * @param[in] faceManager the face manager
   * @param[in] elemManager the element region manager
   * @param[in] subRegion the element sub-region
   * @param[in] mimeticInnerProductBase the inner product to dispatch
   * @param[in] fluid the fluid model
   * @param[in] permeability the permeability model
   * @param[in] regionFilter the region filter
   * @param[in] dt the time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   localIndex const er,
                   localIndex const esr,
                   real64 const lengthTolerance,
                   string const elemDofKey,
                   string const faceDofKey,
                   string const solverName,
                   NodeManager const & nodeManager,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   CellElementSubRegion const & subRegion,
                   mimeticInnerProduct::MimeticInnerProductBase const & mimeticInnerProductBase,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::PermeabilityBase const & permeability,
                   SortedArrayView< localIndex const > const & regionFilter,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    mimeticInnerProductDispatch( mimeticInnerProductBase,
                                 [&] ( auto const mimeticInnerProduct )
    {
      using IP = TYPEOFREF( mimeticInnerProduct );

      internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
      {
        ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
          elemManager.constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
        dofNumberAccessor.setName( solverName + "/accessors/" + elemDofKey );

        using kernelType = ElementBasedAssemblyKernel< NUM_FACES, IP >;
        typename kernelType::FlowAccessors flowAccessors( elemManager, solverName );

        ElementBasedAssemblyKernel< NUM_FACES, IP >
        kernel( rankOffset, er, esr, lengthTolerance, faceDofKey, nodeManager, faceManager,
                subRegion, dofNumberAccessor, flowAccessors, fluid, permeability,
                regionFilter, dt, localMatrix, localRhs );
        ElementBasedAssemblyKernel< NUM_FACES, IP >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    } );
  }

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY >
  static void
  launch( arrayView1d< real64 const > const & localResidual,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & facePresDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ElementViewConst< arrayView1d< real64 const > > const & elemVolume,
          real64 const & defaultViscosity,
          real64 * localResidualNorm )
  {

    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > sumScaled( 0.0 );

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
            normalizer += elemVolume[er][esr][ei];
            elemCounter++;
          }
        }
        normalizer /= elemCounter;
        normalizer /= defaultViscosity;

        localIndex const lid = LvArray::integerConversion< localIndex >( facePresDofNumber[iface] - rankOffset );
        real64 const val = localResidual[lid] / normalizer; // to get something dimensionless
        sumScaled += val * val;
      }
    } );

    *localResidualNorm = *localResidualNorm + sumScaled.get();
  }

};

} // namespace singlePhaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
