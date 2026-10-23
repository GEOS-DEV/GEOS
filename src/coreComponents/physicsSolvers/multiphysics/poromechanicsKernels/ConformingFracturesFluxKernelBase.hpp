/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ConformingFracturesFluxKernelBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_CONFORMINGFRACTURESFLUXKERNELBASE_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_CONFORMINGFRACTURESFLUXKERNELBASE_HPP

#include "constitutive/permeability/PermeabilityFields.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
{

namespace poromechanicsConformingFracturesKernels
{

/**
 * @class ConformingFracturesFluxKernelBase
 * @brief Aperture-dependent part of the flux kernels on conforming fractures, shared by the single-phase
 *        and multiphase kernels: fracture permeability / hydraulic aperture views, dR_dAper matrix, and the
 *        aperture derivatives of the transmissibilities (incl. at fracture junctions).
 * @tparam FLUX_KERNEL_BASE the flow flux kernel on a SurfaceElementStencilWrapper
 *         (e.g. singlePhaseFVMKernels::FluxComputeKernel, isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel)
 */
template< typename FLUX_KERNEL_BASE >
class ConformingFracturesFluxKernelBase : public FLUX_KERNEL_BASE
{
public:

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using FracturePermeabilityAccessors = StencilMaterialAccessors< constitutive::PermeabilityBase,
                                                                  fields::permeability::dPerm_dDispJump >;
  using FractureApertureAccessors = StencilAccessors< fields::flow::hydraulicAperture,
                                                      fields::flow::dHydraulicAperture_dNormalJump >;

  using Base = FLUX_KERNEL_BASE;
  using Base::numDof;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::m_dt;
  using Base::m_dofNumber;
  using Base::m_permeability;
  using Base::m_dPerm_dPres;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor
   * @param[in] fracturePermeabilityAccessors accessors to the derivative of the fracture permeability wrt the jump
   * @param[in] fractureApertureAccessors accessors to the hydraulic aperture and its derivative wrt the normal jump
   * @param[in] dR_dAper matrix of the derivatives of the flow residual wrt the normal jump
   * @param[in] baseArgs the arguments of the constructor of FLUX_KERNEL_BASE
   */
  template< typename ... BASE_ARGS >
  ConformingFracturesFluxKernelBase( FracturePermeabilityAccessors const & fracturePermeabilityAccessors,
                                     FractureApertureAccessors const & fractureApertureAccessors,
                                     CRSMatrixView< real64, localIndex const > const & dR_dAper,
                                     BASE_ARGS && ... baseArgs )
    : Base( std::forward< BASE_ARGS >( baseArgs )... ),
    m_dR_dAper( dR_dAper ),
    m_dPerm_dDispJump( fracturePermeabilityAccessors.get( fields::permeability::dPerm_dDispJump {} ) ),
    m_hydraulicAperture( fractureApertureAccessors.get( fields::flow::hydraulicAperture {} ) ),
    m_dHydraulicAperture_dNormalJump( fractureApertureAccessors.get( fields::flow::dHydraulicAperture_dNormalJump {} ) )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack, common to the conforming fracture flux kernels
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      localColIndices( numElems )
    {}

    /// Local (subregion) indices of the elements of the stencil, the columns of dR_dAper
    stackArray1d< localIndex, maxNumElems > localColIndices;

    /// t_j * d(transmissibility)/d(t_j) of each connection, for every half-transmissibility t_j = weight_j * permeability_j
    /// of the stencil (unsigned: derivative of entry [0] of the transmissibility). At a fracture junction the
    /// transmissibility of a connection also depends on the elements it does not connect.
    real64 dTrans_dLogHalfTrans[maxNumConns][maxNumElems]{};

    /// d(ln t_j)/d(normalJump_j) for every element j of the stencil: permeability part plus stencil-weight part
    /// (the stencil weight is linear in the hydraulic aperture)
    real64 dLogHalfTrans_dNormalJump[maxNumElems]{};

    /// d(ln hydraulicAperture_j)/d(normalJump_j) for every element j of the stencil (stencil-weight part of
    /// dLogHalfTrans_dNormalJump), kept separately for the derived kernels (e.g. thermal conduction)
    real64 dLogAperture_dNormalJump[maxNumElems]{};
  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[i * numDof + jdof] = offset + jdof;
      }
      stack.localColIndices[ i ] = m_sei( iconn, i );
    }
  }

  /**
   * @brief Compute the transmissibilities, their pressure derivatives and their aperture derivatives
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeTransmissibilityAndApertureDerivatives( localIndex const iconn,
                                                      StackVariables & stack ) const
  {
    // Transmissibilities, their pressure derivatives and their aperture derivatives (wrt the half-transmissibility
    // of every element of the stencil), in a single stencil traversal
    m_stencilWrapper.computeWeightsAndLogDerivatives( iconn,
                                                      m_permeability,
                                                      m_dPerm_dPres,
                                                      stack.transmissibility,
                                                      stack.dTrans_dPres,
                                                      stack.dTrans_dLogHalfTrans );

    // t_j = weight_j(hydraulicAperture_j) * perm_j(hydraulicAperture_j), the weight being linear in the hydraulic aperture
    for( localIndex j = 0; j < stack.stencilSize; ++j )
    {
      localIndex const er  = m_seri[iconn][j];
      localIndex const esr = m_sesri[iconn][j];
      localIndex const ei  = m_sei[iconn][j];

      real64 const hydraulicAperture = m_hydraulicAperture[er][esr][ei];
      stack.dLogAperture_dNormalJump[j] = ( hydraulicAperture > 0.0 )
        ? m_dHydraulicAperture_dNormalJump[er][esr][ei] / hydraulicAperture
        : 0.0;

      real64 dLogHalfTrans_dNormalJump = stack.dLogAperture_dNormalJump[j];
      real64 const perm = m_permeability[er][esr][ei][0][0];
      if( perm > 0.0 )
      {
        dLogHalfTrans_dNormalJump += m_dPerm_dDispJump[er][esr][ei][0][0][0] / perm;
      }
      stack.dLogHalfTrans_dNormalJump[j] = dLogHalfTrans_dNormalJump;
    }
  }

  /**
   * @brief Accumulate the aperture derivative of the flux of one connection, for every element of the stencil
   *        (at a junction, also those this connection does not link)
   * @tparam DFLUX_DAPERTURE the type of the [element][stencil element] slice receiving the derivatives
   * @param[in] stack the stack variables
   * @param[in] k the two elements of the connection
   * @param[in] connectionIndex the connection index
   * @param[in] dFlux_dTrans the derivative of the flux wrt the transmissibility
   * @param[inout] dFlux_dAperture the accumulated derivatives (+ on k[0], - on k[1])
   * @param[out] dFlux_dAper the derivatives of the flux of this connection only (multiplied by dt)
   */
  template< typename DFLUX_DAPERTURE >
  GEOS_HOST_DEVICE
  void accumulateFluxApertureDerivative( StackVariables const & stack,
                                         localIndex const (&k)[2],
                                         localIndex const connectionIndex,
                                         real64 const dFlux_dTrans,
                                         DFLUX_DAPERTURE && dFlux_dAperture,
                                         real64 ( & dFlux_dAper )[maxNumElems] ) const
  {
    for( localIndex j = 0; j < stack.stencilSize; ++j )
    {
      dFlux_dAper[j] = m_dt * dFlux_dTrans * stack.dTrans_dLogHalfTrans[connectionIndex][j] * stack.dLogHalfTrans_dNormalJump[j];
      dFlux_dAperture[k[0]][j] += dFlux_dAper[j];
      dFlux_dAperture[k[1]][j] -= dFlux_dAper[j];
    }
  }

protected:

  /// Derivatives of the flow residual wrt the normal jump
  CRSMatrixView< real64, localIndex const > m_dR_dAper;

  /// View on the derivative of the fracture permeability wrt the displacement jump
  ElementViewConst< arrayView4d< real64 const > > const m_dPerm_dDispJump;

  /// Views on the hydraulic aperture and its derivative wrt the normal jump
  ElementViewConst< arrayView1d< real64 const > > const m_hydraulicAperture;
  ElementViewConst< arrayView1d< real64 const > > const m_dHydraulicAperture_dNormalJump;
};

} // namespace poromechanicsConformingFracturesKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_CONFORMINGFRACTURESFLUXKERNELBASE_HPP
