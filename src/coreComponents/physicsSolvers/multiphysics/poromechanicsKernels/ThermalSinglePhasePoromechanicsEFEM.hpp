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
 * @file ThermalSinglePhasePoromechanicsEFEM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_HPP_

#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp"

namespace geos
{

namespace thermoPoromechanicsEFEMKernels
{
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermalSinglePhasePoromechanicsEFEM :
  public poromechanicsEFEMKernels::SinglePhasePoromechanicsEFEM< SUBREGION_TYPE,
                                                                 CONSTITUTIVE_TYPE,
                                                                 FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = poromechanicsEFEMKernels::SinglePhasePoromechanicsEFEM< SUBREGION_TYPE,
                                                                       CONSTITUTIVE_TYPE,
                                                                       FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  /// Compile time value for the number of gotquadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  ThermalSinglePhasePoromechanicsEFEM( NodeManager const & nodeManager,
                                       EdgeManager const & edgeManager,
                                       FaceManager const & faceManager,
                                       localIndex const targetRegionIndex,
                                       SUBREGION_TYPE const & elementSubRegion,
                                       FE_TYPE const & finiteElementSpace,
                                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                                       EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
                                       arrayView1d< globalIndex const > const dispDofNumber,
                                       arrayView1d< globalIndex const > const jumpDofNumber,
                                       string const inputFlowDofKey,
                                       globalIndex const rankOffset,
                                       CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                       arrayView1d< real64 > const inputRhs,
                                       real64 const (&inputGravityVector)[3],
                                       string const fluidModelKey );

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geos::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;


    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables()
    {}
  };
  //*****************************************************************************

  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );
  //END_kernelLauncher

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

private:

};


using ThermalSinglePhasePoromechanicsEFEMKernelFactory = finiteElement::KernelFactory< ThermalSinglePhasePoromechanicsEFEM,
                                                                                       EmbeddedSurfaceSubRegion const &,
                                                                                       arrayView1d< globalIndex const > const,
                                                                                       arrayView1d< globalIndex const > const,
                                                                                       string const,
                                                                                       globalIndex const,
                                                                                       CRSMatrixView< real64, globalIndex const > const,
                                                                                       arrayView1d< real64 > const,
                                                                                       real64 const (&)[3],
                                                                                       string const >;

} // namespace thermoPoromechanicsEFEMKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_HPP_
