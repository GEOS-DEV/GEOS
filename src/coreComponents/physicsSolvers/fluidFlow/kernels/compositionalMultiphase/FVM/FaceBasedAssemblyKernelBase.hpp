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
 * @file FaceBasedAssemblyKernelBase.hpp
 */

#ifndef GEOSX_FACEBASEDASSEMBLYKERNELBASE_HPP
#define GEOSX_FACEBASEDASSEMBLYKERNELBASE_HPP

#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/fields/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

using namespace constitutive;

/**
 * @brief Base class for FaceBasedAssemblyKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of components/dofs).
 */
class FaceBasedAssemblyKernelBase
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = geos::ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = geos::ElementRegionManager::ElementViewAccessor< geos::arrayView1d< geos::globalIndex const > >;

  using CompFlowAccessors =
    geos::StencilAccessors< geos::fields::ghostRank,
                            geos::fields::flow::gravityCoefficient,
                            geos::fields::flow::pressure,
                            geos::fields::flow::dGlobalCompFraction_dGlobalCompDensity,
                            geos::fields::flow::dPhaseVolumeFraction,
                            geos::fields::flow::phaseMobility,
                            geos::fields::flow::dPhaseMobility >;
  using MultiFluidAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::MultiFluidBase,
                                    geos::fields::multifluid::phaseMassDensity,
                                    geos::fields::multifluid::dPhaseMassDensity,
                                    geos::fields::multifluid::phaseCompFraction,
                                    geos::fields::multifluid::dPhaseCompFraction >;

  using CapPressureAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::CapillaryPressureBase,
                                    geos::fields::cappres::phaseCapPressure,
                                    geos::fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using PermeabilityAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::PermeabilityBase,
                                    geos::fields::permeability::permeability,
                                    geos::fields::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] compFlowAccessors accessor for wrappers registered by the solver
   * @param[in] multiFluidAccessors accessor for wrappers registered by the multifluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( geos::integer const numPhases,
                               geos::globalIndex const rankOffset,
                               geos::integer const hasCapPressure,
                               DofNumberAccessor const & dofNumberAccessor,
                               CompFlowAccessors const & compFlowAccessors,
                               MultiFluidAccessors const & multiFluidAccessors,
                               CapPressureAccessors const & capPressureAccessors,
                               PermeabilityAccessors const & permeabilityAccessors,
                               geos::real64 const & dt,
                               geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                               geos::arrayView1d< geos::real64 > const & localRhs );

protected:

  /// Number of fluid phases
  geos::integer const m_numPhases;

  /// Offset for my MPI rank
  geos::globalIndex const m_rankOffset;

  /// Flag to specify whether capillary pressure is used or not
  geos::integer const m_hasCapPressure;

  /// Time step size
  geos::real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< geos::arrayView1d< geos::globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< geos::arrayView3d< geos::real64 const > > const m_permeability;
  ElementViewConst< geos::arrayView3d< geos::real64 const > > const m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< geos::arrayView1d< geos::integer const > > const m_ghostRank;
  ElementViewConst< geos::arrayView1d< geos::real64 const > > const m_gravCoef;

  // Primary and secondary variables

  /// Views on pressure
  ElementViewConst< geos::arrayView1d< geos::real64 const > > const m_pres;

  /// Views on derivatives of phase volume fractions and comp fractions
  ElementViewConst< arrayView3d< geos::real64 const, geos::compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac;

  /// Views on phase mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob;

  /// Views on phase mass densities
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseMassDens;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dPhaseMassDens;

  /// Views on phase component fractions
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac;
  ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac;

  /// Views on phase capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
};

}

}

#endif //GEOSX_FACEBASEDASSEMBLYKERNELBASE_HPP
