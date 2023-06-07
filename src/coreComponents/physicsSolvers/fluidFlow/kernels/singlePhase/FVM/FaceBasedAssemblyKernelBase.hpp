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


namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** FaceBasedAssemblyKernelBase ********************************/

/**
 * @brief Base class for FaceBasedAssemblyKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of dofs).
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

  using SinglePhaseFlowAccessors =
    geos::StencilAccessors< geos::fields::ghostRank,
                            geos::fields::flow::pressure,
                            geos::fields::flow::pressure_n,
                            geos::fields::flow::gravityCoefficient,
                            geos::fields::flow::mobility,
                            geos::fields::flow::dMobility_dPressure >;

  using SinglePhaseFluidAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::SingleFluidBase,
                                    geos::fields::singlefluid::density,
                                    geos::fields::singlefluid::dDensity_dPressure >;

  using SlurryFluidAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::SlurryFluidBase,
                                    geos::fields::singlefluid::density,
                                    geos::fields::singlefluid::dDensity_dPressure >;

  using PermeabilityAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::PermeabilityBase,
                                    geos::fields::permeability::permeability,
                                    geos::fields::permeability::dPerm_dPressure >;

  using ProppantPermeabilityAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::PermeabilityBase,
                                    geos::fields::permeability::permeability,
                                    geos::fields::permeability::dPerm_dPressure,
                                    geos::fields::permeability::dPerm_dDispJump,
                                    geos::fields::permeability::permeabilityMultiplier >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] singleFlowAccessors accessor for wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the singlefluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( geos::globalIndex const rankOffset,
                               DofNumberAccessor const & dofNumberAccessor,
                               SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                               SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                               PermeabilityAccessors const & permeabilityAccessors,
                               geos::real64 const & dt,
                               geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                               geos::arrayView1d< geos::real64 > const & localRhs )
    : m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( geos::fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( geos::fields::permeability::dPerm_dPressure {} ) ),
    m_ghostRank( singlePhaseFlowAccessors.get( geos::fields::ghostRank {} ) ),
    m_gravCoef( singlePhaseFlowAccessors.get( geos::fields::flow::gravityCoefficient {} ) ),
    m_pres( singlePhaseFlowAccessors.get( geos::fields::flow::pressure {} ) ),
    m_mob( singlePhaseFlowAccessors.get( geos::fields::flow::mobility {} ) ),
    m_dMob_dPres( singlePhaseFlowAccessors.get( geos::fields::flow::dMobility_dPressure {} ) ),
    m_dens( singlePhaseFluidAccessors.get( geos::fields::singlefluid::density {} ) ),
    m_dDens_dPres( singlePhaseFluidAccessors.get( geos::fields::singlefluid::dDensity_dPressure {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  { }

protected:

  /// Offset for my MPI rank
  geos::globalIndex const m_rankOffset;

  /// Time step size
  geos::real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< geos::arrayView1d< geos::globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< geos::arrayView3d< geos::real64 const > > m_permeability;
  ElementViewConst< geos::arrayView3d< geos::real64 const > > m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< geos::arrayView1d< geos::integer const > > const m_ghostRank;
  ElementViewConst< geos::arrayView1d< geos::real64 const > > const m_gravCoef;

  // Primary and secondary variables
  /// Views on pressure
  ElementViewConst< geos::arrayView1d< geos::real64 const > > const m_pres;

  /// Views on fluid mobility
  ElementViewConst< geos::arrayView1d< geos::real64 const > > const m_mob;
  ElementViewConst< geos::arrayView1d< geos::real64 const > > const m_dMob_dPres;

  /// Views on fluid density
  ElementViewConst< geos::arrayView2d< geos::real64 const > > const m_dens;
  ElementViewConst< geos::arrayView2d< geos::real64 const > > const m_dDens_dPres;

  // Residual and jacobian

  /// View on the local CRS matrix
  geos::CRSMatrixView< geos::real64, geos::globalIndex const > const m_localMatrix;
  /// View on the local RHS
  geos::arrayView1d< geos::real64 > const m_localRhs;
};

}

}

#endif //GEOSX_FACEBASEDASSEMBLYKERNELBASE_HPP
