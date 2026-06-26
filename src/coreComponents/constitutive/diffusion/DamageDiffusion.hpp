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
 * @file DamageDiffusion.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DIFFUSION_DAMAGEDIFFUSION_HPP_
#define GEOS_CONSTITUTIVE_DIFFUSION_DAMAGEDIFFUSION_HPP_

#include "constitutive/diffusion/DiffusionBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Kernel-callable update class for damage-dependent diffusion.
 * Diffusivity scales with cell-averaged damage as D_dim = bulkDiffusivity * exp( alpha_dim * damage ),
 * where alpha_dim is a per-direction (x, y, z) exponent, allowing damage-driven enhancement
 * to be anisotropic (e.g. channel transport along the intended crack-growth direction only).
 */
class DamageDiffusionUpdate : public DiffusionBaseUpdate
{
public:

  DamageDiffusionUpdate( arrayView3d< real64 > const & diffusivity,
                         arrayView3d< real64 > const & dDiffusivity_dTemperature,
                         arrayView1d< real64 const > const & bulkDiffusivity,
                         R1Tensor const & damageDependenceConstant )
    : DiffusionBaseUpdate( diffusivity, dDiffusivity_dTemperature ),
    m_bulkDiffusivity( bulkDiffusivity ),
    m_damageDependenceConstant( damageDependenceConstant )
  {}

  /**
   * @brief Update the three (possibly anisotropic) diffusivity components for element k from cell-averaged damage.
   * @param k element index
   * @param damage cell-averaged damage in [0, 1]
   */
  GEOS_HOST_DEVICE
  void updateDamageDiffusivity( localIndex const k, real64 const & damage ) const
  {
    for( localIndex dim = 0; dim < 3; ++dim )
    {
      m_diffusivity[k][0][dim] = m_bulkDiffusivity[k] * LvArray::math::exp( m_damageDependenceConstant[dim] * damage );
    }
  }

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & temperature ) const override
  { GEOS_UNUSED_VAR( k, q, temperature ); }

private:

  arrayView1d< real64 const > m_bulkDiffusivity;
  R1Tensor m_damageDependenceConstant;

};


/**
 * @brief Damage-dependent diffusion constitutive model.
 * Intended for use in EigenstrainReactiveSolid with DamageElasticIsotropic.
 */
class DamageDiffusion : public DiffusionBase
{
public:

  DamageDiffusion( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "DamageDiffusion"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = DamageDiffusionUpdate;

  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_diffusivity,
                          m_dDiffusivity_dTemperature,
                          m_bulkDiffusivity,
                          m_damageDependenceConstant );
  }

  struct viewKeyStruct : public DiffusionBase::viewKeyStruct
  {
    static constexpr char const * defaultBulkDiffusivityString() { return "defaultBulkDiffusivity"; }
    static constexpr char const * damageDependenceConstantString() { return "damageDependenceConstant"; }
  };

private:

  /// Default (uniform) diffusivity of the intact bulk material [m^2/s], read from XML
  real64 m_defaultBulkDiffusivity;

  /// Per-element diffusivity of the intact (undamaged) bulk material [m^2/s]
  array1d< real64 > m_bulkDiffusivity;

  /// Per-direction (x, y, z) exponents controlling sensitivity of diffusivity to damage
  R1Tensor m_damageDependenceConstant;

};

} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_DIFFUSION_DAMAGEDIFFUSION_HPP_
