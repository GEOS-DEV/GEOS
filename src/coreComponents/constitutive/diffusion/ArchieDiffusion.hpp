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
 * @file ArchieDiffusion.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DIFFUSION_ARCHIEDIFFUSION_HPP_
#define GEOS_CONSTITUTIVE_DIFFUSION_ARCHIEDIFFUSION_HPP_

#include "constitutive/diffusion/DiffusionBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Kernel-callable update class for the Archie diffusion law with a damage-dependent cementation exponent.
 * The porosity-weighted effective diffusivity is phi * D = D_w * phi^m(d), with m(d) = m0 - ( m0 - mc ) * d going from
 * the matrix exponent m0 to the fracture exponent mc. The flux kernel applies the porosity weight, so the stored
 * (isotropic) diffusivity is D = D_w * phi^( m(d) - 1 ).
 */
class ArchieDiffusionUpdate : public DiffusionBaseUpdate
{
public:

  ArchieDiffusionUpdate( arrayView3d< real64 > const & diffusivity,
                         arrayView3d< real64 > const & dDiffusivity_dTemperature,
                         real64 const fluidDiffusivity,
                         real64 const matrixCementationExponent,
                         real64 const fractureCementationExponent )
    : DiffusionBaseUpdate( diffusivity, dDiffusivity_dTemperature ),
    m_fluidDiffusivity( fluidDiffusivity ),
    m_matrixCementationExponent( matrixCementationExponent ),
    m_fractureCementationExponent( fractureCementationExponent )
  {}

  /**
   * @brief Update the diffusivity of element k from its cell-averaged damage and porosity.
   * @param k element index
   * @param damage cell-averaged damage in [0, 1]
   * @param porosity porosity of element k
   */
  GEOS_HOST_DEVICE
  void updateFromDamageAndPorosity( localIndex const k, real64 const & damage, real64 const & porosity ) const
  {
    real64 const m = m_matrixCementationExponent - ( m_matrixCementationExponent - m_fractureCementationExponent ) * damage;
    real64 const diffusivity = m_fluidDiffusivity * pow( fmax( porosity, 0.0 ), m - 1.0 );
    for( localIndex dim = 0; dim < 3; ++dim )
    {
      m_diffusivity[k][0][dim] = diffusivity;
    }
  }

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const & temperature ) const override
  { GEOS_UNUSED_VAR( k, q, temperature ); }

private:

  real64 m_fluidDiffusivity;
  real64 m_matrixCementationExponent;
  real64 m_fractureCementationExponent;

};


/**
 * @brief Archie diffusion constitutive model with a damage-dependent cementation exponent.
 * The coupled solid updates the diffusivity from the cell-averaged damage and the porosity;
 * without a damage solid (d = 0) it reduces to the Archie law with the matrix exponent m0.
 */
class ArchieDiffusion : public DiffusionBase
{
public:

  ArchieDiffusion( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "ArchieDiffusion"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ArchieDiffusionUpdate;

  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_diffusivity,
                          m_dDiffusivity_dTemperature,
                          m_fluidDiffusivity,
                          m_matrixCementationExponent,
                          m_fractureCementationExponent );
  }

  struct viewKeyStruct : public DiffusionBase::viewKeyStruct
  {
    static constexpr char const * fluidDiffusivityString() { return "fluidDiffusivity"; }
    static constexpr char const * matrixCementationExponentString() { return "matrixCementationExponent"; }
    static constexpr char const * fractureCementationExponentString() { return "fractureCementationExponent"; }
  };

protected:

  virtual void postInputInitialization() override;

private:

  /// Diffusivity D_w of the pore fluid [m^2/s]
  real64 m_fluidDiffusivity;

  /// Archie cementation exponent of the intact matrix
  real64 m_matrixCementationExponent;

  /// Archie cementation exponent of the fully damaged material (open fracture)
  real64 m_fractureCementationExponent;

};

} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_DIFFUSION_ARCHIEDIFFUSION_HPP_
