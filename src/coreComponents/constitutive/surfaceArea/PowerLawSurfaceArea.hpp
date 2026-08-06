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
 * @file PowerLawSurfaceArea.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SURFACEAREA_POWERLAWSURFACEAREA_HPP_
#define GEOS_CONSTITUTIVE_SURFACEAREA_POWERLAWSURFACEAREA_HPP_

#include "constitutive/surfaceArea/SurfaceAreaBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Update class for the power law reactive surface area model.
 *
 * The specific surface area of the mineral m is scaled from its initial value with the amount of
 * mineral left and with the pore space available to the fluid. On dissolution
 *
 *   A_m = A_m^0 * ( phi_m / phi_m^0 )^n * ( phi / phi_0 )^n'
 *
 * whereas on precipitation the mineral term drops out
 *
 *   A_m = A_m^0 * ( phi / phi_0 )^n'
 *
 * where phi_m is the volume fraction of the mineral m, phi the porosity, and n, n' the mineral
 * and porosity exponents. Both are 2/3 for the classical spherical grain assumption, which keeps
 * the number of reacting mineral grains in a representative elementary volume constant.
 *
 * @note The dissolution form only applies to a primary mineral, with a nonzero initial volume
 * fraction. A mineral appearing from phi_m^0 = 0 only ever takes the precipitation form.
 */
class PowerLawSurfaceAreaUpdates : public SurfaceAreaBaseUpdates
{
public:

  PowerLawSurfaceAreaUpdates( arrayView3d< real64, reactivefluid::USD_SPECIES > const & surfaceArea,
                              arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialSurfaceArea,
                              real64 const & mineralExponent,
                              real64 const & porosityExponent ):
    SurfaceAreaBaseUpdates( surfaceArea, initialSurfaceArea ),
    m_mineralExponent( mineralExponent ),
    m_porosityExponent( porosityExponent )
  {}

  GEOS_HOST_DEVICE
  virtual void updateFromPorosityAndVolumeFractions( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & porosity,
                                                     real64 const & initialPorosity,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & volumeFractions,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & initialVolumeFractions ) const override final
  {
    real64 const porosityRatio = ( initialPorosity > 0.0 ) ?
                                 LvArray::math::max( porosity, 0.0 ) / initialPorosity : 1.0;

    real64 const porosityFactor = pow( porosityRatio, m_porosityExponent );

    for( integer r = 0; r < numKineticReactions(); ++r )
    {
      // the grain size the mineral term scales from is only known for a dissolving primary mineral
      bool const isDissolving = volumeFractions[r] < initialVolumeFractions[r];

      real64 const mineralFactor = ( isDissolving && initialVolumeFractions[r] > 0.0 ) ?
                                   pow( LvArray::math::max( volumeFractions[r], 0.0 ) / initialVolumeFractions[r],
                                        m_mineralExponent ) : 1.0;

      m_surfaceArea[k][q][r] = m_initialSurfaceArea[k][q][r] * mineralFactor * porosityFactor;
    }
  }

private:

  /// Exponent applied to the mineral volume fraction ratio
  real64 const m_mineralExponent;

  /// Exponent applied to the porosity ratio
  real64 const m_porosityExponent;
};


/**
 * @brief Model in which the reactive surface area of each mineral follows a power law of the
 *        mineral volume fraction and of the porosity.
 */
class PowerLawSurfaceArea : public SurfaceAreaBase
{
public:

  PowerLawSurfaceArea( string const & name, Group * const parent );

  static string catalogName() { return "PowerLawSurfaceArea"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public SurfaceAreaBase::viewKeyStruct
  {
    static constexpr char const * mineralExponentString() { return "mineralExponent"; }
    static constexpr char const * porosityExponentString() { return "porosityExponent"; }
  };

  using KernelWrapper = PowerLawSurfaceAreaUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_surfaceArea,
                          m_initialSurfaceArea,
                          m_mineralExponent,
                          m_porosityExponent );
  }

private:

  virtual void postInputInitialization() override;

  /// Exponent applied to the mineral volume fraction ratio
  real64 m_mineralExponent;

  /// Exponent applied to the porosity ratio
  real64 m_porosityExponent;
};

}/* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_SURFACEAREA_POWERLAWSURFACEAREA_HPP_
