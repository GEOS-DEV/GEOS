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
 * @file SubstrateCoverageSurfaceArea.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SURFACEAREA_SUBSTRATECOVERAGESURFACEAREA_HPP_
#define GEOS_CONSTITUTIVE_SURFACEAREA_SUBSTRATECOVERAGESURFACEAREA_HPP_

#include "constitutive/surfaceArea/SurfaceAreaBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Update class for the substrate coverage reactive surface area model.
 *
 *   g        = max( 0, ( phi - phi_crit ) / ( phi_0 - phi_crit ) )
 *   A_self   = A_m^0 * ( phi_m / phi_m^0 )^n
 *   A_sub    = phi_m^0 * A_tot * exp( -phi_m / phi_m^cov )
 *   A_m      = A_self * g                    on dissolution ( phi_m < phi_m^0 )
 *   A_m      = ( A_self + A_sub ) * g        on precipitation
 *
 * where phi_m is the volume fraction of the mineral m, phi the porosity, n the mineral exponent
 * (2/3 for spherical grains), phi_m^cov the precipitate volume fraction coating the substrate, and
 * A_tot the specific area of the substrate, e.g. 2/b_0 for the flow between two parallel plates.
 *
 * @note Both areas are proportional to the initial state of the mineral, so a mineral with a zero
 * initial volume fraction keeps a zero reactive surface area.
 */
class SubstrateCoverageSurfaceAreaUpdates : public SurfaceAreaBaseUpdates
{
public:

  SubstrateCoverageSurfaceAreaUpdates( arrayView3d< real64, reactivefluid::USD_SPECIES > const & surfaceArea,
                                       arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialSurfaceArea,
                                       arrayView1d< real64 const > const & characteristicVolumeFraction,
                                       real64 const & criticalPorosity,
                                       real64 const & substrateSpecificSurfaceArea,
                                       real64 const & mineralExponent ):
    SurfaceAreaBaseUpdates( surfaceArea, initialSurfaceArea ),
    m_characteristicVolumeFraction( characteristicVolumeFraction ),
    m_criticalPorosity( criticalPorosity ),
    m_substrateSpecificSurfaceArea( substrateSpecificSurfaceArea ),
    m_mineralExponent( mineralExponent )
  {}

  GEOS_HOST_DEVICE
  virtual void updateFromPorosityAndVolumeFractions( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & porosity,
                                                     real64 const & initialPorosity,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & volumeFractions,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & initialVolumeFractions ) const override final
  {
    real64 const accessiblePorosity = initialPorosity - m_criticalPorosity;

    real64 const accessibilityFactor = ( accessiblePorosity > 0.0 ) ?
                                       LvArray::math::max( porosity - m_criticalPorosity, 0.0 ) / accessiblePorosity : 0.0;

    for( integer r = 0; r < numKineticReactions(); ++r )
    {
      real64 const volumeFraction = LvArray::math::max( volumeFractions[r], 0.0 );
      real64 const initialVolumeFraction = initialVolumeFractions[r];

      // the grain size the self area scales from is only known for a mineral present initially
      real64 const selfArea = ( initialVolumeFraction > 0.0 ) ?
                              m_initialSurfaceArea[k][q][r] * pow( volumeFraction / initialVolumeFraction,
                                                                   m_mineralExponent ) : 0.0;

      if( volumeFraction < initialVolumeFraction ) // dissolution
      {
        m_surfaceArea[k][q][r] = selfArea * accessibilityFactor;
      }
      else // precipitation, the mineral also grows on the substrate it has not yet coated
      {
        real64 const coverageFactor = LvArray::math::exp( -volumeFraction / m_characteristicVolumeFraction[r] );

        real64 const substrateArea = initialVolumeFraction * m_substrateSpecificSurfaceArea * coverageFactor;

        m_surfaceArea[k][q][r] = ( selfArea + substrateArea ) * accessibilityFactor;
      }
    }
  }

private:

  /// Characteristic precipitate volume fraction coating the substrate, for each mineral
  arrayView1d< real64 const > m_characteristicVolumeFraction;

  /// Porosity below which the rock does not react anymore
  real64 const m_criticalPorosity;

  /// Specific surface area of the substrate available to precipitation
  real64 const m_substrateSpecificSurfaceArea;

  /// Exponent applied to the mineral volume fraction ratio
  real64 const m_mineralExponent;
};


/**
 * @brief Model in which the reactive surface area of each mineral combines the area of the mineral
 *        grains with the area of the substrate they precipitate onto, and vanishes as the porosity
 *        approaches a critical value.
 */
class SubstrateCoverageSurfaceArea : public SurfaceAreaBase
{
public:

  SubstrateCoverageSurfaceArea( string const & name, Group * const parent );

  static string catalogName() { return "SubstrateCoverageSurfaceArea"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public SurfaceAreaBase::viewKeyStruct
  {
    static constexpr char const * characteristicVolumeFractionString() { return "characteristicVolumeFraction"; }
    static constexpr char const * criticalPorosityString() { return "criticalPorosity"; }
    static constexpr char const * substrateSpecificSurfaceAreaString() { return "substrateSpecificSurfaceArea"; }
    static constexpr char const * mineralExponentString() { return "mineralExponent"; }
  };

  using KernelWrapper = SubstrateCoverageSurfaceAreaUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_surfaceArea,
                          m_initialSurfaceArea,
                          m_characteristicVolumeFraction,
                          m_criticalPorosity,
                          m_substrateSpecificSurfaceArea,
                          m_mineralExponent );
  }

private:

  virtual void postInputInitialization() override;

  /// Characteristic precipitate volume fraction coating the substrate, one value per mineral
  array1d< real64 > m_characteristicVolumeFraction;

  /// Porosity below which the rock does not react anymore
  real64 m_criticalPorosity;

  /// Specific surface area of the substrate available to precipitation
  real64 m_substrateSpecificSurfaceArea;

  /// Exponent applied to the mineral volume fraction ratio
  real64 m_mineralExponent;
};

}/* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_SURFACEAREA_SUBSTRATECOVERAGESURFACEAREA_HPP_
