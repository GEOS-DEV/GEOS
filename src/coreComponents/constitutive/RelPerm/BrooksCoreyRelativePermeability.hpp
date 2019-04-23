/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file BrooksCoreyRelativePermeability.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP

#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const brooksCoreyRelativePermeability = "BrooksCoreyRelativePermeability";
}
}

namespace constitutive
{

class BrooksCoreyRelativePermeability : public RelativePermeabilityBase
{
public:

  BrooksCoreyRelativePermeability( std::string const & name, dataRepository::ManagedGroup * const parent );

  virtual ~BrooksCoreyRelativePermeability() override;

  void DeliverClone( string const & name,
                     ManagedGroup * const parent,
                     std::unique_ptr<ConstitutiveBase> & clone ) const override;

  static std::string CatalogName() { return dataRepository::keys::brooksCoreyRelativePermeability; }

  virtual string GetCatalogName() override { return CatalogName(); }


  // RelPerm-specific interface

  /**
   * @brief Function to update state of a single material point.
   * @param[in] phaseVolFraction input phase volume fraction
   * @param[in] k the first index of the storage arrays (elem index)
   * @param[in] q the secound index of the storage arrays (quadrature index)
   *
   * @note This function performs a point update, but should not be called
   *       within a kernel since it is virtual, and the required data is not
   *       guaranteed to be in the target memory space.
   */
  virtual void PointUpdate( arraySlice1d<real64 const> const & phaseVolFraction,
                            localIndex const k,
                            localIndex const q ) override;

   /**
   * @brief Perform a batch constitutive update (all points).
   * @param[in] phaseVolFraction input phase volume fraction
   */
  virtual void BatchUpdate( arrayView2d<real64 const> const & phaseVolumeFraction ) override;

    /**
   * @brief Computes the phase relative permeabilities using the Brooks method
   * @param NP phase index
   * @param[in] phaseVolumeFraction for all phases
   * @param[out] phaseRelPerm the computed relative permeability value vector for all phases
   * @param[out] dPhaseRelPerm_dPhaseVolFrac the computed partial derivative of the relative wrt to the volume fraction of the phases
   * @param[in] phaseOrder vector of phase orders
   * @param[in] phaseMinVolumeFraction vector of minimum phase volume fractions
   * @param[in] phaseRelPermExponentInv vector of exponents used in the computation of the relative permeabilities
   * @param[in] phaseRelPermMaxValue vector of permeability curve end-point values 
   * @param[in] volFracScale scaling factor to apply to the entire relative permeability curve
   * @return (void)
   *
   * This function computes the relative permeabilities for all phases based on the Brooks-Corey method
   */
  inline static void Compute( localIndex const NP,
                              arraySlice1d<real64 const> const & phaseVolFraction,
                              arraySlice1d<real64> const & phaseRelPerm,
                              arraySlice2d<real64> const & dPhaseRelPerm_dPhaseVolFrac,
                              arraySlice1d<real64 const> const & phaseMinVolumeFraction,
                              arraySlice1d<real64 const> const & phaseRelPermExponent,
                              arraySlice1d<real64 const> const & phaseRelPermMaxValue,
                              real64 const & satScale );

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString = "phaseMinVolumeFraction";
    static constexpr auto phaseRelPermExponentString   = "phaseRelPermExponent";
    static constexpr auto phaseRelPermMaxValueString   = "phaseRelPermMaxValue";

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString };
    ViewKey phaseRelPermExponent   = { phaseRelPermExponentString };
    ViewKey phaseRelPermMaxValue   = { phaseRelPermMaxValueString };

  } vieKeysBrooksCoreyRelativePermeability;

protected:
  virtual void PostProcessInput() override;

  array1d<real64> m_phaseMinVolumeFraction;
  array1d<real64> m_phaseRelPermExponent;
  array1d<real64> m_phaseRelPermMaxValue;

  real64 m_satScale;


};


inline void
BrooksCoreyRelativePermeability::Compute( localIndex const NP,
                                          arraySlice1d<real64 const> const & phaseVolFraction,
                                          arraySlice1d<real64> const & phaseRelPerm,
                                          arraySlice2d<real64> const & dPhaseRelPerm_dPhaseVolFrac,
                                          arraySlice1d<real64 const> const & phaseMinVolumeFraction,
                                          arraySlice1d<real64 const> const & phaseRelPermExponent,
                                          arraySlice1d<real64 const> const & phaseRelPermMaxValue,
                                          real64 const & satScale )
{

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    for (localIndex jp = 0; jp < NP; ++jp)
    {
      dPhaseRelPerm_dPhaseVolFrac[ip][jp] = 0.0;
    }
  }
  real64 const satScaleInv = 1.0 / satScale;

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    real64 const satScaled = (phaseVolFraction[ip] - phaseMinVolumeFraction[ip]) * satScaleInv;
    real64 const exponent  = phaseRelPermExponent[ip];
    real64 const scale     = phaseRelPermMaxValue[ip];

    if (satScaled > 0.0 && satScaled < 1.0)
    {
      // intermediate value
      real64 const v = scale * std::pow( satScaled, exponent - 1.0 );

      phaseRelPerm[ip] = v * satScaled;
      dPhaseRelPerm_dPhaseVolFrac[ip][ip] = v * exponent * satScaleInv;
    }
    else
    {
      phaseRelPerm[ip] = (satScaled < 0.0) ? 0.0 : scale;
    }
  }
}

} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
