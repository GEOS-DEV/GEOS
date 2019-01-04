/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::brooksCoreyRelativePermeability; }

  virtual string GetCatalogName() override { return CatalogName(); }


  // RelPerm-specific interface

  virtual void StateUpdatePointRelPerm( arraySlice1d<real64 const> const & phaseVolFraction,
                                        localIndex const k,
                                        localIndex const q ) override;

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

} // namespace constitutive

} // namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
