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
  * @file CompositionalMultiphaseFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/Fluid/MultiFluidPVTPackageWrapper.hpp"

namespace PVTPackage
{
class CompositionalMultiphaseSystem;
}

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const compositionalMultiphaseFluid = "CompositionalMultiphaseFluid";
}
}

namespace constitutive
{

class CompositionalMultiphaseFluid : public MultiFluidPVTPackageWrapper
{
public:

  CompositionalMultiphaseFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~CompositionalMultiphaseFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::compositionalMultiphaseFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }


  struct viewKeyStruct : MultiFluidPVTPackageWrapper::viewKeyStruct
  {
    static constexpr auto equationsOfStateString = "equationsOfState";

    static constexpr auto componentCriticalPressureString    = "componentCriticalPressure";
    static constexpr auto componentCriticalTemperatureString = "componentCriticalTemperature";
    static constexpr auto componentAcentricFactorString      = "componentAcentricFactor";
    static constexpr auto componentVolumeShiftString         = "componentVolumeShift";
    static constexpr auto componentBinaryCoeffString         = "componentBinaryCoeff";
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey equationsOfState = { equationsOfStateString };

    ViewKey componentCriticalPressure    = { componentCriticalPressureString };
    ViewKey componentCriticalTemperature = { componentCriticalTemperatureString };
    ViewKey componentAcentricFactor      = { componentAcentricFactorString };
    ViewKey componentVolumeShift         = { componentVolumeShiftString };
    ViewKey componentBinaryCoeff         = { componentBinaryCoeffString };

  } viewKeysCompositionalMultiphaseFluid;

protected:
  virtual void PostProcessInput() override;

private:

  void createFluid() override;

  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // standard EOS component input
  array1d<real64> m_componentCriticalPressure;
  array1d<real64> m_componentCriticalTemperature;
  array1d<real64> m_componentAcentricFactor;
  array1d<real64> m_componentVolumeShift;
  array2d<real64> m_componentBinaryCoeff;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
