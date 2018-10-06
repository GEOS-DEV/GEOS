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
  * @file CompositionalMultiphaseFluid.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/Fluid/MultiFluidBase.hpp"

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

class CompositionalMultiphaseFluid : public MultiFluidBase
{
public:

  CompositionalMultiphaseFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~CompositionalMultiphaseFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  static std::string CatalogName() { return dataRepository::keys::compositionalMultiphaseFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override final {}

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

  virtual void StateUpdatePointMultiphaseFluid(real64 const & pres,
                                               real64 const & temp,
                                               real64 const * composition,
                                               localIndex const k,
                                               localIndex const q) override;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr auto equationsOfStateString = "equationsOfState";

    static constexpr auto componentCriticalPressureString    = "componentCriticalPressure";
    static constexpr auto componentCriticalTemperatureString = "componentCriticalTemperature";
    static constexpr auto componentAcentricFactorString      = "componentAcentricFactor";
    static constexpr auto componentMolarWeightString         = "componentMolarWeight";
    static constexpr auto componentVolumeShiftString         = "componentVolumeShift";
    static constexpr auto componentBinaryCoeffString         = "componentBinaryCoeff";

    static constexpr auto phaseMoleFractionString                      = "phaseMoleFraction";                      // xi_p
    static constexpr auto dPhaseMoleFraction_dPressureString           = "dPhaseMoleFraction_dPressure";           // dXi_p/dP
    static constexpr auto dPhaseMoleFraction_dTemperatureString        = "dPhaseMoleFraction_dTemperature";        // dXi_p/dT
    static constexpr auto dPhaseMoleFraction_dGlobalCompFractionString = "dPhaseMoleFraction_dGlobalCompFraction"; // dXi_p/dz_c
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey equationsOfState = { equationsOfStateString };

    ViewKey componentCriticalPressure    = { componentCriticalPressureString };
    ViewKey componentCriticalTemperature = { componentCriticalTemperatureString };
    ViewKey componentAcentricFactor      = { componentAcentricFactorString };
    ViewKey componentMolarWeight         = { componentMolarWeightString };
    ViewKey componentVolumeShift         = { componentVolumeShiftString };
    ViewKey componentBinaryCoeff         = { componentBinaryCoeffString };

  } viewKeys;

private:

  void createFluid();

  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // standard EOS component input
  array1d<real64> m_componentCriticalPressure;
  array1d<real64> m_componentCriticalTemperature;
  array1d<real64> m_componentAcentricFactor;
  array1d<real64> m_componentMolarWeight;
  array1d<real64> m_componentVolumeShift;
  array2d<real64> m_componentBinaryCoeff;

  // PVTPackage fluid object
  PVTPackage::CompositionalMultiphaseSystem * m_fluid;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
