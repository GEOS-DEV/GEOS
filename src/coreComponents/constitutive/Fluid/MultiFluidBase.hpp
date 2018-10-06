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

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{

namespace constitutive
{

class MultiFluidBase : public ConstitutiveBase
{
public:

  // define a limit on number of components (used for optimization purposes)
  static constexpr localIndex MAX_NUM_COMPONENTS = 32;

  MultiFluidBase( std::string const & name, ManagedGroup * const parent );

  virtual ~MultiFluidBase() override;

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** MultiphaseFluid-specific interface

  localIndex numFluidComponents() const;

  localIndex numFluidPhases() const;

  bool getMassFractionFlag() const;

  void setMassFractionFlag( bool flag );

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phasesString           = "phases";
    static constexpr auto componentNamesString   = "componentNames";

    static constexpr auto phaseFractionString                            = "phaseFraction";                      // xi_p
    static constexpr auto dPhaseFraction_dPressureString                 = "dPhaseFraction_dPressure";           // dXi_p/dP
    static constexpr auto dPhaseFraction_dTemperatureString              = "dPhaseFraction_dTemperature";        // dXi_p/dT
    static constexpr auto dPhaseFraction_dGlobalCompFractionString       = "dPhaseFraction_dGlobalCompFraction"; // dXi_p/dz_c

    static constexpr auto phaseVolumeFractionString                      = "phaseVolumeFraction";                      // S_p
    static constexpr auto dPhaseVolumeFraction_dPressureString           = "dPhaseVolumeFraction_dPressure";           // dS_p/dP
    static constexpr auto dPhaseVolumeFraction_dTemperatureString        = "dPhaseVolumeFraction_dTemperature";        // dS_p/dT
    static constexpr auto dPhaseVolumeFraction_dGlobalCompFractionString = "dPhaseVolumeFraction_dGlobalCompFraction"; // dS_p/dz_c

    static constexpr auto phaseDensityString                             = "phaseDensity";                       // rho_p
    static constexpr auto dPhaseDensity_dPressureString                  = "dPhaseDensity_dPressure";           // dRho_p/dP
    static constexpr auto dPhaseDensity_dTemperatureString               = "dPhaseDensity_dTemperature";        // dRho_p/dT
    static constexpr auto dPhaseDensity_dGlobalCompFractionString        = "dPhaseDensity_dGlobalCompFraction"; // dRho_p/dz_c

    static constexpr auto phaseComponentFractionString                   = "phaseComponentFraction";                 // x_cp
    static constexpr auto dPhaseCompFraction_dPressureString             = "dPhaseCompFraction_dPressure";           // dx_cp/dP
    static constexpr auto dPhaseCompFraction_dTemperatureString          = "dPhaseCompFraction_dTemperature";        // dx_cp/dT
    static constexpr auto dPhaseCompFraction_dGlobalCompFractionString   = "dPhaseCompFraction_dGlobalCompFraction"; // dx_cp/dz_c

    using ViewKey = dataRepository::ViewKey;

    ViewKey phases           = { phasesString };
    ViewKey componentNames   = { componentNamesString };

    ViewKey phaseFraction                                = { phaseFractionString };                        // xi_p
    ViewKey dPhaseFraction_dPressure                     = { dPhaseFraction_dPressureString };             // dXi_p/dP
    ViewKey dPhaseFraction_dTemperature                  = { dPhaseFraction_dTemperatureString };          // dXi_p/dT
    ViewKey dPhaseFraction_dGlobalCompMoleFraction       = { dPhaseFraction_dGlobalCompFractionString };   // dXi_p/dz_c

    ViewKey phaseVolumeFraction                          = { phaseVolumeFractionString };                      // S_p
    ViewKey dPhaseVolumeFraction_dPressure               = { dPhaseVolumeFraction_dPressureString };           // dS_p/dP
    ViewKey dPhaseVolumeFraction_dTemperature            = { dPhaseVolumeFraction_dTemperatureString };        // dS_p/dT
    ViewKey dPhaseVolumeFraction_dGlobalCompMoleFraction = { dPhaseVolumeFraction_dGlobalCompFractionString }; // dS_p/dz_c

    ViewKey phaseDensity                                 = { phaseDensityString };                             // rho_p
    ViewKey dPhaseDensity_dPressure                      = { dPhaseDensity_dPressureString };                  // dRho_p/dP
    ViewKey dPhaseDensity_dTemperature                   = { dPhaseDensity_dTemperatureString };               // dRho_p/dT
    ViewKey dPhaseDensity_dGlobalCompMoleFraction        = { dPhaseDensity_dGlobalCompFractionString };        // dRho_p/dz_c

    ViewKey phaseCompFraction                            = { phaseComponentFractionString };                   // x_cp
    ViewKey dPhaseCompFraction_dPressure                 = { dPhaseCompFraction_dPressureString };             // dx_cp/dP
    ViewKey dPhaseCompFraction_dTemperature              = { dPhaseCompFraction_dTemperatureString };          // dx_cp/dT
    ViewKey dPhaseCompFraction_dGlobalCompFraction       = { dPhaseCompFraction_dGlobalCompFractionString };   // dx_cp/dz_c

  } viewKeysMultiFluidBase;

protected:

  // flag indicating whether input/output component fractions are treated as mass fractions
  bool m_useMassFractions;

  // general fluid composition information
  string_array m_phases;
  string_array m_componentNames;

  array3d<real64> m_phaseFraction;
  array3d<real64> m_dPhaseFraction_dPressure;
  array3d<real64> m_dPhaseFraction_dTemperature;
  array4d<real64> m_dPhaseFraction_dGlobalCompFraction;

  array3d<real64> m_phaseVolumeFraction;
  array3d<real64> m_dPhaseVolumeFraction_dPressure;
  array3d<real64> m_dPhaseVolumeFraction_dTemperature;
  array4d<real64> m_dPhaseVolumeFraction_dGlobalCompFraction;

  array3d<real64> m_phaseDensity;
  array3d<real64> m_dPhaseDensity_dPressure;
  array3d<real64> m_dPhaseDensity_dTemperature;
  array4d<real64> m_dPhaseDensity_dGlobalCompFraction;

  array4d<real64> m_phaseCompFraction;
  array4d<real64> m_dPhaseCompFraction_dPressure;
  array4d<real64> m_dPhaseCompFraction_dTemperature;
  array5d<real64> m_dPhaseCompFraction_dGlobalCompFraction;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_MULTIFLUIDBASE_HPP
