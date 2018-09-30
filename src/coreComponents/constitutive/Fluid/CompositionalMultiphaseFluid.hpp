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

#include "constitutive/ConstitutiveBase.hpp"

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

class CompositionalMultiphaseFluid : public ConstitutiveBase
{
public:

  static constexpr localIndex MAX_NUM_COMPONENTS = 64;

  CompositionalMultiphaseFluid( std::string const & name, ManagedGroup * const parent );

  virtual ~CompositionalMultiphaseFluid() override;

  std::unique_ptr<ConstitutiveBase> DeliverClone( string const & name,
                                                  ManagedGroup * const parent ) const override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  static std::string CatalogName() { return dataRepository::keys::compositionalMultiphaseFluid; }

  virtual string GetCatalogName() override { return CatalogName(); }

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override final {}

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

  localIndex numFluidComponents() const;

  localIndex numFluidPhases() const;

  virtual void StateUpdatePointMultiphaseFluid(real64 const & pres,
                                               real64 const & temp,
                                               real64 const * composition,
                                               localIndex const k,
                                               localIndex const q) override;

  bool getMassFractionFlag() const      { return m_useMassFractions; }
  void setMassFractionFlag( bool flag ) { m_useMassFractions = flag; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phasesString           = "phases";
    static constexpr auto equationsOfStateString = "equationsOfState";
    static constexpr auto componentNamesString   = "componentNames";

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

    static constexpr auto phaseVolumeFractionString                      = "phaseVolumeFraction";                      // S_p
    static constexpr auto dPhaseVolumeFraction_dPressureString           = "dPhaseVolumeFraction_dPressure";           // dS_p/dP
    static constexpr auto dPhaseVolumeFraction_dTemperatureString        = "dPhaseVolumeFraction_dTemperature";        // dS_p/dT
    static constexpr auto dPhaseVolumeFraction_dGlobalCompFractionString = "dPhaseVolumeFraction_dGlobalCompFraction"; // dS_p/dz_c

    static constexpr auto phaseDensityString                      = "phaseDensity";                       // rho_p
    static constexpr auto dPhaseDensity_dPressureString           = "dPhaseDensity_dPressure";           // dRho_p/dP
    static constexpr auto dPhaseDensity_dTemperatureString        = "dPhaseDensity_dTemperature";        // dRho_p/dT
    static constexpr auto dPhaseDensity_dGlobalCompFractionString = "dPhaseDensity_dGlobalCompFraction"; // dRho_p/dz_c

    static constexpr auto phaseComponentFractionString                 = "phaseComponentFraction";                 // x_cp
    static constexpr auto dPhaseCompFraction_dPressureString           = "dPhaseCompFraction_dPressure";           // dx_cp/dP
    static constexpr auto dPhaseCompFraction_dTemperatureString        = "dPhaseCompFraction_dTemperature";        // dx_cp/dT
    static constexpr auto dPhaseCompFraction_dGlobalCompFractionString = "dPhaseCompFraction_dGlobalCompFraction"; // dx_cp/dz_c
    
    using ViewKey = dataRepository::ViewKey;

    ViewKey phases           = { phasesString };
    ViewKey equationsOfState = { equationsOfStateString };
    ViewKey componentNames   = { componentNamesString };

    ViewKey componentCriticalPressure    = { componentCriticalPressureString };
    ViewKey componentCriticalTemperature = { componentCriticalTemperatureString };
    ViewKey componentAcentricFactor      = { componentAcentricFactorString };
    ViewKey componentMolarWeight         = { componentMolarWeightString };
    ViewKey componentVolumeShift         = { componentVolumeShiftString };
    ViewKey componentBinaryCoeff         = { componentBinaryCoeffString };

    ViewKey phaseMoleFraction                              = { phaseMoleFractionString };                              // xi_p
    ViewKey dPhaseMoleFraction_dPressure                   = { dPhaseMoleFraction_dPressureString };                   // dXi_p/dP
    ViewKey dPhaseMoleFraction_dTemperature                = { dPhaseMoleFraction_dTemperatureString };                // dXi_p/dT
    ViewKey dPhaseMoleFraction_dGlobalCompMoleFraction     = { dPhaseMoleFraction_dGlobalCompFractionString };     // dXi_p/dz_c

    ViewKey phaseVolumeFraction                            = { phaseVolumeFractionString };                            // S_p
    ViewKey dPhaseVolumeFraction_dPressure                 = { dPhaseVolumeFraction_dPressureString };                 // dS_p/dP
    ViewKey dPhaseVolumeFraction_dTemperature              = { dPhaseVolumeFraction_dTemperatureString };              // dS_p/dT
    ViewKey dPhaseVolumeFraction_dGlobalCompMoleFraction   = { dPhaseVolumeFraction_dGlobalCompFractionString };   // dS_p/dz_c

    ViewKey phaseDensity                                   = { phaseDensityString };                                   // rho_p
    ViewKey dPhaseDensity_dPressure                        = { dPhaseDensity_dPressureString };                        // dRho_p/dP
    ViewKey dPhaseDensity_dTemperature                     = { dPhaseDensity_dTemperatureString };                     // dRho_p/dT
    ViewKey dPhaseDensity_dGlobalCompMoleFraction          = { dPhaseDensity_dGlobalCompFractionString };          // dRho_p/dz_c

    ViewKey phaseComponentMassFraction                     = { phaseComponentFractionString };                     // x_cp
    ViewKey dPhaseCompMassFraction_dPressure               = { dPhaseCompFraction_dPressureString };               // dx_cp/dP
    ViewKey dPhaseCompMassFraction_dTemperature            = { dPhaseCompFraction_dTemperatureString };            // dx_cp/dT
    ViewKey dPhaseCompMassFraction_dGlobalCompMoleFraction = { dPhaseCompFraction_dGlobalCompFractionString }; // dx_cp/dz_c

  } viewKeys;

private:

  template<typename T, int DIM>
  struct array_view_helper
  {
    using type = array_view<T, DIM>;
  };

#ifndef GEOSX_USE_ARRAY_BOUNDS_CHECK
  template<typename T>
  struct array_view_helper<T, 1>
  {
    using type = T *;
  };
#endif

  template<int DIM>
  using real_array_view = typename array_view_helper<real64, DIM>::type;

  // helper struct to represent a var and its derivatives
  template<int DIM>
  struct VarContainer
  {
    real_array_view<DIM>   value; // variable value
    real_array_view<DIM>   dPres; // derivative w.r.t. pressure
    real_array_view<DIM>   dTemp; // derivative w.r.t. temperature
    real_array_view<DIM+1> dComp; // derivative w.r.t. composition
  };

  void createFluid();

  // flag indicating whether input/output component fractions are treated as mass fractions
  bool m_useMassFractions;

  // general fluid composition information
  string_array m_phases;
  string_array m_equationsOfState;
  string_array m_componentNames;

  // standard EOS component input
  array1d<real64> m_componentCriticalPressure;
  array1d<real64> m_componentCriticalTemperature;
  array1d<real64> m_componentAcentricFactor;
  array1d<real64> m_componentMolarWeight;
  array1d<real64> m_componentVolumeShift;
  array2d<real64> m_componentBinaryCoeff;

  array3d<real64> m_phaseMoleFraction;
  array3d<real64> m_dPhaseMoleFraction_dPressure;
  array3d<real64> m_dPhaseMoleFraction_dTemperature;
  array4d<real64> m_dPhaseMoleFraction_dGlobalCompFraction;

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

  // PVTPackage fluid object
  PVTPackage::CompositionalMultiphaseSystem * m_fluid;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
