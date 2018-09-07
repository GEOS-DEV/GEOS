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

  localIndex numFluidComponents();

  localIndex numFluidPhases();

  virtual void StateUpdatePointMultiphaseFluid(real64 const & pres,
                                               real64 const & temp,
                                               real64 const * composition,
                                               localIndex const k,
                                               localIndex const q) override;


  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    using ViewKey = dataRepository::ViewKey;

    // inputs
    ViewKey phases           = { "phases" };
    ViewKey equationsOfState = { "equationsOfState" };
    ViewKey componentNames   = { "componentNames" };

    ViewKey componentCriticalPressure    = { "componentCriticalPressure" };
    ViewKey componentCriticalTemperature = { "componentCriticalTemperature" };
    ViewKey componentAcentricFactor      = { "componentAcentricFactor" };
    ViewKey componentMolarWeight         = { "componentMolarWeight" };
    ViewKey componentVolumeShift         = { "componentVolumeShift" };
    ViewKey componentBinaryCoeff         = { "componentBinaryCoeff" };

    // constitutive data
    ViewKey phaseMoleFraction            = { "phaseMoleFraction" };          // xi_p
    ViewKey phaseVolumeFraction          = { "phaseVolumeFraction" };        // S_p
    ViewKey phaseDensity                 = { "phaseDensity" };               // rho_p
    ViewKey phaseComponentMoleFraction   = { "phaseComponentMoleFraction" }; // x_cp
    ViewKey phaseComponentDensity        = { "phaseComponentDensity" };      // rho_cp

    // derivatives
    ViewKey dPhaseMoleFraction_dPressure                   = { "dPhaseMoleFraction_dPressure" };                   // dXi_p/dP
    ViewKey dPhaseMoleFraction_dTemperature                = { "dPhaseMoleFraction_dTemperature" };                // dXi_p/dT
    ViewKey dPhaseMoleFraction_dGlobalCompMoleFraction     = { "dPhaseMoleFraction_dGlobalCompMoleFraction" };     // dXi_p/dz_c
    ViewKey dPhaseVolumeFraction_dPressure                 = { "dPhaseVolumeFraction_dPressure" };                 // dS_p/dP
    ViewKey dPhaseVolumeFraction_dTemperature              = { "dPhaseVolumeFraction_dTemperature" };              // dS_p/dT
    ViewKey dPhaseVolumeFraction_dGlobalCompMoleFraction   = { "dPhaseVolumeFraction_dGlobalCompMoleFraction" };   // dS_p/dz_c
    ViewKey dPhaseDensity_dPressure                        = { "dPhaseDensity_dPressure" };                        // dRho_p/dP
    ViewKey dPhaseDensity_dTemperature                     = { "dPhaseDensity_dTemperature" };                     // dRho_p/dT
    ViewKey dPhaseDensity_dGlobalCompMoleFraction          = { "dPhaseDensity_dGlobalCompMoleFraction" };          // dRho_p/dz_c
    ViewKey dPhaseCompMoleFraction_dPressure               = { "dPhaseCompMoleFraction_dPressure" };               // dx_cp/dP
    ViewKey dPhaseCompMoleFraction_dTemperature            = { "dPhaseCompMoleFraction_dTemperature" };            // dx_cp/dT
    ViewKey dPhaseCompMoleFraction_dGlobalCompMoleFraction = { "dPhaseCompMoleFraction_dGlobalCompMoleFraction" }; // dx_cp/dz_c
    ViewKey dPhaseCompDensity_dPressure                    = { "dPhaseCompDensity_dPressure" };                    // dRho_cp/dP
    ViewKey dPhaseCompDensity_dTemperature                 = { "dPhaseCompDensity_dTemperature" };                 // dRho_cp/dT
    ViewKey dPhaseCompDensity_dGlobalCompMoleFraction      = { "dPhaseCompDensity_dGlobalCompMoleFraction" };      // dRho_cp/dz_c

  } viewKeys;

private:

  // helper struct to represent a var and its derivatives
  template<int DIM>
  struct VarContainer
  {
    array_view<real64, DIM>   value; // variable value
    array_view<real64, DIM>   dPres; // derivative w.r.t. pressure
    array_view<real64, DIM>   dTemp; // derivative w.r.t. temperature
    array_view<real64, DIM+1> dComp; // derivative w.r.t. composition
  };

  void createFluid();

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
  array4d<real64> m_dPhaseMoleFraction_dGlobalCompMoleFraction;

  array3d<real64> m_phaseVolumeFraction;
  array3d<real64> m_dPhaseVolumeFraction_dPressure;
  array3d<real64> m_dPhaseVolumeFraction_dTemperature;
  array4d<real64> m_dPhaseVolumeFraction_dGlobalCompMoleFraction;

  array3d<real64> m_phaseDensity;
  array3d<real64> m_dPhaseDensity_dPressure;
  array3d<real64> m_dPhaseDensity_dTemperature;
  array4d<real64> m_dPhaseDensity_dGlobalCompMoleFraction;

  array4d<real64> m_phaseCompMoleFraction;
  array4d<real64> m_dPhaseCompMoleFraction_dPressure;
  array4d<real64> m_dPhaseCompMoleFraction_dTemperature;
  array5d<real64> m_dPhaseCompMoleFraction_dGlobalCompMoleFraction;

  array4d<real64> m_phaseCompDensity;
  array4d<real64> m_dPhaseCompDensity_dPressure;
  array4d<real64> m_dPhaseCompDensity_dTemperature;
  array5d<real64> m_dPhaseCompDensity_dGlobalCompMoleFraction;

  // PVTPackage fluid object
  PVTPackage::CompositionalMultiphaseSystem * m_fluid;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
