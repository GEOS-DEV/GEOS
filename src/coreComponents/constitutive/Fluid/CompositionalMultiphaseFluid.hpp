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
    ViewKey componentMolarWeigth         = { "componentMolarWeigth" };
    ViewKey componentVolumeShift         = { "componentVolumeShift" };
    ViewKey componentBinaryCoeff         = { "componentBinaryCoeff" };

    // constitutive data
    ViewKey phaseVolumeFraction          = { "phaseVolumeFraction" };         // S_p
    ViewKey phaseDensity                 = { "phaseDensity" };                // rho_p
    ViewKey phaseComponentMoleFraction   = { "phaseComponentMoleFraction" };  // x_cp
    ViewKey phaseComponentDensity        = { "phaseComponentDensity" };       // rho_cp

    // derivatives
    ViewKey dPhaseVolumeFraction_dPhasePressure            = { "dPhaseVolumeFraction_dPhasePressure" };            // dS_p/dP_p
    ViewKey dPhaseVolumeFraction_dGlobalCompMoleFraction   = { "dPhaseVolumeFraction_dGlobalCompMoleFraction" };   // dS_p/dz_c
    ViewKey dPhaseDensity_dPhasePressure                   = { "dPhaseDensity_dPhasePressure" };                   // dRho_p/dP_p
    ViewKey dPhaseDensity_dGlobalCompMoleFraction          = { "dPhaseDensity_dGlobalCompMoleFraction" };          // dRho_p/dz_c
    ViewKey dPhaseCompMoleFraction_dPhasePressure          = { "dPhaseCompMoleFraction_dPhasePressure" };          // dx_cp/dP_p
    ViewKey dPhaseCompMoleFraction_dGlobalCompMoleFraction = { "dPhaseCompMoleFraction_dGlobalCompMoleFraction" }; // dx_cp/dz_c
    ViewKey dPhaseCompDensity_dPhasePressure               = { "dPhaseCompDensity_dPhasePressure" };               // dRho_cp/dP_p
    ViewKey dPhaseCompDensity_dGlobalCompMoleFraction      = { "dPhaseCompDensity_dGlobalCompMoleFraction" };      // dRho_cp/dz_c

  } viewKeys;

private:

  void createFluid();

  // general fluid composition information
  string_array m_phases;
  string_array m_equationsOfState;
  string_array m_componentNames;

  // standard EOS component input
  array1d<real64> m_componentCriticalPressure;
  array1d<real64> m_componentCriticalTemperature;
  array1d<real64> m_componentAcentricFactor;
  array1d<real64> m_componentMolarWeigth;
  array1d<real64> m_componentVolumeShift;
  array2d<real64> m_componentBinaryCoeff;

  // data storage
  array3d<real64> m_phaseVolumeFraction;
  array3d<real64> m_phaseDensity;
  array4d<real64> m_phaseCompMoleFraction;
  array4d<real64> m_phaseCompDensity;

  // derivatives
  array3d<real64> m_dPhaseVolumeFraction_dPhasePressure;
  array4d<real64> m_dPhaseVolumeFraction_dGlobalCompMoleFraction;

  array3d<real64> m_dPhaseDensity_dPhasePressure;
  array4d<real64> m_dPhaseDensity_dGlobalCompMoleFraction;

  array4d<real64> m_dPhaseCompMoleFraction_dPhasePressure;
  array5d<real64> m_dPhaseCompMoleFraction_dGlobalCompMoleFraction;

  array4d<real64> m_dPhaseCompDensity_dPhasePressure;
  array5d<real64> m_dPhaseCompDensity_dGlobalCompMoleFraction;

  // PVTPackage fluid object
  PVTPackage::CompositionalMultiphaseSystem * m_fluid;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_COMPOSITIONALMULTIPHASEFLUID_HPP_
