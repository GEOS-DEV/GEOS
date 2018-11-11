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
  * @file RelativePermeabilityBase.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"

namespace geosx
{

namespace constitutive
{

class RelativePermeabilityBase : public ConstitutiveBase
{
public:

  struct PhaseType
  {
    static constexpr integer OIL            = 0;
    static constexpr integer GAS            = 1;
    static constexpr integer WATER          = 2;
    static constexpr integer MAX_NUM_PHASES = 3;
  };

  RelativePermeabilityBase( std::string const & name, ManagedGroup * const parent );

  virtual ~RelativePermeabilityBase() override;

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override;

  virtual void AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                         localIndex const numPts ) override;

  virtual void StateUpdate( dataRepository::ManagedGroup const * const input,
                            dataRepository::ManagedGroup const * const parameters,
                            dataRepository::ManagedGroup * const stateVariables,
                            integer const systemAssembleFlag ) const override {}

  // RelPerm-specific interface

  virtual void StateUpdatePointRelPerm( arraySlice1d<real64> const & phaseVolFraction ) {}

  localIndex numFluidPhases() const;

  string const & phaseName( localIndex ip ) const;

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto phaseNamesString     = "phaseNames";
    static constexpr auto phaseTypesString     = "phaseTypes";
    static constexpr auto phaseOrderString     = "phaseTypes";

    static constexpr auto phaseRelPermString                    = "phaseRelPermString";              // Kr
    static constexpr auto dPhaseRelPerm_dPhaseVolFractionString = "dPhaseRelPerm_dPhaseVolFraction"; // dKr_p/dS_p

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseNames = { phaseNamesString };
    ViewKey phaseTypes = { phaseTypesString };
    ViewKey phaseOrder = { phaseOrderString };

    ViewKey phaseRelPerm                    = { phaseRelPermString };                    // Kr_p
    ViewKey dPhaseRelPerm_dPhaseVolFraction = { dPhaseRelPerm_dPhaseVolFractionString }; // dKr_p/dS_p

  } viewKeysRelativePermeabilityBase;

protected:

  // phase names read from input
  string_array     m_phaseNames;

  // phase ordering info
  array1d<integer> m_phaseTypes;
  array1d<integer> m_phaseOrder;

  // output quantities
  array3d<real64>  m_phaseRelPerm;
  array4d<real64>  m_dPhaseRelPerm_dPhaseVolFrac;

};

} // namespace constitutive

} // namespace geosx


#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_RELATIVEPERMEABILITYBASE_HPP
