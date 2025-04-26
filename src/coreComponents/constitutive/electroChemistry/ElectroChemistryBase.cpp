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
 * @file ElectroChemistryBase.cpp
 */

#include "constitutive/electroChemistry/ElectroChemistryBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ElectroChemistryBase::ElectroChemistryBase(string const& name, Group * const parent)
:
ConstitutiveBase(name, parent),
m_conductivity()
{
  registerWrapper(viewKeyStruct::defaultConductivityString(), &m_defaultConductivity).
    setApplyDefaultValue(1.0).
    setInputFlag(InputFlags::OPTIONAL).
    setDescription("Default Electro Conductivity");

  registerWrapper(viewKeyStruct::conductivityString(), &m_conductivity).
    setApplyDefaultValue(-1.0). // will be overwritten
    setDescription("Electro Conductivity Field");
}

ElectroChemistryBase::~ElectroChemistryBase() {}

void ElectroChemistryBase::postInputInitialization()
{
  this->getWrapper<array1d<real64>>(viewKeyStruct::conductivityString()).
    setApplyDefaultValue(m_defaultConductivity);
}

REGISTER_CATALOG_ENTRY(ConstitutiveBase, ElectroChemistryBase, string const&, Group * const)
} // namespace constitutive

} // namespace geos