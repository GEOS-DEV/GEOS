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
 * @file ButlerVolmerReaction.cpp
 */

#include "constitutive/electroChemistry/ButlerVolmerReaction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ButlerVolmerInterface::ButlerVolmerInterface(string const& name, Group * const parent)
:
ConstitutiveBase(name, parent),
m_krxn()
{
  registerWrapper(viewKeyStruct::defaultReactivityCoefficientString(), &m_defaultKrxn).
    setApplyDefaultValue(1.0).
    setInputFlag(InputFlags::REQUIRED).
    setDescription("Default BV Reactivity Coefficient");

  registerWrapper(viewKeyStruct::reactivityCoefficientString(), &m_krxn).
    setApplyDefaultValue(-1.0). // will be overwritten
    setPlotLevel(PlotLevel::LEVEL_0).
    setDescription("Reactivity Coefficient Field");
}

ButlerVolmerInterface::~ButlerVolmerInterface() {}

void ButlerVolmerInterface::postInputInitialization()
{
  this->getWrapper<array1d<real64>>(viewKeyStruct::reactivityCoefficientString()).
    setApplyDefaultValue(m_defaultKrxn);
}

REGISTER_CATALOG_ENTRY(ConstitutiveBase, ButlerVolmerInterface, string const&, Group * const)
} // namespace constitutive

} // namespace geos