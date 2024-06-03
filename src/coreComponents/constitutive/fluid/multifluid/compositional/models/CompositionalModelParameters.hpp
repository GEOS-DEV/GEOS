/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalModelParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALMODELPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALMODELPARAMETERS_HPP_

#include "ModelParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ConstantViscosity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
struct CompositionalModelParameters
{
  static std::unique_ptr< ModelParameters > createModelParameters()
  {
    return std::make_unique< ModelParameters >();
  }
};

// Constant viscosity parameters
template<>
struct CompositionalModelParameters<
  CompositionalTwoPhasePengRobinsonConstantViscosity::FlashModel,
  CompositionalTwoPhasePengRobinsonConstantViscosity::Phase1Model,
  CompositionalTwoPhasePengRobinsonConstantViscosity::Phase2Model,
  CompositionalTwoPhasePengRobinsonConstantViscosity::Phase3Model >
{
  static std::unique_ptr< ModelParameters > createModelParameters()
  {
    return std::make_unique< ConstantViscosity::Parameters >();
  }
};

template<>
struct CompositionalModelParameters<
  CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity::FlashModel,
  CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity::Phase1Model,
  CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity::Phase2Model,
  CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity::Phase3Model >
{
  static std::unique_ptr< ModelParameters > createModelParameters()
  {
    return std::make_unique< ConstantViscosity::Parameters >();
  }
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALMODELPARAMETERS_HPP_
