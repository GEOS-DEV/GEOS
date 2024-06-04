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
 * @file ModelParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_MODELPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_MODELPARAMETERS_HPP_

namespace geos
{

namespace constitutive
{

class MultiFluidBase;

namespace compositional
{

class ComponentProperties;

class ModelParameters
{
public:
  ModelParameters() = default;
  virtual ~ModelParameters() = default;

  virtual void registerParameters( MultiFluidBase * fluid )
  {
    GEOS_UNUSED_VAR( fluid );
  }

  virtual void postProcessInput( MultiFluidBase const * fluid, ComponentProperties const & componentProperties )
  {
    GEOS_UNUSED_VAR( fluid );
    GEOS_UNUSED_VAR( componentProperties );
  }
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_MODELPARAMETERS_HPP_
