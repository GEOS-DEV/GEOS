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

/*
 * SourceFluxBoundaryCondition.hpp
 *
 */

#ifndef GEOS_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_
#define GEOS_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_

#include "FieldSpecification.hpp"

namespace geos
{

/**
 * @class SourceFluxBoundaryCondition
 * A class to manage Neumann boundary conditions
 */
class SourceFluxBoundaryCondition : public FieldSpecification
{
public:
  /**
   * @brief constructor
   * @param name the name of the FieldSpecification in the data repository
   * @param parent the parent group of this group.
   */
  SourceFluxBoundaryCondition( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  SourceFluxBoundaryCondition() = delete;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "SourceFlux"; }

  virtual const string getCatalogName() const override
  {
    return SourceFluxBoundaryCondition::catalogName();
  }

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecification::viewKeyStruct
  {
    /// @return The key for the temperature of the injected stream
    constexpr static char const * injectionTemperatureString() { return "injectionTemperature"; }
  };

  /**
   * @brief Getter for the temperature of the injected stream
   * @return the injection temperature [K], negative when it is not set
   */
  real64 getInjectionTemperature() const { return m_injectionTemperature; }

private:

  /// Temperature of the injected stream [K], used by the thermal flow solvers (negative: not set)
  real64 m_injectionTemperature;

};

} /* namespace geos */

#endif /* GEOS_FIELDSPECIFICATION_SOURCEFLUXBOUNDARYCONDITION_HPP_ */
