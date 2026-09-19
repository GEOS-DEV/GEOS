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
 * @file MixedMimeticDiscretization.hpp
 */

#ifndef GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATION_HPP_
#define GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATION_HPP_

#include "dataRepository/Group.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

namespace geos
{

/**
 * @class MixedMimeticDiscretization
 *
 * Cell-wise inner product of the mixed mimetic finite difference solvers and parameters of
 * its selection eta (ConsistencyAdaptation): consistency tolerance of the residual layer and
 * degeneracy tolerance; without adaptation the selected inner product is used in every cell.
 */
class MixedMimeticDiscretization : public dataRepository::Group
{
public:

  /// Alias for CatalogInterface, necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface< MixedMimeticDiscretization, string const &, Group * const >;

  /**
   * @brief Return the data type in the data repository.
   * @return the data type in the data repository
   */
  static typename CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static string catalogName() { return "MixedMimeticDiscretization"; }

  MixedMimeticDiscretization() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the MixedMimeticDiscretization in the data repository
   * @param parent the parent group of this group.
   */
  MixedMimeticDiscretization( string const & name, dataRepository::Group * const parent );

  /**
   * @brief View keys.
   */
  struct viewKeyStruct
  {
    /// @return The key for the type of inner product
    static constexpr char const * innerProductTypeString() { return "innerProductType"; }

    /// @return The key for the inner product
    static constexpr char const * innerProductString() { return "innerProduct"; }

    /// @return The key for the adaptive consistency flag
    static constexpr char const * adaptiveConsistencyString() { return "adaptiveConsistency"; }

    /// @return The key for the consistency tolerance
    static constexpr char const * consistencyToleranceString() { return "consistencyTolerance"; }

    /// @return The key for the nominal gradient of the projection probe
    static constexpr char const * nominalGradientString() { return "nominalGradient"; }

    /// @return The key for the degeneracy tolerance
    static constexpr char const * degeneracyToleranceString() { return "degeneracyTolerance"; }
  };

  /**
   * @brief @return Whether the consistency layer is enabled
   */
  bool isAdaptiveConsistency() const { return m_adaptiveConsistency == 1; }

  /**
   * @brief @return Whether the selected inner product is the (diagonal) TPFA inner product
   */
  bool isTpfaInnerProduct() const;

  /**
   * @brief @return The tolerance of the consistency layer: eta = 1 where the indicator exceeds it
   */
  real64 getConsistencyTolerance() const { return m_consistencyTolerance; }

  /**
   * @brief @return The nominal gradient used to build the projected admissible flow field
   */
  R1Tensor getNominalGradient() const { return m_nominalGradient; }

  /**
   * @brief @return The degeneracy tolerance in percent of the node-star volume
   */
  real64 getDegeneracyTolerance() const { return m_degeneracyTolerance; }

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  /// type of inner product used in the mixed mimetic solver
  string m_innerProductType;

  /// flag enabling the consistency layer (1 = adaptive, 0 = selected inner product everywhere)
  integer m_adaptiveConsistency;

  /// tolerance of the consistency layer
  real64 m_consistencyTolerance;

  /// nominal gradient inducing the projected admissible flow field
  R1Tensor m_nominalGradient;

  /// cells whose volume is below this percentage of the volume of their node star use the diagonal product
  real64 m_degeneracyTolerance;

  /**
   * @brief Factory method to instantiate a type of mimetic inner product.
   * @return A unique_ptr< MimeticInnerProductBase > which contains the new instantiation.
   */
  std::unique_ptr< mimeticInnerProduct::MimeticInnerProductBase > factory( string const & mimeticInnerProductType ) const;

};

}

#endif //GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATION_HPP_
