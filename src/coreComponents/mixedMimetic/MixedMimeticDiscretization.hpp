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
#include "common/format/EnumStrings.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

namespace geos
{

/**
 * @class MixedMimeticDiscretization
 *
 * Provides management of the cell-wise inner products and of the residual-based
 * Global Adaptation (GA) parameters used by the mixed mimetic finite difference solvers.
 */
class MixedMimeticDiscretization : public dataRepository::Group
{
public:

  /// Adaptation paradigm: no adaptation (single inner product everywhere) or Global Adaptation
  enum class AdaptationType : integer
  {
    None,   ///< the inner product selected by innerProductType is used in every cell
    Global  ///< face-assembled residual indicators select TPFA or MFD cell-wise
  };

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

    /// @return The key for the adaptation type
    static constexpr char const * adaptationTypeString() { return "adaptationType"; }

    /// @return The key for the residual tolerance
    static constexpr char const * residualToleranceString() { return "residualTolerance"; }

    /// @return The key for the nominal gradient of the projection probe
    static constexpr char const * nominalGradientString() { return "nominalGradient"; }
  };

  /**
   * @brief @return The adaptation type
   */
  AdaptationType getAdaptationType() const { return m_adaptationType; }

  /**
   * @brief @return The residual tolerance used in the marking criterion
   */
  real64 getResidualTolerance() const { return m_residualTolerance; }

  /**
   * @brief @return The nominal gradient used to build the projected admissible flow field
   */
  R1Tensor getNominalGradient() const { return m_nominalGradient; }

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  /// type of inner product used in the mixed mimetic solver
  string m_innerProductType;

  /// adaptation paradigm (none or global)
  AdaptationType m_adaptationType;

  /// user-prescribed tolerance for the marking criterion
  real64 m_residualTolerance;

  /// nominal gradient inducing the projected admissible flow field
  R1Tensor m_nominalGradient;

  /**
   * @brief Factory method to instantiate a type of mimetic inner product.
   * @return A unique_ptr< MimeticInnerProductBase > which contains the new instantiation.
   */
  std::unique_ptr< mimeticInnerProduct::MimeticInnerProductBase > factory( string const & mimeticInnerProductType ) const;

};

ENUM_STRINGS( MixedMimeticDiscretization::AdaptationType,
              "none",
              "global" );

}

#endif //GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATION_HPP_
