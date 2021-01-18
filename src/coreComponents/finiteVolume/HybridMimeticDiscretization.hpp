/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HybridMimeticDiscretization.hpp
 */

#ifndef GEOSX_FINITEVOLUME_HYBRIDMIMETICDISCRETIZATION_HPP_
#define GEOSX_FINITEVOLUME_HYBRIDMIMETICDISCRETIZATION_HPP_

#include "dataRepository/Group.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

namespace geosx
{

/**
 * @class HybridMimeticDiscretization
 *
 * Provides management of the inner product when using a hybrid FVM solver
 */
class HybridMimeticDiscretization : public dataRepository::Group
{
public:

  /// Alias for CatalogInterface, necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface< HybridMimeticDiscretization, string const &, Group * const >;

  /**
   * @brief Return the data type in the data repository.
   * @return the data type in the data repository
   */
  static typename CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static std::string CatalogName() { return "HybridMimeticDiscretization"; }

  HybridMimeticDiscretization() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the HybridMimeticDiscretization in the data repository
   * @param parent the parent group of this group.
   */
  HybridMimeticDiscretization( std::string const & name, dataRepository::Group * const parent );

  /**
   * @brief View keys.
   */
  struct viewKeyStruct
  {
    /// The key for coefficientName
    static constexpr auto coeffNameString        = "coefficientName";
    /// The key for transMultiplier
    static constexpr auto transMultiplierString  = "TransMultiplier";
    /// The key for the type of inner product
    static constexpr auto innerProductTypeString = "innerProductType";
    /// The key for the inner product
    static constexpr auto innerProductString     = "innerProduct";
  };

protected:

  virtual void RegisterDataOnMesh( Group * const meshBodies ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const rootGroup ) override;

private:

  /// name of the coefficient field
  string m_coeffName;

  /// type of of inner product used in the hybrid FVM solver
  string m_innerProductType;

  /**
   * @brief Factory method to instantiate a type of mimetic inner product.
   * @return A unique_ptr< MimeticInnerProductBase > which contains the new
   *   instantiation.
   */
  std::unique_ptr< mimeticInnerProduct::MimeticInnerProductBase > factory( string const & mimeticInnerProductType ) const;

};

}


#endif //GEOSX_FINITEVOLUME_HYBRIDMIMETICDISCRETIZATION_HPP_
