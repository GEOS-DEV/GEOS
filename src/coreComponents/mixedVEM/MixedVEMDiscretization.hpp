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
 * @file MixedVEMDiscretization.hpp
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMDISCRETIZATION_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMDISCRETIZATION_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @class MixedVEMDiscretization
 *
 * Input object of the lowest-order mixed virtual element discretization of elasticity.
 *
 * The element operators K_E and B_E are the same either way. What the hybridization flag
 * selects is how the global problem is formed and solved: the mixed form assembles the
 * indefinite saddle point system in the face tractions and the element displacements,
 * whereas the hybridized form breaks the stress space, condenses both element unknowns
 * out and solves the symmetric positive definite interface problem for the multiplier,
 * recovering (sigma_E, u_E) elementwise afterwards.
 */
class MixedVEMDiscretization : public dataRepository::Group
{
public:

  /// Alias for CatalogInterface, necessary declarations for factory instantiation of derived classes
  using CatalogInterface = dataRepository::CatalogInterface< MixedVEMDiscretization, string const &, Group * const >;

  /**
   * @brief Return the data type in the data repository.
   * @return the data type in the data repository
   */
  static typename CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static string catalogName() { return "MixedVEMDiscretization"; }

  MixedVEMDiscretization() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the MixedVEMDiscretization in the data repository
   * @param parent the parent group of this group.
   */
  MixedVEMDiscretization( string const & name, dataRepository::Group * const parent );

  /**
   * @brief View keys.
   */
  struct viewKeyStruct
  {
    /// @return The key for hybridization
    static constexpr char const * hybridizationString() { return "hybridization"; }
  };

  /**
   * @brief @return whether the hybridized form is requested.
   */
  bool useHybridization() const { return m_hybridization > 0; }

private:

  /// flag selecting the hybridized form instead of the mixed saddle point form
  integer m_hybridization;

};

}

#endif // GEOS_MIXEDVEM_MIXEDVEMDISCRETIZATION_HPP_
