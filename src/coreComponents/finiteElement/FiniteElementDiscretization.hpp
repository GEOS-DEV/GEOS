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
 * @file FiniteElementSpace.hpp
 */

#ifndef GEOS_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_
#define GEOS_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "mesh/ElementType.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "FiniteElementDispatch.hpp"

namespace geos
{

class FiniteElementDiscretization : public dataRepository::Group
{
public:

  /// Enumerator of available interpolation types
  enum class Formulation : integer
  {
    Default,
    SEM,
    DG
  };


  FiniteElementDiscretization() = delete;

  explicit FiniteElementDiscretization( string const & name, Group * const parent );

  ~FiniteElementDiscretization() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static string catalogName() { return "FiniteElementSpace"; }

  ///@}

  /**
   * @brief Factory method to instantiate a type of finite element formulation.
   * @param parentElementShape String key that indicates the type of
   *   element/basis/formulation that should be instantiated.
   * @return A unique_ptr< FinteElementBase > which contains the new
   *   instantiation.
   */
  std::unique_ptr< finiteElement::FiniteElementBase >
  factory( ElementType const parentElementShape ) const;

  int getOrder() const { return m_order; }

  Formulation getFormulation() const { return m_formulation; }

private:

  struct viewKeyStruct
  {
    static constexpr char const * orderString() { return "order"; }
    static constexpr char const * formulationString() { return "formulation"; }
    static constexpr char const * useVemString() { return "useVirtualElements"; }
    static constexpr char const * useHighOrderQuadratureRuleString() { return "useHighOrderQuadratureRule"; }
  };

  /// The order of the finite element basis
  int m_order;

  /// Optional string indicating any specialized formulation type.
  Formulation m_formulation;

  /// Optional parameter indicating if the class should use Virtual Elements.
  int m_useVem;

  /// Optional parameter indicating if the class should use a high order quadrature rule.
  int m_useHighOrderQuadratureRule;

  void postInputInitialization() override final;

};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( FiniteElementDiscretization::Formulation,
              "default",
              "SEM",
              "DG" );

} /* namespace geos */

#endif /* GEOS_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_ */
