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
 * @file StencilAccessors.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STENCILACCESSORS_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STENCILACCESSORS_HPP_

#include "mesh/ElementRegionManager.hpp"
#include "codingUtilities/traits.hpp"
#include "codingUtilities/Utilities.hpp"

#include <tuple>

namespace geosx
{

/**
 * @brief A struct to automatically construct and store element view accessors
 * @struct StencilAccessors
 * @tparam TRAITS the pack containing the types of the fields
 */
template< typename ... TRAITS >
class StencilAccessors
{
public:

  /**
   * @brief @return reference to the accessor for the given field
   * @tparam TRAIT the field trait type
   */
  template< typename TRAIT >
  auto get( TRAIT ) const
  {
    constexpr std::size_t idx = traits::type_list_index< TRAIT, std::tuple< TRAITS ... > >;
    static_assert( idx != std::tuple_size< std::tuple< TRAITS... > >::value, "input trait/stencil does not match the available traits/stencils." );
    return std::get< idx >( m_accessors ).toNestedViewConst();
  }

  template< typename TRAIT >
  auto get() const
  {
    return get( TRAIT{} );
  }

  /**
   * @brief Constructor for the struct
   * @param[in] elemManager a reference to the elemRegionManager
   * @param[in] solverName the name of the solver creating the view accessors
   */
  StencilAccessors( ElementRegionManager const & elemManager,
                    string const & solverName )
  {
    forEachArgInTuple( std::tuple< TRAITS ... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      using TRAIT = TYPEOFREF( t );

      auto & acc = std::get< idx() >( m_accessors );
      acc = elemManager.constructExtrinsicAccessor< TRAIT >();
      acc.setName( solverName + "/accessors/" + TRAIT::key() );
    } );
  }

protected:

  /// the tuple storing all the accessors
  std::tuple< ElementRegionManager::ElementViewAccessor< traits::ViewTypeConst< typename TRAITS::type > > ... > m_accessors;

  /**
   * @brief Constructor for the struct
   */
  StencilAccessors() = default;
};

/**
 * @brief A struct to automatically construct and store element view accessors
 * @struct StencilMaterialAccessors
 * @tparam MATERIAL_TYPE the type of the material model
 * @tparam TRAITS the pack containing the types of the fields
 */
template< typename MATERIAL_TYPE, typename ... TRAITS >
class StencilMaterialAccessors : public StencilAccessors< TRAITS ... >
{
public:

  /**
   * @brief Constructor for the struct
   * @tparam MATERIAL_TYPE  type of the constitutive model
   * @param[in] elemManager a reference to the elemRegionManager
   * @param[in] solverName the name of the solver creating the view accessors
   */
  StencilMaterialAccessors( ElementRegionManager const & elemManager,
                            string const & solverName ):
    StencilAccessors< TRAITS ... >()
  {
    forEachArgInTuple( std::tuple< TRAITS ... >{}, [&]( auto t, auto idx )
    {
      GEOSX_UNUSED_VAR( t );
      using TRAIT = TYPEOFREF( t );

      auto & acc = std::get< idx() >( m_accessors );
      bool const allowMissingViews = false;
      acc = elemManager.constructMaterialExtrinsicAccessor< MATERIAL_TYPE, TRAIT >( allowMissingViews );
      acc.setName( solverName + "/accessors/" + TRAIT::key() );
    } );
  }
};


}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_STENCILACCESSORS_HPP_
