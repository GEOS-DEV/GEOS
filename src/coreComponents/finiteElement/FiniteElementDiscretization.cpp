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
 * @file FiniteElementSpace.cpp
 */

#include "FiniteElementDiscretization.hpp"

// TODO make this not dependent on this header...need better key implementation

namespace geosx
{
using namespace dataRepository;
using namespace finiteElement;


FiniteElementDiscretization::FiniteElementDiscretization( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::orderString(), &m_order ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The order of the finite element basis." );

  registerWrapper( viewKeyStruct::formulationString(), &m_formulation ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "default" ).
    setDescription( "Specifier to indicate any specialized formuations. "
                    "For instance, one of the many enhanced assumed strain "
                    "methods of the Hexahedron parent shape would be indicated "
                    "here" );

  registerWrapper( viewKeyStruct::useVemString(), &m_useVem ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Specifier to indicate whether to force the use of VEM" );
}

FiniteElementDiscretization::~FiniteElementDiscretization()
{}


void FiniteElementDiscretization::postProcessInput()
{
//  GEOSX_ERROR_IF_NE_MSG( m_order, 1, "Higher order finite element spaces are currently not supported." );
  GEOSX_ERROR_IF_NE_MSG( m_formulation, "default", "Only standard element formulations are currently supported." );
  GEOSX_ERROR_IF_GT_MSG( m_useVem, 1, "The flag useVirtualElements can be either 0 or 1" );
}

std::unique_ptr< FiniteElementBase >
FiniteElementDiscretization::factory( ElementType const parentElementShape ) const
{
  if( m_order==1 )
  {
    switch( parentElementShape )
    {
      case ElementType::Triangle:      return std::make_unique< H1_TriangleFace_Lagrange1_Gauss1 >();
      case ElementType::Quadrilateral: return std::make_unique< H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
      case ElementType::Tetrahedron:
      {
        if( m_useVem == 1 )
        {
          return std::make_unique< H1_Tetrahedron_VEM_Gauss1 >();
        }
        else
        {
          return std::make_unique< H1_Tetrahedron_Lagrange1_Gauss1 >();
        }
      }
      case ElementType::Pyramid:
      {
        if( m_useVem == 1 )
        {
          return std::make_unique< H1_Pyramid_VEM_Gauss1 >();
        }
        else
        {
          return std::make_unique< H1_Pyramid_Lagrange1_Gauss5 >();
        }
      }
      case ElementType::Wedge:
      {
        if( m_useVem == 1 )
        {
          return std::make_unique< H1_Wedge_VEM_Gauss1 >();
        }
        else
        {
          return std::make_unique< H1_Wedge_Lagrange1_Gauss6 >();
        }
      }
      case ElementType::Hexahedron:
      {
        if( m_useVem == 1 )
        {
          return std::make_unique< H1_Hexahedron_VEM_Gauss1 >();
        }
        else
        {
          return std::make_unique< H1_Hexahedron_Lagrange1_GaussLegendre2 >();
        }
      }
      case ElementType::Prism5:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism5 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism5_VEM_Gauss1 >();
      }
      case ElementType::Prism6:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism6 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism6_VEM_Gauss1 >();
      }
      case ElementType::Prism7:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism7 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism7_VEM_Gauss1 >();
      }
      case ElementType::Prism8:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism8 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism8_VEM_Gauss1 >();
      }
      case ElementType::Prism9:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism9 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism9_VEM_Gauss1 >();
      }
      case ElementType::Prism10:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism10 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism10_VEM_Gauss1 >();
      }
      case ElementType::Prism11:
      {
        GEOSX_ERROR_IF( m_useVem != 1,
                        "Element type Prism11 available only when using the Virtual Element Method" );
        return std::make_unique< H1_Prism11_VEM_Gauss1 >();
      }
      default:
      {
        GEOSX_ERROR( "Element type " << parentElementShape << " does not have an associated element formulation." );
      }
    }
    return {};
  }

  if( m_order==3 )
  {
    switch( parentElementShape )
    {
      case ElementType::Hexahedron:
        return std::make_unique< Q3_Hexahedron_Lagrange_GaussLobatto >();
      default:
      {
        GEOSX_ERROR( "Element type " << parentElementShape << " does not have an associated element formulation." );
      }
    }
    return {};
  }

  GEOSX_ERROR( "Element type " << parentElementShape << " does not have an associated element formulation." );
  return {};
}

REGISTER_CATALOG_ENTRY( Group, FiniteElementDiscretization, string const &, Group * const )

} /* namespace geosx */
