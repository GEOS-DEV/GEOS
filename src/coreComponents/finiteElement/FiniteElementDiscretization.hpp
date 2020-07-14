/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FiniteElementSpace.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEELEMENT_FINITEELEMENTSPACE_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "ElementLibrary/FiniteElement.h"
#include "FiniteElementDispatch.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

class NodeManager;
class CellBlockManager;
class ElementSubRegionBase;

namespace dataRepository
{
namespace keys
{
string const finiteElementSpace = "FiniteElementSpace";
string const basis = "basis";
string const quadrature = "quadrature";
string const dNdX = "dNdX";
string const detJ = "detJ";
string const parentSpace="parentSpace";
}
}

class BasisBase;
class QuadratureBase;
class FiniteElementBase;

class FiniteElementDiscretization : public dataRepository::Group
{
public:

  FiniteElementDiscretization() = delete;

  explicit FiniteElementDiscretization( std::string const & name, Group * const parent );

  ~FiniteElementDiscretization() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static string CatalogName() { return dataRepository::keys::finiteElementSpace; }

  ///@}


  std::unique_ptr< FiniteElementBase > getFiniteElement( string const & catalogName ) const;


  template< typename SUBREGION_TYPE >
  void CalculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                        SUBREGION_TYPE * const elementSubRegion ) const
  {
    GEOSX_MARK_FUNCTION;

    array4d< real64 > & dNdX = elementSubRegion->dNdX();
    array2d< real64 > & detJ = elementSubRegion->detJ();
    auto const & elemsToNodes = elementSubRegion->nodeList().toViewConst();

#if 1
    string const elementTypeString = elementSubRegion->GetElementTypeString();
    finiteElement::dispatch( elementTypeString,
                             [&] ( auto const finiteElement )
    {
      localIndex const numNodesPerElem = finiteElement.numNodes;
      localIndex const numQuadraturePointsPerElem = finiteElement.numQuadraturePoints;
      dNdX.resizeWithoutInitializationOrDestruction( elementSubRegion->size(), numQuadraturePointsPerElem, numNodesPerElem, 3 );
      detJ.resize( elementSubRegion->size(), numQuadraturePointsPerElem );

      for( localIndex k = 0; k < elementSubRegion->size(); ++k )
      {
        real64 xLocal[numNodesPerElem][3];
        for( localIndex a=0; a< numNodesPerElem; ++a )
        {
          localIndex const nodeIndex = elemsToNodes[ k][ a ];
          for( int i=0; i<3; ++i )
          {
            xLocal[ a ][ i ] = X[ nodeIndex ][ i ];
          }
        }



        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          real64 dNdXLocal[numNodesPerElem][3];
          detJ( k, q ) = finiteElement.shapeFunctionDerivatives( q, xLocal, dNdXLocal );

          for( localIndex b = 0; b < numNodesPerElem; ++b )
          {
            LvArray::tensorOps::copy< 3 >( dNdX[ k ][ q ][ b ], dNdXLocal[b] );
          }
        }
      }

    } );
#else

    {
      std::unique_ptr< FiniteElementBase > fe = getFiniteElement( m_parentSpace );
      dNdX.resizeWithoutInitializationOrDestruction( elementSubRegion->size(), m_quadrature->size(), fe->dofs_per_element(), 3 );
      detJ.resize( elementSubRegion->size(), m_quadrature->size() );
    }

    // We can't use a simple RAJA loop here because each thread needs it's own FiniteElementBase, and we can't capture
    // fe
    // in the lambda because it's a unique_ptr.
    PRAGMA_OMP( "omp parallel" )
    {
      std::unique_ptr< FiniteElementBase > fe = getFiniteElement( m_parentSpace );

      PRAGMA_OMP( "omp for" )
      for( localIndex k = 0; k < elementSubRegion->size(); ++k )
      {
        fe->reinit( X, elemsToNodes[k] );

        for( localIndex q = 0; q < fe->n_quadrature_points(); ++q )
        {
          detJ( k, q ) = fe->JxW( q );
          for( localIndex b = 0; b < fe->dofs_per_element(); ++b )
          {
            LvArray::tensorOps::copy< 3 >( dNdX[ k ][ q ][ b ], fe->gradient( b, q ) );
          }
        }
      }
    }
#endif
  }

  localIndex getNumberOfQuadraturePoints() const;

  string m_basisName;
  string m_quadratureName;
  string m_parentSpace;

  BasisBase const *    m_basis    = nullptr;
  QuadratureBase const * m_quadrature = nullptr;
protected:
  void PostProcessInput() override final;

};

} /* namespace geosx */

#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_ */
