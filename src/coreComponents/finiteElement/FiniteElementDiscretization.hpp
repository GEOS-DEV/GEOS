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

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "LvArray/src/tensorOps.hpp"
//#include "ElementLibrary/FiniteElement.h"
#include "FiniteElementDispatch.hpp"


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

  template< typename SUBREGION_TYPE,
            typename FE_TYPE >
  void CalculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                        SUBREGION_TYPE * const elementSubRegion,
                                        FE_TYPE & fe ) const;


  std::unique_ptr<FiniteElementShapeFunctionKernelBase>
  factory( string const & parentElementShape ) const;

  struct viewKeyStruct
  {
    /// String key to adjacency list
    static constexpr auto orderString = "order";
    static constexpr auto formulationString = "formulation";
  };

  string m_basisName;
  string m_quadratureName;
  string m_parentSpace;

  int m_order;
  string m_formulation;

protected:
  void PostProcessInput() override final;

};

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
void
FiniteElementDiscretization::
CalculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                 SUBREGION_TYPE * const elementSubRegion,
                                 FE_TYPE & finiteElement ) const
{
  GEOSX_MARK_FUNCTION;

  array4d< real64 > & dNdX = elementSubRegion->dNdX();
  array2d< real64 > & detJ = elementSubRegion->detJ();
  auto const & elemsToNodes = elementSubRegion->nodeList().toViewConst();

  string const elementTypeString = elementSubRegion->GetElementTypeString();

  constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
  constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
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

}


} /* namespace geosx */

#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_ */
