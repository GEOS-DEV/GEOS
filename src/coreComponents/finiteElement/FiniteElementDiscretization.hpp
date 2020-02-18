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


  std::unique_ptr<FiniteElementBase> getFiniteElement( string const & catalogName ) const;

  void ApplySpaceToTargetCells( ElementSubRegionBase * const group ) const;



  template< typename SUBREGION_TYPE >
  void CalculateShapeFunctionGradients( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X,
                                        SUBREGION_TYPE * const elementSubRegion ) const
  {
    GEOSX_MARK_FUNCTION;

    arrayView3d<R1Tensor> const & dNdX = elementSubRegion->template getReference< array3d< R1Tensor > >(dataRepository::keys::dNdX);
    arrayView2d<real64> const & detJ = elementSubRegion->template getReference< array2d<real64> >(dataRepository::keys::detJ);
    auto const & elemsToNodes = elementSubRegion->nodeList().toViewConst();

    PRAGMA_OMP( omp parallel )
    {
      std::unique_ptr<FiniteElementBase> fe = getFiniteElement( m_parentSpace );

      PRAGMA_OMP( omp for )
      for (localIndex k = 0 ; k < elementSubRegion->size() ; ++k)
      {
        fe->reinit(X, elemsToNodes[k]);

        for( localIndex q = 0 ; q < fe->n_quadrature_points() ; ++q )
        {
          detJ(k, q) = fe->JxW(q);
          for (localIndex b = 0 ; b < fe->dofs_per_element() ; ++b)
          {
            dNdX[k][q][b] = fe->gradient(b, q);
          }
        }
      }
    }
  }

  localIndex getNumberOfQuadraturePoints() const;

  string m_basisName;
  string m_quadratureName;
  string m_parentSpace;

  BasisBase const *    m_basis    = nullptr;
  QuadratureBase const * m_quadrature = nullptr;
  FiniteElementBase * m_finiteElement = nullptr;
protected:
  void PostProcessInput() override final;

};

} /* namespace geosx */

#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_ */
