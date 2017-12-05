/*
 * FiniteElementSpace.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#include "FiniteElementSpace.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mesh/NodeManager.hpp"
#include "FiniteElementManager.hpp"
#include "basis/BasisBase.hpp"
#include "quadrature/QuadratureBase.hpp"
#include "ElementLibrary/FiniteElement.h"
#include "codingUtilities/Utilities.hpp"

// TODO make this not dependent on this header...need better key implementation
#include "mesh/CellBlockSubRegion.hpp"

namespace geosx
{
using namespace dataRepository;



FiniteElementSpace::FiniteElementSpace( std::string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{}

FiniteElementSpace::~FiniteElementSpace()
{}


void FiniteElementSpace::BuildDataStructure( dataRepository::ManagedGroup * const parent )
{}

void FiniteElementSpace::FillDocumentationNode( dataRepository::ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();


  docNode->AllocateChildNode( keys::basis,
                              keys::basis,
                              -1,
                              "string",
                              "string",
                              "name of basis function object.",
                              "name of the basis function object.",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );
  docNode->AllocateChildNode( keys::quadrature,
                              keys::quadrature,
                              -1,
                              "string",
                              "string",
                              "name of basis function object.",
                              "name of the basis function object.",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );



}

void FiniteElementSpace::ApplySpaceToTargetCells( dataRepository::ManagedGroup * const cellBlock ) const
{

  // TODO THis crap needs to get cleaned up and worked out in the data structure
  // much better than this.
  // Need to provide some mechanism to set the sizedFromParent during the
  // registration, or only allow documentation node
  // registration.


  auto dNdXView        = cellBlock->RegisterViewWrapper< array< Array2dT<R1Tensor> > >(keys::dNdX);
  dNdXView->setSizedFromParent(1);
  dNdXView->resize();
  auto & dNdX            = dNdXView->reference();

  for( auto & entry : dNdX )
  {
    entry.resize( m_finiteElement->dofs_per_element(), m_quadrature->size() );
  }



  auto & constitutiveMapView = *(cellBlock->RegisterViewWrapper< std::pair< Array2dT< localIndex >, Array2dT< localIndex > > >(keys::constitutiveMap) );
  constitutiveMapView.setSizedFromParent(1);
  auto & constitutiveMap = constitutiveMapView.reference();
  constitutiveMap.first.resize(cellBlock->size(), m_quadrature->size() );
  constitutiveMap.second.resize(cellBlock->size(), m_quadrature->size() );



  auto & detJView = *(cellBlock->RegisterViewWrapper< Array2dT< real64 > >(keys::detJ));
  detJView.setSizedFromParent(1);
  auto & detJ = detJView.reference();
  detJ.resize(cellBlock->size(), m_quadrature->size() );


}

void FiniteElementSpace::CalculateShapeFunctionGradients( r1_array const &  X,
                                                          dataRepository::ManagedGroup * const cellBlock ) const
{
  auto & dNdX            = cellBlock->getReference< array< Array2dT<R1Tensor> > >(keys::dNdX);
  auto & detJ            = cellBlock->getReference< Array2dT<real64> >(keys::detJ);
  lArray2d const & elemsToNodes = cellBlock->getWrapper<lArray2d>(std::string("nodeList"))->reference();// getData<lArray2d>(keys::nodeList);

  array<R1Tensor> X_elemLocal( m_finiteElement->dofs_per_element() );


  for (localIndex k = 0 ; k < cellBlock->size() ; ++k)
  {
    const localIndex* const elemToNodeMap = elemsToNodes[k];

    CopyGlobalToLocal(elemToNodeMap, X, X_elemLocal);

    m_finiteElement->reinit(X_elemLocal);

    for( localIndex q = 0 ; q < m_finiteElement->n_quadrature_points() ; ++q )
    {

      detJ(k, q) = m_finiteElement->JxW(q);
      for (localIndex b = 0 ; b < m_finiteElement->dofs_per_element() ; ++b)
      {
        dNdX(k)(q, b) = m_finiteElement->gradient(b, q);
      }
      //std::cout<<"Element, ip, dNdX :"<<k<<", "<<a<<",
      // "<<m_dNdX(k)(a)[0]<<std::endl;

    }
  }
}

void FiniteElementSpace::ReadXML_PostProcess()
{
  auto const & basisName = this->getData<string>(keys::basis);
  auto const & quadratureName = this->getData<string>(keys::quadrature);

  // TODO find a better way to do this that doesn't involve getParent(). We
  // shouldn't really use that unless there is no
  // other choice.
  m_basis = this->getParent()->GetGroup(keys::basisFunctions)->getData<BasisBase>(basisName);
  m_quadrature = this->getParent()->GetGroup(keys::quadratureRules)->getData<QuadratureBase>(quadratureName);
  m_finiteElement = new FiniteElement<3>( *m_basis, *m_quadrature, 0);

}

void FiniteElementSpace::InitializePreSubGroups( ManagedGroup * const group )
{
//  auto const & basisName = this->getData<string>(keys::basis) ;
//  auto const & quadratureName = this->getData<string>(keys::quadrature) ;
//
//  // TODO find a better way to do this that doesn't involve getParent(). We
// shouldn't really use that unless there is no
//  // other choice.
//  m_basis =
// this->getParent()->GetGroup(keys::basisFunctions)->getData<BasisBase>(basisName);
//  m_quadrature =
// this->getParent()->GetGroup(keys::quadratureRules)->getData<QuadratureBase>(quadratureName);
//  m_finiteElement = new FiniteElement<3>( *m_basis, *m_quadrature, 0);
}


REGISTER_CATALOG_ENTRY( ManagedGroup, FiniteElementSpace, std::string const &, ManagedGroup * const )

} /* namespace geosx */
