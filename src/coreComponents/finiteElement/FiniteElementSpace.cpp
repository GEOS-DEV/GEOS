/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
#include "managers/NumericalMethodsManager.hpp"
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
{
  delete m_finiteElement;
}


void FiniteElementSpace::BuildDataStructure( dataRepository::ManagedGroup * const parent )
{}

void FiniteElementSpace::FillDocumentationNode()
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

  //Ensure data is contiguous
  auto dNdXView       = cellBlock->RegisterViewWrapper< multidimensionalArray::ManagedArray< R1Tensor, 3 > >(keys::dNdX);
  dNdXView->setSizedFromParent(1);
  multidimensionalArray::ManagedArray< R1Tensor, 3 > &  dNdX = dNdXView->reference();
  dNdX.resize( cellBlock->size(), m_quadrature->size(), m_finiteElement->dofs_per_element() );

  auto & constitutiveMapView = *(cellBlock->getWrapper< std::pair< array2d< localIndex >, array2d< localIndex > > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString));
  constitutiveMapView.setSizedFromParent(1);
  auto & constitutiveMap = constitutiveMapView.reference();
  constitutiveMap.first.resize(cellBlock->size(), m_quadrature->size() );
  constitutiveMap.second.resize(cellBlock->size(), m_quadrature->size() );

  auto & detJView = *(cellBlock->RegisterViewWrapper< array2d< real64 > >(keys::detJ));
  detJView.setSizedFromParent(1);
  auto & detJ = detJView.reference();
  detJ.resize(cellBlock->size(), m_quadrature->size() );


}

void FiniteElementSpace::CalculateShapeFunctionGradients( r1_array const &  X,
                                                          dataRepository::ManagedGroup * const cellBlock ) const
{
  auto & dNdX            = cellBlock->getReference< multidimensionalArray::ManagedArray< R1Tensor, 3 > >(keys::dNdX);
  auto & detJ            = cellBlock->getReference< array2d<real64> >(keys::detJ);
  FixedOneToManyRelation const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(std::string("nodeList"))->reference();// getData<array2d<localIndex>>(keys::nodeList);

  array1d<R1Tensor> X_elemLocal( m_finiteElement->dofs_per_element() );


  for (localIndex k = 0 ; k < cellBlock->size() ; ++k)
  {
    localIndex const * const elemToNodeMap = elemsToNodes[k];

    CopyGlobalToLocal(elemToNodeMap, X, X_elemLocal);

    m_finiteElement->reinit(X_elemLocal);

    for( localIndex q = 0 ; q < m_finiteElement->n_quadrature_points() ; ++q )
    {

      detJ(k, q) = m_finiteElement->JxW(q);
      for (localIndex b = 0 ; b < m_finiteElement->dofs_per_element() ; ++b)
      {
        dNdX[k][q][b] =  m_finiteElement->gradient(b, q);
//        std::cout<<"dNdX["<<k<<"]["<<q<<"]["<<b<<"] :"<<dNdX[k][q][b]<<std::endl;
      }

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
  ManagedGroup const *  numericalMethods = this->getParent()->getParent();
  ManagedGroup const *  basisManager = numericalMethods->GetGroup(keys::basisFunctions);
  ManagedGroup const *  quadratureManager = numericalMethods->GetGroup(keys::quadratureRules);
  
  m_basis = basisManager->getData<BasisBase>(basisName);
  m_quadrature = quadratureManager->getData<QuadratureBase>(quadratureName);
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
