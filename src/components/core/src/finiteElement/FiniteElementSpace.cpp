// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
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

#include "../../../PhysicsSolverPackage1/src/SolidMechanicsLagrangianFEM-MiniApp/Layout.hpp"

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

  auto dNdXView        = cellBlock->RegisterViewWrapper< array< Array2dT<R1Tensor> > >(keys::dNdX);
  dNdXView->setSizedFromParent(1);
  static_cast< ViewWrapperBase * >(dNdXView)->resize();
  auto & dNdX            = dNdXView->reference();

  for( auto & entry : dNdX )
  {
    entry.resize( m_finiteElement->dofs_per_element(), m_quadrature->size() );
  }

  auto & constitutiveMapView = *(cellBlock->getWrapper< std::pair< Array2dT< localIndex >, Array2dT< localIndex > > >(CellBlockSubRegion::viewKeyStruct::constitutiveMapString));
  constitutiveMapView.setSizedFromParent(1);
  auto & constitutiveMap = constitutiveMapView.reference();
  constitutiveMap.first.resize(cellBlock->size(), m_quadrature->size() );
  constitutiveMap.second.resize(cellBlock->size(), m_quadrature->size() );

#if !defined(EXTERNAL_KERNELS) || defined(ARRAY_OF_OBJECTS_LAYOUT)
  auto dNdXView_all       = cellBlock->RegisterViewWrapper< multidimensionalArray::ManagedArray< real64, 4 > >("dNdX_all");  
  dNdXView_all->setSizedFromParent(1);
  multidimensionalArray::ManagedArray< real64, 4 > &  dNdX_all = dNdXView_all->reference();
  dNdX_all.resize( cellBlock->size(), m_quadrature->size(), m_finiteElement->dofs_per_element(), 3);
#endif  
  
#if defined(EXTERNAL_KERNELS) && defined(OBJECT_OF_ARRAYS_LAYOUT)
  auto dNdXView_x        = cellBlock->RegisterViewWrapper< multidimensionalArray::ManagedArray< real64, 3 > >("dNdX_x");
  auto dNdXView_y        = cellBlock->RegisterViewWrapper< multidimensionalArray::ManagedArray< real64, 3 > >("dNdX_y");
  auto dNdXView_z        = cellBlock->RegisterViewWrapper<multidimensionalArray:: ManagedArray< real64, 3 > >("dNdX_z");
  
  dNdXView_x->setSizedFromParent(1);
  dNdXView_y->setSizedFromParent(1);
  dNdXView_z->setSizedFromParent(1);
  
  multidimensionalArray::ManagedArray< real64, 3 > &  dNdX_x = dNdXView_x->reference();
  multidimensionalArray::ManagedArray< real64, 3 > &  dNdX_y = dNdXView_y->reference();
  multidimensionalArray::ManagedArray< real64, 3 > &  dNdX_z = dNdXView_z->reference();

  dNdX_x.resize( cellBlock->size(), m_quadrature->size(), m_finiteElement->dofs_per_element() );
  dNdX_y.resize( cellBlock->size(), m_quadrature->size(), m_finiteElement->dofs_per_element() );
  dNdX_z.resize( cellBlock->size(), m_quadrature->size(), m_finiteElement->dofs_per_element() );
#endif
  
  auto & detJView = *(cellBlock->RegisterViewWrapper< Array2dT< real64 > >(keys::detJ));
  detJView.setSizedFromParent(1);
  auto & detJ = detJView.reference();
  detJ.resize(cellBlock->size(), m_quadrature->size() );


}

void FiniteElementSpace::CalculateShapeFunctionGradients( r1_array const &  X,
                                                          dataRepository::ManagedGroup * const cellBlock ) const
{
  auto & dNdX            = cellBlock->getReference< array< Array2dT<R1Tensor> > >(keys::dNdX);

#if !defined(EXTERNAL_KERNELS) || defined(ARRAY_OBJECTS_LAYOUT)  
  auto & dNdX_all          = cellBlock->getReference< multidimensionalArray::ManagedArray<real64, 4> >("dNdX_all");
#endif  

#if defined(EXTERNAL_KERNELS) && defined(OBJECT_OF_ARRAYS_LAYOUT)
  auto & dNdX_x          = cellBlock->getReference< multidimensionalArray::ManagedArray<real64, 3> >("dNdX_x");
  auto & dNdX_y          = cellBlock->getReference< multidimensionalArray::ManagedArray<real64, 3> >("dNdX_y");
  auto & dNdX_z          = cellBlock->getReference< multidimensionalArray::ManagedArray<real64, 3> >("dNdX_z");
#endif  
  
  auto & detJ            = cellBlock->getReference< Array2dT<real64> >(keys::detJ);
  FixedOneToManyRelation const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(std::string("nodeList"))->reference();// getData<lArray2d>(keys::nodeList);

  array<R1Tensor> X_elemLocal( m_finiteElement->dofs_per_element() );


  for (localIndex k = 0 ; k < cellBlock->size() ; ++k)
  {
    arrayView1d<localIndex const> const elemToNodeMap = elemsToNodes[k];

    CopyGlobalToLocal(elemToNodeMap, X, X_elemLocal);

    m_finiteElement->reinit(X_elemLocal);

    for( localIndex q = 0 ; q < m_finiteElement->n_quadrature_points() ; ++q )
    {

      detJ(k, q) = m_finiteElement->JxW(q);
      for (localIndex b = 0 ; b < m_finiteElement->dofs_per_element() ; ++b)
      {
        R1Tensor temp = m_finiteElement->gradient(b, q);
        
        dNdX(k)(q, b) = temp;
        
#if !defined(EXTERNAL_KERNELS) || defined(ARRAY_OBJECTS_LAYOUT)        
        dNdX_all[k][q][b][0] = temp[0];
        dNdX_all[k][q][b][1] = temp[1];
        dNdX_all[k][q][b][2] = temp[2];
#endif        

#if defined(EXTERNAL_KERNELS) && defined(OBJECT_OF_ARRAYS_LAYOUT)        
        dNdX_x[k][q][b] = temp[0];
        dNdX_y[k][q][b] = temp[1];
        dNdX_z[k][q][b] = temp[2];
#endif        

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
