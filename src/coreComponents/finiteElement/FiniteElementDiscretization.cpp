/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

#include "FiniteElementDiscretization.hpp"

#include "../mesh/CellElementSubRegion.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mesh/NodeManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "basis/BasisBase.hpp"
#include "quadrature/QuadratureBase.hpp"
#include "ElementLibrary/FiniteElement.h"
#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"

// TODO make this not dependent on this header...need better key implementation

namespace geosx
{
using namespace dataRepository;



FiniteElementDiscretization::FiniteElementDiscretization( std::string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);

  RegisterViewWrapper( keys::basis, &m_basisName, false )->setInputFlag(InputFlags::REQUIRED);
  RegisterViewWrapper( keys::quadrature, &m_quadratureName, false )->setInputFlag(InputFlags::REQUIRED);
}

FiniteElementDiscretization::~FiniteElementDiscretization()
{
  delete m_finiteElement;
}

localIndex FiniteElementDiscretization::getNumberOfQuadraturePoints() const
{
  return m_quadrature->size();
}

std::unique_ptr<FiniteElementBase> FiniteElementDiscretization::getFiniteElement( string const & catalogName ) const
{
  return FiniteElementBase::CatalogInterface::Factory( catalogName,
                                                       *m_basis,
                                                       *m_quadrature,
                                                       0 );
}

void FiniteElementDiscretization::ApplySpaceToTargetCells( ElementSubRegionBase * const cellBlock ) const
{
  GEOSX_MARK_FUNCTION;

  // TODO THis crap needs to get cleaned up and worked out in the data structure
  // much better than this.
  // Need to provide some mechanism to set the sizedFromParent during the
  // registration, or only allow documentation node
  // registration.

  std::unique_ptr<FiniteElementBase> fe = getFiniteElement( cellBlock->GetElementTypeString() );

  // dNdX holds a lot of POD data and it gets set in the method below so there's no need to zero initialize it.
  array3d< R1Tensor > &  dNdX = cellBlock->RegisterViewWrapper< array3d< R1Tensor > >(keys::dNdX)->reference();
  dNdX.resizeWithoutInitializationOrDestruction( cellBlock->size(), m_quadrature->size(), fe->dofs_per_element() );

  auto & constitutiveMap = cellBlock->getWrapper< std::pair< array2d< localIndex >, array2d< localIndex > > >(CellElementSubRegion::viewKeyStruct::constitutiveMapString)->reference();
  constitutiveMap.first.resize(cellBlock->size(), m_quadrature->size() );
  constitutiveMap.second.resize(cellBlock->size(), m_quadrature->size() );

  array2d< real64 > & detJ = cellBlock->RegisterViewWrapper< array2d< real64 > >(keys::detJ)->reference();
  detJ.resize(cellBlock->size(), m_quadrature->size() );
}

void FiniteElementDiscretization::CalculateShapeFunctionGradients( arrayView1d<R1Tensor const> const & X,
                                                                   ElementSubRegionBase * const elementSubRegion ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView3d<R1Tensor> const & dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);
  arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);
  FixedOneToManyRelation const & elemsToNodes = elementSubRegion->getWrapper<FixedOneToManyRelation>(std::string("nodeList"))->reference();

  PRAGMA_OMP( omp parallel )
  {
    std::unique_ptr<FiniteElementBase> fe = getFiniteElement( elementSubRegion->GetElementTypeString() );

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

void FiniteElementDiscretization::PostProcessInput()
{
  auto const & basisName = this->getReference<string>(keys::basis);
  auto const & quadratureName = this->getReference<string>(keys::quadrature);

  // TODO find a better way to do this that doesn't involve getParent(). We
  // shouldn't really use that unless there is no
  // other choice.
  ManagedGroup const *  numericalMethods = this->getParent()->getParent();
  ManagedGroup const *  basisManager = numericalMethods->GetGroup(keys::basisFunctions);
  ManagedGroup const *  quadratureManager = numericalMethods->GetGroup(keys::quadratureRules);
  
  m_basis = basisManager->GetGroup<BasisBase>(basisName);
  m_quadrature = quadratureManager->GetGroup<QuadratureBase>(quadratureName);
  m_finiteElement = new FiniteElement<3>( *m_basis, *m_quadrature, 0);
}



REGISTER_CATALOG_ENTRY( ManagedGroup, FiniteElementDiscretization, std::string const &, ManagedGroup * const )

} /* namespace geosx */
