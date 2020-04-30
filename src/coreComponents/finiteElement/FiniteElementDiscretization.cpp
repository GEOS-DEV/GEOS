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
 * @file FiniteElementSpace.cpp
 */

#include "FiniteElementDiscretization.hpp"

#include "mesh/CellElementSubRegion.hpp"
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



FiniteElementDiscretization::FiniteElementDiscretization( std::string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( keys::basis, &m_basisName )->setInputFlag( InputFlags::REQUIRED );
  registerWrapper( keys::quadrature, &m_quadratureName )->setInputFlag( InputFlags::REQUIRED );
  registerWrapper( keys::parentSpace, &m_parentSpace )->setInputFlag( InputFlags::REQUIRED );
}

FiniteElementDiscretization::~FiniteElementDiscretization()
{
  delete m_finiteElement;
}

localIndex FiniteElementDiscretization::getNumberOfQuadraturePoints() const
{
  return m_quadrature->size();
}

std::unique_ptr< FiniteElementBase > FiniteElementDiscretization::getFiniteElement( string const & ) const
{
  return FiniteElementBase::CatalogInterface::Factory( m_parentSpace,
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

  //TODO: wu40: Temporarily use the parent space (read from xml) to assign element type for finite element calculation
  // (for C3D6 mesh).
  //Need to do this in a more natural way.
  std::unique_ptr< FiniteElementBase > fe = getFiniteElement( m_parentSpace );

  // dNdX holds a lot of POD data and it gets set in the method below so there's no need to zero initialize it.
  array3d< R1Tensor > & dNdX = cellBlock->registerWrapper< array3d< R1Tensor > >( keys::dNdX )->reference();
  dNdX.resizeWithoutInitializationOrDestruction( cellBlock->size(), m_quadrature->size(), fe->dofs_per_element() );

  array2d< real64 > & detJ = cellBlock->registerWrapper< array2d< real64 > >( keys::detJ )->reference();
  detJ.resize( cellBlock->size(), m_quadrature->size() );
}

void FiniteElementDiscretization::PostProcessInput()
{
  auto const & basisName = this->getReference< string >( keys::basis );
  auto const & quadratureName = this->getReference< string >( keys::quadrature );

  // TODO find a better way to do this that doesn't involve getParent(). We
  // shouldn't really use that unless there is no
  // other choice.
  Group const *  numericalMethods = this->getParent()->getParent();
  Group const *  basisManager = numericalMethods->GetGroup( keys::basisFunctions );
  Group const *  quadratureManager = numericalMethods->GetGroup( keys::quadratureRules );

  m_basis = basisManager->GetGroup< BasisBase >( basisName );
  m_quadrature = quadratureManager->GetGroup< QuadratureBase >( quadratureName );
  m_finiteElement = new FiniteElement< 3 >( *m_basis, *m_quadrature, 0 );
}



REGISTER_CATALOG_ENTRY( Group, FiniteElementDiscretization, std::string const &, Group * const )

} /* namespace geosx */
