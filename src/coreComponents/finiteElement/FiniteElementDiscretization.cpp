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
{}

//localIndex FiniteElementDiscretization::getNumberOfQuadraturePoints() const
//{
//  return m_quadrature->size();
//}

//std::unique_ptr< FiniteElementBase > FiniteElementDiscretization::getFiniteElement( string const & ) const
//{
//  return FiniteElementBase::CatalogInterface::Factory( m_parentSpace,
//                                                       *m_basis,
//                                                       *m_quadrature,
//                                                       0 );
//}

void FiniteElementDiscretization::PostProcessInput()
{
//  auto const & basisName = this->getReference< string >( keys::basis );
//  auto const & quadratureName = this->getReference< string >( keys::quadrature );

  // TODO find a better way to do this that doesn't involve getParent(). We
  // shouldn't really use that unless there is no
  // other choice.
//  NumericalMethodsManager const & numericalMethods = *(this->getParent()->getParent()->group_cast< NumericalMethodsManager const * >());
//  Group const & basisManager = numericalMethods.getBasisFunctions();
//  Group const & quadratureManager = numericalMethods.getQuadratureRules();

//  m_basis = basisManager.GetGroup< BasisBase >( basisName );
//  m_quadrature = quadratureManager.GetGroup< QuadratureBase >( quadratureName );
}



REGISTER_CATALOG_ENTRY( Group, FiniteElementDiscretization, std::string const &, Group * const )

} /* namespace geosx */
