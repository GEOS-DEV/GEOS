/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KValueFlashParameters.hpp
 */

#include "KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
KValueFlashParameters< NUM_PHASE >::KValueFlashParameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

template< integer NUM_PHASE >
std::unique_ptr< ModelParameters > KValueFlashParameters< NUM_PHASE >::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< KValueFlashParameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< KValueFlashParameters >( std::move( parameters ) );
}

template< integer NUM_PHASE >
void KValueFlashParameters< NUM_PHASE >::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::kValueTablesString(), &m_kValueTables ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "List of k-value tables for each phase." );
}

template< integer NUM_PHASE >
void KValueFlashParameters< NUM_PHASE >::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                      ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numFluidPhase = fluid->numFluidPhases();
  integer const numFluidComponent = fluid->numFluidComponents();

  GEOS_THROW_IF_NE_MSG( numPhases, numFluidPhase,
                        GEOS_FMT( "{}: invalid number of phases for the fluid.", fluid->getFullName() ),
                        InputError );

  integer const numPhaseTable = m_kValueTables.size( 0 );
  integer const numComponentTable = m_kValueTables.size( 1 );

  GEOS_THROW_IF_NE_MSG( numPhases-1, numPhaseTable,
                        GEOS_FMT( "{}: invalid number of phases provided for k-value tables. "
                                  "First dimension of table should be {}", fluid->getFullName(), numPhases-1 ),
                        InputError );

  GEOS_THROW_IF_NE_MSG( numFluidComponent, numComponentTable,
                        GEOS_FMT( "{}: invalid number of components provided for k-value tables. "
                                  "Second dimension of table should be {}", fluid->getFullName(), numFluidComponent ),
                        InputError );
}

// Instantiate
template class KValueFlashParameters< 2 >;
template class KValueFlashParameters< 3 >;

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
