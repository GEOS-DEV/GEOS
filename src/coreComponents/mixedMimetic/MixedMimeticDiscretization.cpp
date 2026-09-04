/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MixedMimeticDiscretization.cpp
 */

#include "MixedMimeticDiscretization.hpp"

#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"

namespace geos
{

using namespace dataRepository;
using namespace mimeticInnerProduct;

MixedMimeticDiscretization::MixedMimeticDiscretization( string const & name,
                                                        Group * const parent )
  : Group( name, parent ),
  m_isAdaptive( 1 ),
  m_residualTolerance( 1e-3 ),
  m_nominalGradient( { 1.0, 1.0, 1.0 } )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::innerProductTypeString(), &m_innerProductType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Type of inner product used in the MFD-compatible cells of the mixed mimetic solver" );

  registerWrapper( viewKeyStruct::adaptiveString(), &m_isAdaptive ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag enabling the residual-based Global Adaptation: when enabled (1, default), the "
                    "cell-wise inner product is selected between TPFA and innerProductType according to the "
                    "consistency indicator; when disabled (0), innerProductType is used in every cell" );

  registerWrapper( viewKeyStruct::residualToleranceString(), &m_residualTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e-3 ).
    setDescription( "Tolerance on the Global Adaptation residual indicator used in the cell marking criterion" );

  registerWrapper( viewKeyStruct::nominalGradientString(), &m_nominalGradient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_nominalGradient ).
    setDescription( "Nominal pressure gradient inducing the projected admissible flow field used by the residual indicators" );
}

void MixedMimeticDiscretization::postInputInitialization()
{
  Group::postInputInitialization();

  GEOS_THROW_IF_LT_MSG( m_residualTolerance, 0.0,
                        GEOS_FMT( "{}: the residual tolerance cannot be negative",
                                  getWrapperDataContext( viewKeyStruct::residualToleranceString() ) ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( LvArray::tensorOps::l2Norm< 3 >( m_nominalGradient.data ), LvArray::NumericLimits< real64 >::epsilon,
                        GEOS_FMT( "{}: the nominal gradient cannot be the zero vector",
                                  getWrapperDataContext( viewKeyStruct::nominalGradientString() ) ),
                        InputError );

  // degenerate combination: with a TPFA inner product, the adaptive blend mixes TPFA with TPFA
  // and the operator is TPFA everywhere regardless of the marking; the consistency indicator
  // fields are still computed and output, which makes this combination a useful diagnostic mode
  if( isAdaptive() && m_innerProductType == mimeticInnerProduct::MimeticInnerProductTypeStrings::TPFA )
  {
    GEOS_WARNING( GEOS_FMT( "{}: 'adaptive' is enabled but 'innerProductType' is TPFA: the adaptation has no "
                            "effect on the discretization (both operators of the adaptive blend coincide), and "
                            "the scheme reduces to full TPFA. The consistency indicators are still computed and "
                            "output. Select an MFD inner product (e.g. quasiTPFA) to activate the adaptation.",
                            getDataContext() ) );
  }
}

void MixedMimeticDiscretization::initializePostInitialConditionsPreSubGroups()
{
  Group::initializePostInitialConditionsPreSubGroups();

  std::unique_ptr< MimeticInnerProductBase > newMimeticIP = factory( m_innerProductType );

  registerWrapper< MimeticInnerProductBase >( viewKeyStruct::innerProductString(), std::move( newMimeticIP ) ).
    setRestartFlags( dataRepository::RestartFlags::NO_WRITE );
}

bool MixedMimeticDiscretization::isTpfaInnerProduct() const
{
  return m_innerProductType == mimeticInnerProduct::MimeticInnerProductTypeStrings::TPFA;
}

MixedMimeticDiscretization::CatalogInterface::CatalogType &
MixedMimeticDiscretization::getCatalog()
{
  static MixedMimeticDiscretization::CatalogInterface::CatalogType catalog;
  return catalog;
}

std::unique_ptr< MimeticInnerProductBase >
MixedMimeticDiscretization::factory( string const & mimeticInnerProductType ) const
{
  std::unique_ptr< MimeticInnerProductBase > rval;
  if( mimeticInnerProductType == MimeticInnerProductTypeStrings::TPFA )
  {
    rval = std::make_unique< TPFAInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::QuasiTPFA )
  {
    rval = std::make_unique< QuasiTPFAInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::Simple )
  {
    rval = std::make_unique< SimpleInnerProduct >();
  }
  else if( mimeticInnerProductType == MimeticInnerProductTypeStrings::BdVLM )
  {
    rval = std::make_unique< BdVLMInnerProduct >();
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Key value of {} does not have an associated mimetic inner product implementing the mixed-form "
                          "mass matrix (valid options: {}, {}, {}, {}).",
                          mimeticInnerProductType,
                          MimeticInnerProductTypeStrings::TPFA,
                          MimeticInnerProductTypeStrings::QuasiTPFA,
                          MimeticInnerProductTypeStrings::Simple,
                          MimeticInnerProductTypeStrings::BdVLM ),
                getDataContext() );
  }
  return rval;
}

REGISTER_CATALOG_ENTRY( MixedMimeticDiscretization, MixedMimeticDiscretization, string const &, Group * const )

}
