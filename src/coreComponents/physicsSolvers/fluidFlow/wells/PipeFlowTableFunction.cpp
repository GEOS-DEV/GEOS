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
 * @file PipeFlowTableFunction.cpp
 */

#include "PipeFlowTableFunction.hpp"

#include "common/DataTypes.hpp"
#include <algorithm>

namespace geos
{

using namespace dataRepository;

PipeFlowTableFunction::PipeFlowTableFunction( const string & name,
                                              Group * const parent ):
  MultivariableTableFunction( name, parent )
{

  registerWrapper( viewKeyStruct::rateType(), &m_rateType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of rate entered in the rates array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::rateArray(), &m_rate ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of rates" );

  registerWrapper( viewKeyStruct::waterFractionType(), &m_waterFractionType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of water fraction  entered in the wfr array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::waterFractionArray(), &m_wfr ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of water fractions " );

  registerWrapper( viewKeyStruct::gasFractionType(), &m_gasFractionType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of gas fraction  entered in the wfr array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::gasFractionArray(), &m_gfr ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of gas fractions " );

  registerWrapper( viewKeyStruct::wellHeadPressureArray(), &m_whp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of well head pressures " );

  registerWrapper( viewKeyStruct::bottomHolePressureArray(), &m_bhp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of gas fractions " );

}

void PipeFlowTableFunction::initializeFunction()
{}
REGISTER_CATALOG_ENTRY( FunctionBase, PipeFlowTableFunction, string const &, Group * const )

} // end of namespace geos
