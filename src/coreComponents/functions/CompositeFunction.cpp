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

#include "FunctionManager.hpp"
#include "CompositeFunction.hpp"

namespace geos
{

using namespace dataRepository;


CompositeFunction::CompositeFunction( const string & name,
                                      Group * const parent ):
  FunctionBase( name, parent ),
  m_numSubFunctions(),
  m_subFunctions()
{
  registerWrapper( viewKeyStruct::functionNamesString(), &m_functionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of source functions. The order must match the variableNames argument." );

  registerWrapper( viewKeyStruct::operationTypeString(), &m_operationType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Operation to apply to the functions. Valid options:\n* " + EnumStrings< OperationType >::concat( "\n* " ) );
}

CompositeFunction::~CompositeFunction()
{}

void CompositeFunction::initializeFunction()
{
  // Grab pointers to sub functions
  FunctionManager & functionManager = FunctionManager::getInstance();
  m_numSubFunctions = LvArray::integerConversion< localIndex >( m_functionNames.size());
  for( localIndex ii=0; ii<m_numSubFunctions; ++ii )
  {
    m_subFunctions.emplace_back( &functionManager.getGroup< FunctionBase >( m_functionNames[ii] ) );
  }
}

 // Function to dynamically set the operation based on a string
void CompositeFunction::setOperation() 
{
  switch( m_operationType )
  {
    case OperationType::sum:
      m_operation = math::Sum{};
      break;
    case OperationType::product:
      m_operation = math::Product{};
      break;
    case OperationType::difference:
      m_operation = math::Difference{};
      break;
    case OperationType::division:
      m_operation = math::Division{};
      break;
    default:
      GEOS_THROW( GEOS_FMT("Unsupported operation {}", m_operationType ), std::invalid_argument );
  }
}

void CompositeFunction::evaluate( dataRepository::Group const & group,
                                  real64 const time,
                                  SortedArrayView< localIndex const > const & set,
                                  arrayView1d< real64 > const & result ) const
{
  // Evaluate each of the subFunctions independently and place the results into
  // a temporary field
  array1d< real64_array > subFunctionResults;
  for( localIndex ii=0; ii<m_numSubFunctions; ++ii )
  {
    real64_array tmp( result.size());
    m_subFunctions[ii]->evaluate( group, time, set, tmp );
    subFunctionResults.emplace_back( std::move( tmp ));
  }

  // Evaluate the symbolic math
  forAll< serialPolicy >( set.size(), [&, result, set]( localIndex const i )
  {
    localIndex const ii = set[ i ];
    std::vector<real64> functionResults(m_numSubFunctions);
    for( localIndex jj=0; jj<m_numSubFunctions; ++jj )
    {
      functionResults[jj] = subFunctionResults[jj][ii];
    }

    result[ii] = math::apply( m_operation, functionResults );
  } );
}

real64 CompositeFunction::evaluate( real64 const * const input ) const
{
  std::vector<real64> functionResults(m_numSubFunctions);

  for( localIndex ii=0; ii<m_numSubFunctions; ++ii )
  {
    functionResults[ii] = m_subFunctions[ii]->evaluate( input );
  }

  return math::apply( m_operation, functionResults );
}

REGISTER_CATALOG_ENTRY( FunctionBase, CompositeFunction, string const &, Group * const )

} // namespace geos
