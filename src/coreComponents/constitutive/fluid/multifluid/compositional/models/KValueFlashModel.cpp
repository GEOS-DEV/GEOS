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
 * @file KValueFlashModel.cpp
 */

#include "KValueFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.hpp"

#include "functions/FunctionManager.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< >
string KValueFlashModel< 2 >::catalogName()
{
  return "TwoPhaseKValue";
}

template< >
string KValueFlashModel< 3 >::catalogName()
{
  return "ThreePhaseKValue";
}

template< integer NUM_PHASE >
KValueFlashModel< NUM_PHASE >::KValueFlashModel( string const & name,
                                                 ComponentProperties const & componentProperties,
                                                 ModelParameters const & modelParameters,
                                                 arrayView1d< integer const > const phaseTypes ):
  FunctionBase( name, componentProperties ),
  m_phaseTypes( phaseTypes )
{
  m_parameters = modelParameters.get< KValueFlashParameters< NUM_PHASE > >();
  integer const numComps = m_componentProperties.getNumberOfComponents();
  m_presentComponents.resize( numComps );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    m_presentComponents[ic] = ic;
  }
}

template< integer NUM_PHASE >
typename KValueFlashModel< NUM_PHASE >::KernelWrapper
KValueFlashModel< NUM_PHASE >::createKernelWrapper() const
{
  array1d< integer > phaseOrder( KernelWrapper::getNumberOfPhases());
  calculatePhaseOrdering( m_phaseTypes, phaseOrder );
  integer const liquidIndex = phaseOrder[0];
  integer const vapourIndex = phaseOrder[1];
  integer const aqueousIndex = (NUM_PHASE == 3) ? phaseOrder[2] : -1;

  string pressureTableName;
  string temperatureTableName;
  m_parameters->createTables( functionName(),
                              pressureTableName,
                              temperatureTableName );


  FunctionManager const & functionManager = FunctionManager::getInstance();
  TableFunction const * pressureTable = functionManager.getGroupPointer< TableFunction >( pressureTableName );
  TableFunction const * temperatureTable = functionManager.getGroupPointer< TableFunction >( temperatureTableName );

  return KernelWrapper( m_componentProperties.getNumberOfComponents(),
                        liquidIndex,
                        vapourIndex,
                        aqueousIndex,
                        *pressureTable,
                        *temperatureTable,
                        m_presentComponents,
                        m_parameters->m_kValueHyperCube );
}

template< integer NUM_PHASE >
KValueFlashModelUpdate< NUM_PHASE >::KValueFlashModelUpdate( integer const numComponents,
                                                             integer const liquidIndex,
                                                             integer const vapourIndex,
                                                             integer const aqueousIndex,
                                                             TableFunction const & pressureTable,
                                                             TableFunction const & temperatureTable,
                                                             arrayView1d< integer const > const & presentComponents,
                                                             arrayView4d< real64 const > const & kValues ):
  m_liquidIndex( liquidIndex ),
  m_vapourIndex( vapourIndex ),
  m_aquoesIndex( aqueousIndex ),
  m_numComponents( numComponents ),
  m_numPressurePoints( pressureTable.getCoordinates()[0].size() ),
  m_numTemperaturePoints( temperatureTable.getCoordinates()[0].size() ),
  m_pressureTable( pressureTable.createKernelWrapper()),
  m_temperatureTable( temperatureTable.createKernelWrapper()),
  m_presentComponents( presentComponents ),
  m_kValues( kValues )
{}

template< integer NUM_PHASE >
std::unique_ptr< ModelParameters >
KValueFlashModel< NUM_PHASE >::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return KValueFlashParameters< NUM_PHASE >::create( std::move( parameters ) );
}

template< integer NUM_PHASE >
void KValueFlashModel< NUM_PHASE >::calculatePhaseOrdering( arrayView1d< integer const > const & phaseTypes,
                                                            arrayView1d< integer > const & phaseOrder )
{
  if constexpr (NUM_PHASE == 2)
  {
    NegativeTwoPhaseFlashModel::calculatePhaseOrdering( phaseTypes, phaseOrder );
  }
  else if constexpr (NUM_PHASE == 3)
  {
    ImmiscibleWaterFlashModel::calculatePhaseOrdering( phaseTypes, phaseOrder );
  }
}

// Instantiate
template class KValueFlashModel< 2 >;
template class KValueFlashModel< 3 >;

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
