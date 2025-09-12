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
#include "functions/FunctionManager.hpp"
#include "functions/MultivariableNonuniformTableFunctionKernels.hpp"
#include "common/DataTypes.hpp"
#include <algorithm>

namespace geos
{

using namespace dataRepository;

PipeFlowTableFunction::PipeFlowTableFunction( const string & name,
                                              Group * const parent ):
  MultivariableNonuniformTableFunction( name, parent ),
  m_tableFunction( nullptr )
{
  registerWrapper( viewKeyStruct::tableName(), &m_tableName ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Flow table name. Table is associated with a Min/MaxWHPConstraint with the constraints flowTableName field" );

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
    setDescription( "Array of bottom hole pressures representing the dependent variable for the table function.\n"
                    "The quantities must be entered such that the rate index increases first, followed by whp, waterFraction, and finally gasFraction.\n"
                    " For example  bhp[1]= f(rate[1],whp[1],,waterFraction[1],gasFraction[1]) , bhp[2]= f(rate[2],whp[1],waterFraction[1],gasFraction[1])" );

}

void PipeFlowTableFunction::postInputInitialization()
{
  // Validate independent table values are not decreasing
  auto checkNotDecreasing = []( const auto & arr, const std::string & name )
  {
    for( auto i = 1; i < arr.size(); ++i )
    {
      if( arr[i] < arr[i-1] )
      {
        GEOS_ERROR( name << " array values must be non-decreasing, but arr[" << i << "] < arr[" << i-1 << "]" );
      }
    }
  };

  checkNotDecreasing( m_rate, getRateType());
  checkNotDecreasing( m_whp, "wellHeadPressure" );
  checkNotDecreasing( m_wfr, getWaterFractionType());
  checkNotDecreasing( m_gfr, getGasFractionType());

  // fluid model phase naming convention associated with rate type
  if( m_rateType == "liq" )
  {
    m_ratePhases.resize( 2 );
    m_ratePhases[0] = "oil";
    m_ratePhases[1] = "water";
  }
  else if( m_rateType == "oil" )
  {
    m_ratePhases.resize( 1 );
    m_ratePhases[0] = "oil";
  }
  else if( m_rateType == "gas" )
  {
    m_ratePhases.resize( 1 );
    m_ratePhases[0] = "gas";
  }

  initializeFunction();
}

void PipeFlowTableFunction::initializeFunction()
{
  localIndex constexpr nDims = 4;
  localIndex constexpr nOps = 1;

  FunctionManager * functionManager = &FunctionManager::getInstance();
  // Create nonuniform version of uniformly spaced table
  m_tableFunction = &(dynamicCast< MultivariableNonuniformTableFunction & >( *functionManager->createChild( MultivariableNonuniformTableFunction::catalogName(), getTableName()  ) ));

  // find max number if independent vars
  int maxVar = std::max( m_rate.size(), m_whp.size());
  maxVar = std::max( maxVar, m_wfr.size());
  maxVar = std::max( maxVar, m_gfr.size());

  // Copy inputs into formats needed by table function

  integer_array axisPoints( nDims );

  axisPoints[0] = m_whp.size();
  axisPoints[1] = m_gfr.size();
  axisPoints[2] = m_wfr.size();
  axisPoints[3] = m_rate.size();

  array2d< real64 > axisCoordinates( nDims, maxVar );

  for( int i=0; i<m_whp.size(); i++ )
  {
    axisCoordinates[0][i] = m_whp[i];
  }

  for( int i=0; i<m_gfr.size(); i++ )
  {
    axisCoordinates[1][i] = m_gfr[i];
  }
  for( int i=0; i<m_wfr.size(); i++ )
  {
    axisCoordinates[2][i] = m_wfr[i];
  }
  for( int i=0; i< m_rate.size(); i++ )
  {
    axisCoordinates[3][i] = m_rate[i];
  }
  setTableCoordinates( nDims, nOps, axisCoordinates, axisPoints );
  setTableValues( m_bhp );
  MultivariableNonuniformTableFunction::initializeFunction();


}
void PipeFlowTableFunction::calculateBHP( array1d< real64 > const & phaseRates, real64 const & whp, real64 & bhp, integer & solveStat ) const
{

  MultivariableNonuniformTableFunctionStaticKernel< 4, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );

  solveStat = 0;                                                                  // Assume success
  // liq(oil)=0 vap = 1 wat = 2
  real64 totalLiquedRate = (phaseRates[0] + phaseRates[2]);
  //for( int i = 0; i < phaseRates.size(); ++i )
  //{
  //  totalVolumeRate += phaseRates[i];
//  }
  std::cout << bhp << " " << phaseRates << " " << whp << std::endl;

  real64 wct = phaseRates[2]/totalLiquedRate;
  real64 gor = m_gfr[0]; //phaseRates[1]/phaseRates[0];
  real64 liq = (phaseRates[0] + phaseRates[1] );  // liquid rate
  array1d< real64 > table_coords( 4 );
  table_coords[0]=whp;  // liquid rate
  table_coords[1]=gor; // well head pressure
  table_coords[2]=wct; // water cut
  table_coords[3]=liq; // gas oil ratio
  std::cout << "PipeFlowTableFunction::calculateWHP input bhp = " << bhp << " liq = " << liq << " whp " << whp << " wct = " << wct << " gor = " << gor << std::endl;

  array1d< real64 > table_bhp( 1 );
  array2d< real64 > derivs( 1, 4 );
  kernel.compute( table_coords, table_bhp, derivs );
  bhp = table_bhp[0];
}

void PipeFlowTableFunction::calculateWHP( real64 const & bhp, array1d< real64 > const & phaseRates, real64 & whp, integer & solveStat ) const
{

  MultivariableNonuniformTableFunctionStaticKernel< 4, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );

  solveStat = 0;  // Assume success
  // liq(oil)=0 vap = 1 wat = 2
  //real64 totalVolumeRate = 0.0;
  //for( int i = 0; i < phaseRates.size(); ++i )
  //{
  //  totalVolumeRate += phaseRates[i];
//  }
  std::cout << bhp << " " << phaseRates << " " << whp << std::endl;
  integer m_sign=-1;
  real64 wct = 0; // phaseRates[2]/totalVolumeRate*m_sign;
  real64 gor = 0; //phaseRates[1]/phaseRates[0];
  real64 liq = (phaseRates[0] /* + phaseRates[1]*/)*m_sign;  // liquid rate
  array1d< real64 > table_coords( 4 );

  table_coords[0]=whp; // well head pressure
  table_coords[1]=gor; // gas oil ratio
  table_coords[2]=wct; // water cut
  table_coords[3]=liq;  // liquid rate
  std::cout << "PipeFlowTableFunction::calculateWHP input bhp = " << bhp << " liq = " << liq << " whp " << whp << " wct = " << wct << " gor = " << gor << std::endl;

  array1d< real64 > table_bhp( 1 );
  array2d< real64 > derivs( 1, 4 );
  kernel.compute( table_coords, table_bhp, derivs );
  std::cout << " residual bhp = " <<  table_bhp[0] - bhp << std::endl;
  integer const maxIters=20;
  real64 const tol = 1e-6;
  integer iter = 0;
  while( iter < maxIters && std::abs( table_bhp[0] - bhp ) > tol )
  {
    // update whp
    table_coords[0] -= ( table_bhp[0] - bhp ) / derivs[0][0];
    kernel.compute( table_coords, table_bhp, derivs );
    std::cout << " PipeFlowTableFunction::calculateWHP iter = " << iter << " whp = " << table_coords[0] << " residual = " << table_bhp[0] - bhp << std::endl;

    ++iter;
  }
  whp = table_coords[0];
}

void PipeFlowTableFunction::writeTable() const
{

  std::ofstream of;
  std::string filename = getTableName() + ".csv";
  of.open( filename );



  MultivariableNonuniformTableFunctionStaticKernel< 4, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );

  of << "wellHeadPressure,gasFraction,waterFraction,rate,bottomHolePressure" << std::endl;
  array1d< real64 > table_coords( 4 );

  for( integer j=0; j < m_whp.size(); ++j )
  {
    table_coords[0]=m_whp[j];   // well head pressure

    for( integer l=0; l < m_gfr.size(); ++l )
    {
      table_coords[1]=m_gfr[l];     // gas oil ratio
      for( integer k=0; k < m_wfr.size(); ++k )
      {
        table_coords[2]=m_wfr[k]; // water cut
        for( integer i=0; i < m_rate.size(); ++i )
        {
          table_coords[3]=m_rate[i]; // liquid rate

          array1d< real64 > table_bhp( 1 );
          array2d< real64 > derivs( 1, 4 );
          kernel.compute( table_coords, table_bhp, derivs );
          of << table_coords[0] << "," << table_coords[1] << "," << table_coords[2] << "," << table_coords[3] << "," << table_bhp[0] << std::endl;
        }
      }
    }
  }
  of.close();


}

REGISTER_CATALOG_ENTRY( FunctionBase, PipeFlowTableFunction, string const &, Group * const )

} // end of namespace geos
