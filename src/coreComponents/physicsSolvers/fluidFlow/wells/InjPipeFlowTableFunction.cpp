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
 * @file InjPipeFlowTableFunction.cpp
 */

#include "InjPipeFlowTableFunction.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/MultivariableNonuniformTableFunctionKernels.hpp"
#include "common/DataTypes.hpp"
#include <algorithm>
#include <fstream>
#include <vector>

namespace geos
{

using namespace dataRepository;

InjPipeFlowTableFunction::InjPipeFlowTableFunction( const string & name,
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


  registerWrapper( viewKeyStruct::wellHeadPressureArray(), &m_whp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of well head pressures " );



  registerWrapper( viewKeyStruct::bottomHolePressureArray(), &m_bhp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of bottom hole pressures representing the dependent variable for the table function.\n"
                    "The quantities must be entered such that the rate index increases first, followed by whp.\n"
                    " For example  bhp[1]= f(rate[1],whp[1]) , bhp[2]= f(rate[2],whp[1])" );

}

void InjPipeFlowTableFunction::postInputInitialization()
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
  initializeFunction();
}

void InjPipeFlowTableFunction::initializeFunction()
{
  localIndex constexpr nDims = 2;
  localIndex constexpr nOps = 1;

  FunctionManager * functionManager = &FunctionManager::getInstance();
  // Create nonuniform version of uniformly spaced table
  m_tableFunction = &(dynamicCast< MultivariableNonuniformTableFunction & >( *functionManager->createChild( MultivariableNonuniformTableFunction::catalogName(), getTableName()  ) ));
  m_tableFunction0 = &dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_b"+getTableName() ) );

  // find max number if independent vars
  int maxVar = std::max( m_rate.size(), m_whp.size());

  // Copy inputs into formats needed by table function

  integer_array axisPoints( nDims );

  axisPoints[0] = m_rate.size();
  axisPoints[1] = m_whp.size();

  integer_array axisPointsr( nDims );
  for( int i=0; i<nDims; i++ )
  {
    axisPointsr[nDims-i-1]= axisPoints[i];
  }

  array1d< array1d< real64 > >axisCoordinates;
  axisCoordinates.resize( nDims );
  array2d< real64 > axisCoord ( nDims, maxVar );
  axisCoordinates[0].resize( m_rate.size());
  for( int i=0; i< m_rate.size(); i++ )
  {
    axisCoordinates[0][i] = m_rate[i];
    axisCoord( 1, i )  = m_rate[i];
  }
  axisCoordinates[1].resize( m_whp.size());
  for( int i=0; i<m_whp.size(); i++ )
  {
    axisCoordinates[1][i] = m_whp[i];
    axisCoord( 0, i )  = m_whp[i];
  }

  setTableCoordinates( nDims, nOps, axisCoord, axisPointsr );
  setTableValues( m_bhp );
  MultivariableNonuniformTableFunction::initializeFunction();
  GEOS_LOG_RANK( GEOS_FMT( "{}: flow table min and max BHP  {} Pa   {} Pa", getName(), *std::min_element( m_bhp.begin(), m_bhp.end()), *std::max_element( m_bhp.begin(), m_bhp.end()) ) );

  m_tableFunction0->setTableCoordinates( axisCoordinates, { units::Dimensionless } );
  m_tableFunction0->setTableValues( m_bhp, units::Dimensionless );
  m_tableFunction0->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  //m_tableFunction0.setInputVarNames( inputVarNames );
  m_tableFunction0->reInitializeFunction();
}
integer
InjPipeFlowTableFunction::getRateBracket( real64 const & rate, integer & b0, integer & b1 ) const
{
  integer nRate = m_rate.size();
  integer stat=0;
  if( rate < m_rate[0] )
  {
    b0=0;
    b1=1;
    stat = 0;
  }
  else if( rate > m_rate[nRate-1] )
  {
    b0=std::max( 0, nRate-2 );
    b1=nRate-1;
    stat= 2;
  }
  else
  {
    integer ifnd=0;
    for( integer j=1; j<nRate-1; j++ )
    {
      if( m_rate[j] <= rate && rate <= m_rate[j+1] )
      {
        ifnd=j;
        break;
      }
    }
    b0=ifnd;
    b1=ifnd+1;
    stat=1;
  }
  return stat;
}
real64
InjPipeFlowTableFunction::calculatedPdQ( real64 const & volRate, real64 const & whp, real64 & ql0, real64 & ql1, real64 & bhp0, real64 & bhp1 ) const
{
  // what is the phase index
  // Calculate ratios that are assumed fixed for table lookup
  real64 lrate = -volRate;

  integer b0, b1;
  integer bStat = getRateBracket( lrate, b0, b1 );
  GEOS_UNUSED_VAR( bStat );
  // Get bhp at lower bracket rates
  ql0 = -m_rate[b0];
  integer solveStat=1;
  calculateBHP( ql0, whp, bhp0, solveStat );

  // Get bhp at upper bracket rates
  ql1 = -m_rate[b1];
  solveStat=1;
  calculateBHP( ql1, whp, bhp1, solveStat );

  real64 dP_dQ = ( bhp1 - bhp0 ) / ( m_rate[b1] - m_rate[b0] );
  return dP_dQ;
}
void InjPipeFlowTableFunction::calculateBHP( real64 const & volRate, real64 const & whp, real64 & bhp, integer & solveStat ) const
{

  MultivariableNonuniformTableFunctionStaticKernel< 2, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );

  TableFunction::KernelWrapper kernelWrapper = m_tableFunction0->createKernelWrapper();

  real64 const m_sign=1.0;
  real64 liq = volRate;

  array1d< real64 > table_coords( 2 );

  array2d< real64 > table_derv( 1, 2 );
  table_coords[0]=liq*m_sign; // gas oil ratio
  table_coords[1]=whp;  // well head pressure

  array1d< real64 > table_coordsr( 2 );

  for( integer i=0; i<2; i++ )
  {
    table_coordsr[2-i-1]=table_coords[i];
  }
  if( solveStat == 0 )
  {
    array1d< real64 > table_bhp( 1 );
    kernel.compute( table_coordsr, table_bhp, table_derv );
    real64 bhpbt = table_bhp[0];
    bhp=bhpbt;
  }
  else
  {
    array1d< real64 > table_bhp( 1 );
    real64 derivatives[2]{};
    table_bhp[0] = kernelWrapper.compute( table_coords, derivatives );
    bhp = table_bhp[0];
  }
  if( whp >m_whp[m_whp.size()-1] )
    solveStat=2;
  else if( whp < m_whp[0] )
    solveStat=0;
  else
    solveStat=1;

}

void InjPipeFlowTableFunction::calculateWHP( const std::string & wellName, real64 const & bhp, real64 const & totalVolumeRate, real64 & whp, integer & solveStat ) const
{

  MultivariableNonuniformTableFunctionStaticKernel< 2, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );
  real64 const m_sign=1.0;
  real64 liq = totalVolumeRate;

  array1d< real64 > table_coords( 2 );
  array1d< real64 > table_bhp( 1 );
  array2d< real64 > table_derv( 1, 2 );
  table_coords[1]=liq*m_sign; // gas oil rat

  table_coords[0]=whp;    // well head pressure
  kernel.compute( table_coords, table_bhp, table_derv );

  integer nWHP = m_whp.size();
  for( integer i=0; i<nWHP; i++ )
  {
    table_coords[0]=m_whp[i];    // well head pressure
    kernel.compute( table_coords, table_bhp, table_derv );
  }
  integer foundBracket= 0;
  table_coords[0]=m_whp[nWHP-1];  // well head pressure
  kernel.compute( table_coords, table_bhp, table_derv );
  double bhpN = table_bhp[0];
  //double bhpN = kernelWrapper.compute( table_coords, derivatives );
  double bhp0, whp0;
  if( bhpN < bhp )
  {
    solveStat=2;
    table_coords[0]=m_whp[nWHP-2];
    kernel.compute( table_coords, table_bhp, table_derv );
    bhp0 = table_bhp[0];
    double dwhp_dp = ( m_whp[nWHP-1] - m_whp[nWHP-2] )/( bhpN - bhp0 );
    whp0 = m_whp[nWHP-1] + ( bhp - bhpN )*dwhp_dp;
    whp=whp0;
    solveStat=2;
  }
  else
  {
    // check low end
    table_coords[0]=m_whp[0];  // well head pressure
    kernel.compute( table_coords, table_bhp, table_derv );
    bhp0 = table_bhp[0];
    //bhp0 = kernelWrapper.compute( table_coords, derivatives );
    whp0  = m_whp[0];
    double whpN;
    if( bhp <  bhp0 )
    {
      solveStat=0;
      table_coords[0]=m_whp[1];
      kernel.compute( table_coords, table_bhp, table_derv );
      bhpN = table_bhp[0];
      whpN = m_whp[0] + ( bhp - bhp0 )*( m_whp[1] - whp0 )/( bhpN - bhp0 );
      whp=whpN;
      solveStat=0;
      return;
    }
    // search for bracketing whp
    for( integer i=1; i<nWHP; ++i )
    {
      table_coords[0]=m_whp[i]; // well head pressure
      kernel.compute( table_coords, table_bhp, table_derv );
      bhpN = table_bhp[0];
      //bhpN = kernelWrapper.compute( table_coords, derivatives );
      if( bhp0 <= bhp && bhp <= bhpN )
      {
        foundBracket=1;
        whpN = m_whp[i];
        break;
      }
      bhp0 = bhpN;
      whp0 = m_whp[i];
    }
    GEOS_THROW_IF( foundBracket ==0,
                   wellName << ": InjPipeFlowTableFunction::calculateWHP failed to find bracketing whp values",
                   InputError, getDataContext());
    if( foundBracket ==0 )
    {
      throw; std::runtime_error( "InjPipeFlowTableFunction::calculateWHP failed to find bracketing whp values" );
    }
    else
    {
      solveStat=1;
      whp  = whp0 + ( bhp - bhp0  )*( whpN - whp0  )/( bhpN - bhp0 );
      solveStat=1;
    }
  }
}

REGISTER_CATALOG_ENTRY( FunctionBase, InjPipeFlowTableFunction, string const &, Group * const )

} // end of namespace geos
