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
#include <fstream>
#include <vector>

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
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Type of water fraction  entered in the wfr array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::waterFractionArray(), &m_wfr ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Array of water fractions " );

  registerWrapper( viewKeyStruct::gasFractionType(), &m_gasFractionType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Type of gas fraction  entered in the wfr array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::gasFractionArray(), &m_gfr ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Array of gas fractions " );

  registerWrapper( viewKeyStruct::wellHeadPressureArray(), &m_whp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of well head pressures " );

  registerWrapper( viewKeyStruct::gasLiftArray(), &m_gasLift ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Array of well gas lift rates " );

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
  std::cout << " rate " << m_rate.size() << " whp " << m_whp.size() << " wfr " << m_wfr.size() << " gfr " << m_gfr.size() <<std::endl;
  checkNotDecreasing( m_rate, getRateType());
  checkNotDecreasing( m_whp, "wellHeadPressure" );
  checkNotDecreasing( m_wfr, getWaterFractionType());
  checkNotDecreasing( m_gfr, getGasFractionType());
  checkNotDecreasing( m_gasLift, getGasFractionType());
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
  localIndex constexpr nDims = 5;
  localIndex constexpr nOps = 1;

  FunctionManager * functionManager = &FunctionManager::getInstance();
  // Create nonuniform version of uniformly spaced table
  m_tableFunction = &(dynamicCast< MultivariableNonuniformTableFunction & >( *functionManager->createChild( MultivariableNonuniformTableFunction::catalogName(), getTableName()  ) ));
  m_tableFunction0 = &dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_b" ) );

  // find max number if independent vars
  int maxVar = std::max( m_rate.size(), m_whp.size());
  maxVar = std::max( maxVar, m_wfr.size());
  maxVar = std::max( maxVar, m_gfr.size());
  maxVar = std::max( maxVar, m_gasLift.size());
  // Copy inputs into formats needed by table function

  integer_array axisPoints( nDims );

  axisPoints[0] = m_rate.size();
  axisPoints[1] = m_whp.size();
  axisPoints[2] = m_wfr.size();
  axisPoints[3] = m_gfr.size();
  axisPoints[4] = m_gasLift.size();

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
    axisCoord( 4, i )  = m_rate[i];
  }
  axisCoordinates[1].resize( m_whp.size());
  for( int i=0; i<m_whp.size(); i++ )
  {
    axisCoordinates[1][i] = m_whp[i];
    axisCoord( 3, i )  = m_whp[i];
  }
  axisCoordinates[2].resize( m_wfr.size());
  for( int i=0; i<m_wfr.size(); i++ )
  {
    axisCoordinates[2][i] = m_wfr[i];
    axisCoord( 2, i )  = m_wfr[i];
  }
  axisCoordinates[3].resize( m_gfr.size());
  for( int i=0; i<m_gfr.size(); i++ )
  {
    axisCoordinates[3][i] = m_gfr[i];
    axisCoord( 1, i )  = m_gfr[i];
  }
  axisCoordinates[4].resize( m_gasLift.size());
  for( int i=0; i< m_gasLift.size(); i++ )
  {
    axisCoordinates[4][i] = m_gasLift[i];
    axisCoord( 0, i )  = m_gasLift[i];
  }



  setTableCoordinates( nDims, nOps, axisCoord, axisPointsr );
  setTableValues( m_bhp );
  MultivariableNonuniformTableFunction::initializeFunction();
  std::cout << " bhp table min max " << *std::min_element( m_bhp.begin(), m_bhp.end()) << " " << *std::max_element( m_bhp.begin(), m_bhp.end()) << std::endl;

  m_tableFunction0->setTableCoordinates( axisCoordinates, { units::Dimensionless } );
  m_tableFunction0->setTableValues( m_bhp, units::Dimensionless );
  m_tableFunction0->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  //m_tableFunction0.setInputVarNames( inputVarNames );
  m_tableFunction0->reInitializeFunction();
  writeTable();
}
integer
PipeFlowTableFunction::getRateBracket( real64 const & rate, integer & b0, integer & b1 ) const
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
PipeFlowTableFunction::calculatedPdQ( array1d< real64 > const & phaseRates, real64 const & whp, real64 & ql0, real64 & ql1, real64 & bhp0, real64 & bhp1 ) const
{
  // Calculate ratios that are assumed fixed for table lookup
  real64 lrate = -( phaseRates[0] + phaseRates[2] );
  real64 wfr = -phaseRates[2]/lrate;
  real64 gor = 0;

  if( !isZero( phaseRates[0] ))
    gor =  phaseRates[1]/( phaseRates[0] );

  integer b0, b1;
  integer bStat = getRateBracket( lrate, b0, b1 );
  GEOS_UNUSED_VAR( bStat );
  // Get bhp at lower bracket rates
  ql0 = -m_rate[b0];
  phaseRates[0] = ql0*(1-wfr);
  phaseRates[1] = phaseRates[0]*gor;
  phaseRates[2] = ql0* wfr;
  integer solveStat=1;
  calculateBHP( phaseRates, whp, bhp0, solveStat );

  // Get bhp at upper bracket rates
  ql1 = -m_rate[b1];
  phaseRates[0] = ql1*(1-wfr);
  phaseRates[1] = phaseRates[0]*gor;
  phaseRates[2] = ql1* wfr;
  solveStat=1;
  calculateBHP( phaseRates, whp, bhp1, solveStat );

  real64 dP_dQ = ( bhp1 - bhp0 ) / ( m_rate[b1] - m_rate[b0] );
  return dP_dQ;
}
void PipeFlowTableFunction::calculateBHP( array1d< real64 > const & phaseRates, real64 const & whp, real64 & bhp, integer & solveStat ) const
{

  MultivariableNonuniformTableFunctionStaticKernel< 5, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );

  TableFunction::KernelWrapper kernelWrapper = m_tableFunction0->createKernelWrapper();
  // Assume success
  // liq(oil)=0 vap = 1 wat = 2
  real64 const gasLift=0.0;
  real64 const m_sign=-1.0;
  real64 liq = (phaseRates[0] + phaseRates[2]);
  //for( int i = 0; i < phaseRates.size(); ++i )
  //{
  //  totalVolumeRate += phaseRates[i];
//  }

  std::cout << bhp << " " << phaseRates << " " << whp << std::endl;

  real64 wct = 0;
  if( phaseRates[0]*m_sign > 0 )
  {
    wct = phaseRates[2]/(phaseRates[0]+phaseRates[2]);
  }
#if 0
  if( wct < m_wfr[0] )
  {
    wct = m_wfr[0] + 0.00000001;
  }
  else if( wct > m_wfr[ m_wfr.size() -1 ] )
  {
    wct = m_wfr[ m_wfr.size() -1 ]- +0.00000001;
  }
#endif
  real64 gor = 0;
  if( phaseRates[1]*m_sign  > 0 )
    gor =  phaseRates[1]/(phaseRates[0] );
#if 0
  if( gor < m_gfr[0] )
  {
    gor = m_gfr[0]+   0.00000001;
  }
  else if( gor > m_gfr[ m_gfr.size() -1 ] )
  {
    gor = m_gfr[ m_gfr.size() -1 ]-  0.00000001;
  }
#endif
  //gor = m_gfr[0];
  array1d< real64 > table_coords( 5 );

  array2d< real64 > table_derv( 1, 5 );
  table_coords[0]=liq*m_sign; // gas oil ratio
  table_coords[1]=whp;  // well head pressure
  table_coords[2]=wct;  // water cut
  table_coords[3]=gor; // well head pressure
  table_coords[4]=gasLift; // gas lift rate
  array1d< real64 > table_coordsr( 5 );

  for( integer i=0; i<5; i++ )
  {
    table_coordsr[5-i-1]=table_coords[i];
  }
  if( solveStat == 0 )
  {
    std::cout << "PipeFlowTableFunction::calculateBHP input coords = " << table_coordsr << std::endl;

    array1d< real64 > table_bhp( 1 );
    kernel.compute( table_coordsr, table_bhp, table_derv );
    real64 bhpbt = table_bhp[0];
    // std::cout << " PipeFlowTableFunction::calculateBHP initial bhp = " << bhp << " bhp calc " << table_bhp[0] << " " << table_coordsr
    // <<std::endl;
    std::cout << "PipeFlowTableFunction::calculateBHP 0  output bhp = " << bhp <<    " liq = " << liq << " whp " << whp << " wct = " << wct << " gor = " << gor << " " << table_derv[0][3]<< std::endl;
    bhp=bhpbt;
  }
  else
  {
    array1d< real64 > table_bhp( 1 );
    real64 derivatives[5]{};
    table_bhp[0] = kernelWrapper.compute( table_coords, derivatives );
    bhp = table_bhp[0];
    std::cout << "PipeFlowTableFunction::calculateBHP 1 output bhp = " << bhp <<  " " <<  " liq = " << liq << " whp " << whp << " wct = " << wct << " gor = " << gor << " " << derivatives[1] <<
      std::endl;
  }
  if( whp >m_whp[m_whp.size()-1] )
    solveStat=2;
  else if( whp < m_whp[0] )
    solveStat=0;
  else
    solveStat=1;

}

void PipeFlowTableFunction::calculateWHP( const std::string & wellName, real64 const & bhp, array1d< real64 > const & phaseRates, real64 & whp, integer & solveStat ) const
{

  MultivariableNonuniformTableFunctionStaticKernel< 5, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );
  TableFunction::KernelWrapper kernelWrapper = m_tableFunction0->createKernelWrapper();
  //solveStat = 0;  // Assume success
  // liq(oil)=0 vap = 1 wat = 2
  //real64 totalVolumeRate = 0.0;
  //for( int i = 0; i < phaseRates.size(); ++i )
  //{
  //  totalVolumeRate += phaseRates[i];
//  }
  real64 const gasLift=0.0;
  std::cout << bhp << " " << phaseRates << " " << whp  << std::endl;
  real64 const m_sign=-1.0;
  real64 liq = (phaseRates[0] + phaseRates[2]);
  //for( int i = 0; i < phaseRates.size(); ++i )
  //{
  //  totalVolumeRate += phaseRates[i];
//  }
  std::cout << bhp << " " << phaseRates << " " << whp  << std::endl;

  real64 wct = 0;
  if( phaseRates[0]*m_sign > 0 )
  {
    wct = phaseRates[2]/(phaseRates[0]+phaseRates[2]);
  }
#if 0
  if( wct < m_wfr[0] )
  {
    wct = m_wfr[0] + 0.00000001;
  }
  else if( wct > m_wfr[ m_wfr.size() -1 ] )
  {
    wct = m_wfr[ m_wfr.size() -1 ]- +0.00000001;
  }
#endif
  real64 gor = 0;
  if( phaseRates[1]*m_sign  > 0 )
    gor =  phaseRates[1]/(phaseRates[0] );
#if 0
  if( gor <  m_gfr[0] )
  {
    gor = m_gfr[0]+   0.00000001;
  }
  else if( gor > m_gfr[ m_gfr.size() -1 ] )
  {
    gor = m_gfr[ m_gfr.size() -1 ]-  0.00000001;
  }
#endif
  std::cout << wellName << " PipeFlowTableFunction::calculateWHP  bhp " << bhp << " liq " << liq << " wct " << wct << " gor " << gor << std::endl;
#if 1
  array1d< real64 > table_coords( 5 );
  real64 derivatives[5]{};
  array1d< real64 > table_bhp( 1 );
  array2d< real64 > table_derv( 1, 5 );
  table_coords[4]=liq*m_sign; // gas oil ratio
  table_coords[2]=wct;  // water cut
  table_coords[1]=gor; // well head pressure
  table_coords[0]=gasLift; // gas lift rat

  table_coords[3]=whp;    // well head pressure
  kernel.compute( table_coords, table_bhp, table_derv );
  std::cout << " PipeFlowTableFunction::calculateWHP initial whp = " << whp << " bhp calc " << table_bhp[0] << " " << table_coords <<std::endl;

  integer nWHP = m_whp.size();
  for( integer i=0; i<nWHP; i++ )
  {
    table_coords[3]=m_whp[i];    // well head pressure
    kernel.compute( table_coords, table_bhp, table_derv );
    std::cout << table_coords[3] << " " << table_bhp[0] << std::endl;
  }
  integer foundBracket= 0;
  table_coords[3]=m_whp[nWHP-1];  // well head pressure
  kernel.compute( table_coords, table_bhp, table_derv );
  double bhpN = table_bhp[0];
  //double bhpN = kernelWrapper.compute( table_coords, derivatives );
  double bhp0, whp0;
  if( bhpN < bhp )
  {
    solveStat=2;
    table_coords[3]=m_whp[nWHP-2];
    kernel.compute( table_coords, table_bhp, table_derv );
    bhp0 = table_bhp[0];
    double dwhp_dp = ( m_whp[nWHP-1] - m_whp[nWHP-2] )/( bhpN - bhp0 );
    //bhp0 = kernelWrapper.compute( table_coords, derivatives );
    whp0 = m_whp[nWHP-1] + ( bhp - bhpN )*dwhp_dp;
    //whp0 = m_whp[nWHP-1] + ( bhp - bhpN )*table_derv[0][3];
    std::cout << table_derv << " " << dwhp_dp << " " <<  m_whp[nWHP-1] + ( bhp - bhpN )*table_derv[0][3]<<std::endl;
    whp=whp0;

    if( std::isnan( whp0 ) )
    {
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP bhpN " << bhpN << " bhp0 " << bhp0 << " bhp " << bhp << std::endl;
    }
    std::cout << wellName << " PipeFlowTableFunction::calculateWHP extrapolate at high whp = " << whp << " " << bhp<< " " <<m_whp[nWHP-1] << " " << m_whp[nWHP-2] << " " << bhpN <<" " << bhp0 <<
      " liq = " << liq*m_sign  << std::endl;
    solveStat=2;
  }
  else
  {
    // check low end
    table_coords[3]=m_whp[0];  // well head pressure
    kernel.compute( table_coords, table_bhp, table_derv );
    bhp0 = table_bhp[0];
    //bhp0 = kernelWrapper.compute( table_coords, derivatives );
    whp0  = m_whp[0];
    double whpN;
    if( bhp <  bhp0 )
    {
      solveStat=0;
      table_coords[3]=m_whp[1];
      kernel.compute( table_coords, table_bhp, table_derv );
      bhpN = table_bhp[0];
      //bhpN = kernelWrapper.compute( table_coords, derivatives );
      whpN = m_whp[0] + ( bhp - bhp0 )*( m_whp[1] - whp0 )/( bhpN - bhp0 );
      //whpN = m_whp[0] + ( bhp - bhp0 )*table_derv[0][3];
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP extrapolate at low whp = " << whpN << " liq = " << liq*m_sign << " bhp "<< bhp << " bhpN " << bhpN << " wct = " << wct <<
        " gor = " << gor << std::endl;
      whp=whpN;
      solveStat=0;
      return;
    }
    // search for bracketing whp
    for( integer i=1; i<nWHP; ++i )
    {
      table_coords[3]=m_whp[i]; // well head pressure
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

    if( foundBracket ==0 )
    {
      throw; std::runtime_error( "PipeFlowTableFunction::calculateWHP failed to find bracketing whp values" );
    }
    else
    {
      solveStat=1;
      whp  = whp0 + ( bhp - bhp0 )*( whpN - whp0  )/( bhpN - bhp0 );
      whp  = whp0 + ( bhp - bhp0  )*( whpN - whp0  )/( bhpN - bhp0 );
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP found bracketing " << whpN << " " << whp0 << " " << bhpN << " " << bhp0 << std::endl;
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP found bracketing whp = " << whp0 << " liq = " << liq*m_sign << " bhp "<< bhp << " bhpN " << bhpN << " wct = " << wct <<
        " gor = " << gor << std::endl;
      solveStat=1;
    }
  }
#else
  real64 table_coords[5]{};
  real64 derivatives[5]{};
  table_coords[0]=liq*m_sign; // gas oil ratio
  table_coords[2]=wct;  // water cut
  table_coords[3]=gor; // well head pressure
  table_coords[4]=gasLift; // gas lift rat

  integer nWHP = m_whp.size();
  integer foundBracket= 0;
  table_coords[1]=m_whp[nWHP-1];  // well head pressure
  double bhpN = kernelWrapper.compute( table_coords, derivatives );
  double bhp0, whp0;
  if( bhpN < bhp )
  {
    table_coords[1]=m_whp[nWHP-2];
    bhp0 = kernelWrapper.compute( table_coords, derivatives );
    whp0 = m_whp[nWHP-1] + ( bhp - bhpN )*( m_whp[nWHP-1] - m_whp[nWHP-2] )/( bhpN - bhp0 );
    whp=whp0;

    if( std::isnan( whp0 ) )
    {
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP bhpN " << bhpN << " bhp0 " << bhp0 << " bhp " << bhp << std::endl;
    }
    std::cout << wellName << " PipeFlowTableFunction::calculateWHP extrapolate at high whp = " << whp << " " << bhp<< " " <<m_whp[nWHP-1] << " " << m_whp[nWHP-2] << " " << bhpN <<" " << bhp0 <<
      " liq = " << liq*m_sign  << std::endl;

  }
  else
  {
    // check low end
    table_coords[1]=m_whp[0];  // well head pressure
    bhp0 = kernelWrapper.compute( table_coords, derivatives );
    whp0  = m_whp[0];
    double whpN;
    if( bhp <  bhp0 )
    {
      table_coords[1]=m_whp[1];
      bhpN = kernelWrapper.compute( table_coords, derivatives );
      whpN = m_whp[0] + ( bhp - bhp0 )*( m_whp[1] - whp0 )/( bhpN - bhp0 );
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP extrapolate at low whp = " << whpN << " liq = " << liq*m_sign << " bhp "<< bhp << " bhpN " << bhpN << " wct = " << wct <<
        " gor = " << gor << std::endl;
      whp=whpN;
      return;
    }
    // search for bracketing whp
    for( integer i=1; i<nWHP; ++i )
    {
      table_coords[1]=m_whp[i]; // well head pressure
      bhpN = kernelWrapper.compute( table_coords, derivatives );
      if( bhp0 <= bhp && bhp <= bhpN )
      {
        foundBracket=1;
        whpN = m_whp[i];
        break;
      }
      bhp0 = bhpN;
      whp0 = m_whp[i];

    }

    if( foundBracket ==0 )
    {
      throw; std::runtime_error( "PipeFlowTableFunction::calculateWHP failed to find bracketing whp values" );
    }
    else
    {
      whp  = whp0 + ( bhp - bhp0 )*( whpN - whp0  )/( bhpN - bhp0 );
      whp  = whp0 + ( bhp - bhp0  )*( whpN - whp0  )/( bhpN - bhp0 );
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP found bracketing " << whpN << " " << whp0 << " " << bhpN << " " << bhp0 << std::endl;
      std::cout << wellName << " PipeFlowTableFunction::calculateWHP found bracketing whp = " << whp0 << " liq = " << liq*m_sign << " bhp "<< bhp << " bhpN " << bhpN << " wct = " << wct <<
        " gor = " << gor << std::endl;
    }
  }
#endif
  return;

  std::cout << "PipeFlowTableFunction::calculateWHP input bhp = " << bhp << " liq = " << liq*m_sign << " whp " << whp0 << " wct = " << wct << " gor = " << gor << std::endl;

  // array1d< real64 > table_bhp( 1 );

  //kernel.compute( table_coords, table_bhp, derivs );
  bhp0 = kernelWrapper.compute( table_coords, derivatives );
  double whpn = whp0 + 2e5;
  table_coords[1]= whpn;

  double bhpn =  kernelWrapper.compute( table_coords, derivatives );
  std::cout << " residual bhp = " <<  bhp0 - bhp << std::endl;
  double dpdwhp = ( bhpn - bhp0 )/2e5;
  integer const maxIters=20;
  real64 const tol = 1e-6;
  integer iter = 0;

  while( iter < maxIters && std::abs( bhpn - bhp ) > tol )
  {
    // update whp
    dpdwhp = ( bhpn - bhp0 ) / (whpn - whp0);
    table_coords[1] -= dpdwhp;
    bhp0=bhpn;
    whp0 =whpn;
    whpn = table_coords[1];
    bhpn = kernelWrapper.compute( table_coords, derivatives );
    std::cout << " PipeFlowTableFunction::calculateWHP iter = " << iter << " bhp " << bhpn << " whp = " << whpn <<  " derive " << dpdwhp<< " residual = " << bhpn - bhp0 << std::endl;

    ++iter;
  }
  whp0 = table_coords[1];
  std::cout << "PipeFlowTableFunction::calculateWHP output whp = " << whp0 << " liq = " << liq*m_sign << " bhp " << bhpn << " wct = " << wct << " gor = " << gor << std::endl;
}

void PipeFlowTableFunction::writeTable() const
{

  std::ofstream of;
  std::string filename12 = getTableName() + ".csv";
  of.open( filename12 );



  MultivariableNonuniformTableFunctionStaticKernel< 5, 1 > kernel( getAxisCoordinates(),
                                                                   getAxisPoints(),
                                                                   getAxisSteps(),
                                                                   getAxisStepInvs(),
                                                                   getAxisHypercubeMults(),
                                                                   getHypercubeData()
                                                                   );

  TableFunction::KernelWrapper kernelWrapper = m_tableFunction0->createKernelWrapper();

  of << "rate,wellHeadPressure,waterFraction,gasFraction,bottomHolePressure" << std::endl;
  //array1d< real64 > table_coords( 5 );
  real64 table_coords[5]{};
  real64 derivatives[5]{};
  table_coords[4] = 0.0; //gas lift
  for( integer i=0; i < m_gfr.size(); ++i )
  {
    table_coords[3]=m_gfr[i];   // well head pressure
    // gas oil ratio
    for( real64 j=m_wfr[0]; j < m_wfr[m_wfr.size()-1]; j=j+0.01 )
    {
      table_coords[2]=j;   // water cut
      for( real64 k=m_whp[0]; k < m_whp[m_whp.size()-1]; k=k+10e5 )
      {
        table_coords[1]=k;
        for( integer l=0; l < m_rate.size(); ++l )
        {
          table_coords[0]=m_rate[l]; // liquid rate

          // array2d< real64 > derivs( 1, 5 );
          //kernel.compute( table_coords, table_bhp, derivs );
          real64 table_bhp = kernelWrapper.compute( table_coords, derivatives );
          of << table_coords[0] << "," << table_coords[1] << "," << table_coords[2] << "," << table_coords[3] << "," << table_bhp << std::endl;
        }
      }
    }
  }
  of.close();


}

REGISTER_CATALOG_ENTRY( FunctionBase, PipeFlowTableFunction, string const &, Group * const )

} // end of namespace geos
