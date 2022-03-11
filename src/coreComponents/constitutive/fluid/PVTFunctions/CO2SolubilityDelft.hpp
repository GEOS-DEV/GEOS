/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CO2SolubilityDelft.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITYDELFT_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITYDELFT_HPP_
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "FlashModelBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "functions/TableFunction.hpp"
#include "delftflash/flash/flash.h"

#include <array>
#include <vector>
#include <type_traits>
#include <numeric>
#include <cmath>
namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

constexpr real64 minForDivision1 = 1e-10;

class CO2SolubilityDelftUpdate final : public FlashModelBaseUpdate
{
public:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;

  CO2SolubilityDelftUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                       TableFunction const & CO2SolubilityDelftTable,
                       integer const CO2Index,
                       integer const waterIndex,
                       integer const phaseGasIndex,
                       integer const phaseLiquidIndex )
    : FlashModelBaseUpdate( componentMolarWeight ),
    m_CO2SolubilityDelftTable( CO2SolubilityDelftTable.createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex ),
    m_phaseGasIndex( phaseGasIndex ),
    m_phaseLiquidIndex( phaseLiquidIndex )
  {}

  template< int USD1, int USD2, int USD3 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice1d< real64, USD2 > const & phaseFraction,
                arraySlice2d< real64, USD3 > const & phaseCompFraction ) const;

  template< int USD1 >
  GEOSX_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    FlashModelBaseUpdate::move( space, touch );
    m_CO2SolubilityDelftTable.move( space, touch );
  }

protected:
  template< typename T >
  void 
  doflash(real64 const & pressure, real64 const & temperature,
                 T const & compFraction,
                 std::vector<double> & phaseFraction,
                 std::vector<std::vector<double> > & phaseCompFraction, int & phase_state) const;

  /// Table with CO2 solubility tabulated as a function (P,T)
  TableFunction::KernelWrapper m_CO2SolubilityDelftTable;

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

  /// Index of the gas phase
  integer m_phaseGasIndex;

  /// Index of the liquid phase
  integer m_phaseLiquidIndex;

};

class CO2SolubilityDelft : public FlashModelBase
{
public:

  CO2SolubilityDelft( string const & name,
                 string_array const & inputParams,
                 string_array const & phaseNames,
                 string_array const & componentNames,
                 array1d< real64 > const & componentMolarWeight );

  static string catalogName() { return "CO2SolubilityDelft"; }

  virtual string getCatalogName() const final { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CO2SolubilityDelftUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  /// Table to compute solubility as a function of pressure and temperature
  TableFunction const * m_CO2SolubilityDelftTable;

  /// Index of the CO2 component
  integer m_CO2Index;

  /// Index of the water component
  integer m_waterIndex;

  /// Index of the gas phase
  integer m_phaseGasIndex;

  /// Index of the liquid phase
  integer m_phaseLiquidIndex;
};

template< int USD1, int USD2, int USD3 >
GEOSX_HOST_DEVICE
inline void
CO2SolubilityDelftUpdate::compute( real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD1 > const & compFraction,
                              arraySlice1d< real64, USD2 > const & phaseFraction,
                              arraySlice2d< real64, USD3 > const & phaseCompFraction ) const
{
  GEOSX_ASSERT( 0 );
  // solubility mol/kg(water)  X = Csat/W
  real64 const input[2] = { pressure, temperature };
  real64 solubility = m_CO2SolubilityDelftTable.compute( input );
  solubility *= m_componentMolarWeight[m_waterIndex];

  // Y = C/W = z/(1-z)
  real64 Y;

  if( compFraction[m_CO2Index] > 1.0 - minForDivision1 )
  {
    Y = compFraction[m_CO2Index] / minForDivision1;
  }
  else
  {
    real64 const oneMinusCompFracInv = 1.0 / (1.0 - compFraction[m_CO2Index]);
    Y = compFraction[m_CO2Index] * oneMinusCompFracInv;
  }

  if( Y < solubility )
  {
    // liquid phase only

    // 1) Compute phase fractions

    phaseFraction[m_phaseLiquidIndex] = 1.0;
    phaseFraction[m_phaseGasIndex] = 0.0;

    // 2) Compute phase component fractions

    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction[m_phaseLiquidIndex][ic] = compFraction[ic];
      // the two following lines are not present in Yue's code, unclear if this will have some consequences
      phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
      phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
    }
  }
  else
  {
    // two-phase flow

    // 1) Compute phase fractions

    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    real64 const onePlusYInv = 1.0 / ( 1.0 + Y );
    phaseFraction[m_phaseLiquidIndex] = (solubility + 1.0) * onePlusYInv;
    phaseFraction[m_phaseGasIndex] = 1.0 - phaseFraction[m_phaseLiquidIndex];

    // 2) Compute phase component fractions

    // liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)
    real64 const onePlusSolubilityInv = 1.0 / ( 1.0 + solubility );
    phaseCompFraction[m_phaseLiquidIndex][m_CO2Index] = solubility * onePlusSolubilityInv;

    phaseCompFraction[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction[m_phaseLiquidIndex][m_CO2Index];

    // gas phase composition  CO2 = 1.0
    phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
  }
}

template< typename Array >
typename Array::value_type sum_array( const Array & args )
{
  return std::accumulate( args.begin(), args.end(), typename Array::value_type( 0 ), std::plus< typename Array::value_type >() );
}
template< typename T >
std::vector< T > Normalize( const std::vector< T > & xin )
{
  std::vector< T > xout = xin;
  auto sum = sum_array( xin );
  for( std::size_t n = 0; n != xout.size(); ++n )
  {
    xout[n] = xout[n] / sum;
  }
  return xout;
}

template< int USD1 >
GEOSX_HOST_DEVICE
inline void
CO2SolubilityDelftUpdate::compute( real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD1 > const & compFraction,
                              PhaseProp::SliceType const phaseFraction,
                              PhaseComp::SliceType const phaseCompFraction ) const
{
  using Deriv = multifluid::DerivativeOffset;
  

  // solubility mol/kg(water)  X = Csat/W
  
  real64 const input[2] = { pressure, temperature };
  real64 solubilityDeriv[2]{};
  real64 solubility = m_CO2SolubilityDelftTable.compute( input, solubilityDeriv );

  solubility *= m_componentMolarWeight[m_waterIndex];
  for( integer ic = 0; ic < 2; ++ic )
  {
    solubilityDeriv[ic] *= m_componentMolarWeight[m_waterIndex];
  }

  // Y = C/W = z/(1-z)
  /*
  real64 Y = 0.0;
  real64 dY_dCompFrac[2]{};
  
  if( compFraction[m_CO2Index] > 1.0 - minForDivision1 )
  {
    Y = compFraction[m_CO2Index] / minForDivision1;
    dY_dCompFrac[m_CO2Index] = 1.0 / minForDivision1;
    dY_dCompFrac[m_waterIndex] = 0.0;
  }
  else
  {
    real64 const oneMinusCompFracInv = 1.0 / (1.0 - compFraction[m_CO2Index]);
    Y = compFraction[m_CO2Index] * oneMinusCompFracInv;
    dY_dCompFrac[m_CO2Index] = oneMinusCompFracInv * oneMinusCompFracInv;
    dY_dCompFrac[m_waterIndex] = 0.0;
  }
  */

  auto setZero = []( real64 & val ){ val = 0.0; };
  LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );
/*
  if( Y < solubility )
  {
    // liquid phase only

    // 1) Compute phase fractions

    phaseFraction.value[m_phaseLiquidIndex] = 1.0;
    phaseFraction.value[m_phaseGasIndex] = 0.0;

    // 2) Compute phase component fractions

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction.value[m_phaseLiquidIndex][ic] = compFraction[ic];

      for( localIndex jc = 0; jc < 2; ++jc )
      {
        phaseCompFraction.derivs[m_phaseLiquidIndex][ic][Deriv::dC+jc] = (ic == jc ) ? 1.0 : 0.0;
        phaseCompFraction.derivs[m_phaseGasIndex][ic][Deriv::dC+jc] = 0.0;
      }
    }
    return;
  }
  else
  {
    // two-phase flow

    // 1) Compute phase fractions

    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    real64 const onePlusYInv = 1.0 / ( 1.0 + Y );
    phaseFraction.value[m_phaseLiquidIndex] = (solubility + 1.0) * onePlusYInv;

    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP] = solubilityDeriv[0] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT] = solubilityDeriv[1] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index] =
      -dY_dCompFrac[m_CO2Index] * phaseFraction.value[m_phaseLiquidIndex] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex] =
      -dY_dCompFrac[m_waterIndex] * phaseFraction.value[m_phaseLiquidIndex] * onePlusYInv;

    phaseFraction.value[m_phaseGasIndex] = 1.0 - phaseFraction.value[m_phaseLiquidIndex];

    phaseFraction.derivs[m_phaseGasIndex][Deriv::dP] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dT] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_CO2Index] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_waterIndex] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex];

    // 2) Compute phase component fractions

    // liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)
    real64 const onePlusSolubilityInv = 1.0 / ( 1.0 + solubility );
    phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index] = solubility * onePlusSolubilityInv;

    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP] = solubilityDeriv[0] * (onePlusSolubilityInv*onePlusSolubilityInv);
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT] = solubilityDeriv[1] * (onePlusSolubilityInv*onePlusSolubilityInv);

    phaseCompFraction.value[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index];

    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dP] = -phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP];
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dT] = -phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT];

    // gas phase composition  CO2 = 1.0

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index]   = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;

    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP]   = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP]   = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = 0.0;
    // phaseCompFraction does not depend on globalComponentFraction

  }
  */
  int np = 2;
  int nc = 2; 
  int phase_state_delft = -1;
  std::vector<double> phaseFractionDelft(np);
  std::vector<std::vector<double> > phaseCompFractionDelft(np);
  for (int i = 0; i < nc; i++)
    phaseCompFractionDelft[i].assign(nc,0.);
 
  doflash(pressure, temperature, compFraction, phaseFractionDelft, phaseCompFractionDelft, phase_state_delft);

  if (phase_state_delft == 0)
  {
    
    /*if (phaseFraction.value[0] < 1-1e-2)
    {
      std::cout<<"PURE GAS \n";
      std::cout<<"INPUT P, T, z: "<<" "<<pressure<<" "<<temperature<<" "<<" " <<compFraction[0]<< " "<< compFraction[1]<<"\n"<<std::flush; 
      int i;
      std::cout<<"GEOSX Flash\n";
      std::cout<<phaseFraction.value[0]<<" "<<phaseFraction.value[1]<<"\n";
      std::cout<<"Phase Frac0 "<<phaseCompFraction.value[0][0] << " "<< phaseCompFraction.value[0][1]<<"\n";
      std::cout<<"Phase Frac1 "<<phaseCompFraction.value[1][0] << " "<< phaseCompFraction.value[1][1]<<"\n"<<std::flush;      
      std::cin>>i;
    }*/
    // gas phase only

    // 1) Compute phase fractions

    phaseFraction.value[m_phaseLiquidIndex] = .0;
    phaseFraction.value[m_phaseGasIndex] = 1.0;

    // 2) Compute phase component fractions

    phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index] = .0;
    phaseCompFraction.value[m_phaseLiquidIndex][m_waterIndex] = 1.0;
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction.value[m_phaseGasIndex][ic] = compFraction[ic];

      for( localIndex jc = 0; jc < 2; ++jc )
      {
        phaseCompFraction.derivs[m_phaseGasIndex][ic][Deriv::dC+jc] = (ic == jc ) ? 1.0 : 0.0;
        phaseCompFraction.derivs[m_phaseLiquidIndex][ic][Deriv::dC+jc] = 0.0;
      }
    }        
  }
  else if (phase_state_delft == 1)
  {
   /* if (phaseFraction.value[1] < 1-1e-2)
    {
      std::cout<<"PURE WATER \n";
      std::cout<<"INPUT P, T, z: "<<" "<<pressure<<" "<<temperature<<" "<<" " <<compFraction[0]<< " "<< compFraction[1]<<"\n"<<std::flush; 
      std::cout<<"GEOSX Flash\n";
      std::cout<<phaseFraction.value[0]<<" "<<phaseFraction.value[1]<<"\n";
      std::cout<<"Phase Frac0 "<<phaseCompFraction.value[0][0] << " "<< phaseCompFraction.value[0][1]<<"\n";
      std::cout<<"Phase Frac1 "<<phaseCompFraction.value[1][0] << " "<< phaseCompFraction.value[1][1]<<"\n"<<std::flush;
      int i;
      std::cin>>i;
    }*/
    // liquid phase only

    // 1) Compute phase fractions

    phaseFraction.value[m_phaseLiquidIndex] = 1.0;
    phaseFraction.value[m_phaseGasIndex] = 0.0;

    // 2) Compute phase component fractions

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction.value[m_phaseLiquidIndex][ic] = compFraction[ic];

      for( localIndex jc = 0; jc < 2; ++jc )
      {
        phaseCompFraction.derivs[m_phaseLiquidIndex][ic][Deriv::dC+jc] = (ic == jc ) ? 1.0 : 0.0;
        phaseCompFraction.derivs[m_phaseGasIndex][ic][Deriv::dC+jc] = 0.0;
      }
    }    
  }
  else
  { 
    //std::cout<<"TWO PHASE\n";
    /*double delta = (fabs(phaseFraction.value[0] - phaseFractionDelft[0]) + fabs(phaseFraction.value[1] - phaseFractionDelft[1]))/2*100;
    if (delta > 5)
    {
      std::cout<<std::flush<<"Delft Flash\n";
      std::cout<<"INPUT P, T, z: "<<" "<<pressure<<" "<<temperature<<" "<<" " <<compFraction[0]<< " "<< compFraction[1]<<"\n"<<std::flush; 
      std::cout<<"PhasesState "<<phase_state_delft<<"\n";
	    std::cout << "V: " << phaseFractionDelft[0] << " " << phaseFractionDelft[1] <<  "\n";
      std::cout<<"Phase Frac0 "<<phaseCompFractionDelft[0][0] << " "<< phaseCompFractionDelft[0][1]<<"\n";
      std::cout<<"Phase Frac1 "<<phaseCompFractionDelft[1][0] << " "<< phaseCompFractionDelft[1][1]<<"\n"<<std::flush;
      std::cout<<"GEOSX Flash\n";
      std::cout<<fabs(phaseFraction.value[0])<<" "<<fabs(phaseFraction.value[1])<<"\n";
      std::cout<<"Phase Frac0 "<<phaseCompFraction.value[0][0] << " "<< phaseCompFraction.value[0][1]<<"\n";
      std::cout<<"Phase Frac1 "<<phaseCompFraction.value[1][0] << " "<< phaseCompFraction.value[1][1]<<"\n"<<std::flush;
      int i = 0;
      std::cin>>i;

    }  
    */
    // two-phase flow

    // 1) Compute phase fractions


    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    for (int ip = 0; ip < np; ip++)
    {
      phaseFraction.value[ip] = phaseFractionDelft[ip];
      for (int ic = 0; ic < nc; ic++)
        phaseCompFraction.value[ip][ic] = phaseCompFractionDelft[ip][ic];
    }
    //needed for derivatives
    int phase_state_delft_der = -1;
    std::vector<double> phaseFractionDelft_der(np);
    std::vector<std::vector<double> > phaseCompFractionDelft_der(np);
    for (int i = 0; i < nc; i++)
    phaseCompFractionDelft_der[i].assign(nc,0.); 
     double const sqrtPrecision = sqrt( std::numeric_limits< double >::epsilon() );   
    //dPress
    {
      const real64 dPressure = sqrtPrecision * ( std::fabs( pressure ) + sqrtPrecision );
      doflash(pressure + dPressure, temperature, compFraction, phaseFractionDelft_der, phaseCompFractionDelft_der, phase_state_delft_der);
      for (int ip = 0; ip < np; ip++)
      {
        phaseFraction.derivs[ip][Deriv::dP] = (phaseFractionDelft_der[ip] - phaseFractionDelft[ip])/dPressure;
        for (int ic = 0; ic < nc; ic++)
        {
          phaseCompFraction.derivs[ip][ic][Deriv::dP] = (phaseCompFractionDelft_der[ip][ic] - phaseCompFractionDelft[ip][ic])/dPressure;        
        }
      }

    }
    //dT
    {
      const real64 dTemp = sqrtPrecision * ( std::fabs( temperature ) + sqrtPrecision );
      doflash(pressure, temperature + dTemp, compFraction, phaseFractionDelft_der, phaseCompFractionDelft_der, phase_state_delft_der);
      for (int ip = 0; ip < np; ip++)
      {
        phaseFraction.derivs[ip][Deriv::dT] = (phaseFractionDelft_der[ip] - phaseFractionDelft[ip])/dTemp;
        for (int ic = 0; ic < nc; ic++)
        {
          phaseCompFraction.derivs[ip][ic][Deriv::dT] = (phaseCompFractionDelft_der[ip][ic] - phaseCompFractionDelft[ip][ic])/dTemp;
        }
      }

    }    
    // Feed
    {
      for( int icder = 0; icder < nc; icder++)
      {
        real64 dz = sqrtPrecision * ( std::fabs( compFraction[icder] ) + sqrtPrecision );
        if( compFraction[icder] + dz > 1 )
        {
          dz = -dz;
        }
        std::vector< real64 > newFeed(nc, 0.);
        for (int ic = 0; ic < nc; ic++)
          newFeed[ic] = compFraction[ic];
        newFeed[icder] += dz;
        newFeed = Normalize( newFeed );
        doflash(pressure, temperature, newFeed, phaseFractionDelft_der, phaseCompFractionDelft_der, phase_state_delft_der);
        for (int ip = 0; ip < np; ip++)
        {
          phaseFraction.derivs[ip][Deriv::dC + icder] = (phaseFractionDelft_der[ip] - phaseFractionDelft[ip])/dz;
          for (int ic2 = 0; ic2 < nc; ic2++)
          {
            phaseCompFraction.derivs[ip][ic2][Deriv::dC + icder] = (phaseCompFractionDelft_der[ip][ic2] - phaseCompFractionDelft[ip][ic2])/dz;
          }
        }
      
      }
    }
    /*
    // phaseCompFraction does not depend on globalComponentFraction
 // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)
    real64 const onePlusYInv = 1.0 / ( 1.0 + Y );
    phaseFraction.value[m_phaseLiquidIndex] = (solubility + 1.0) * onePlusYInv;

    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP] = solubilityDeriv[0] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT] = solubilityDeriv[1] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index] =
      -dY_dCompFrac[m_CO2Index] * phaseFraction.value[m_phaseLiquidIndex] * onePlusYInv;
    phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex] =
      -dY_dCompFrac[m_waterIndex] * phaseFraction.value[m_phaseLiquidIndex] * onePlusYInv;

    phaseFraction.value[m_phaseGasIndex] = 1.0 - phaseFraction.value[m_phaseLiquidIndex];

    phaseFraction.derivs[m_phaseGasIndex][Deriv::dP] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dP];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dT] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dT];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_CO2Index] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_CO2Index];
    phaseFraction.derivs[m_phaseGasIndex][Deriv::dC+m_waterIndex] = -phaseFraction.derivs[m_phaseLiquidIndex][Deriv::dC+m_waterIndex];

    // 2) Compute phase component fractions

    // liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)
    real64 const onePlusSolubilityInv = 1.0 / ( 1.0 + solubility );
    phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index] = solubility * onePlusSolubilityInv;

    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP] = solubilityDeriv[0] * (onePlusSolubilityInv*onePlusSolubilityInv);
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT] = solubilityDeriv[1] * (onePlusSolubilityInv*onePlusSolubilityInv);

    phaseCompFraction.value[m_phaseLiquidIndex][m_waterIndex] = 1.0 - phaseCompFraction.value[m_phaseLiquidIndex][m_CO2Index];

    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dP] = -phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dP];
    phaseCompFraction.derivs[m_phaseLiquidIndex][m_waterIndex][Deriv::dT] = -phaseCompFraction.derivs[m_phaseLiquidIndex][m_CO2Index][Deriv::dT];

    // gas phase composition  CO2 = 1.0

    phaseCompFraction.value[m_phaseGasIndex][m_CO2Index]   = 1.0;
    phaseCompFraction.value[m_phaseGasIndex][m_waterIndex] = 0.0;

    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP]   = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_CO2Index][Deriv::dP]   = 0.0;
    phaseCompFraction.derivs[m_phaseGasIndex][m_waterIndex][Deriv::dT] = 0.0;
    // phaseCompFraction does not depend on globalComponentFraction
    */
   /*
    std::cout<<std::flush<<"Flash\n";
    std::cout<<"INPUT P, T, z: "<<" "<<pressure<<" "<<temperature<<" "<<" " <<compFraction[0]<< " "<< compFraction[1]<<"\n"<<std::flush; 
    std::cout<<"PhasesState "<<phase_state_delft<<"\n";
	  std::cout << "V: " << phaseFraction.value[0] << " " << phaseFraction.value[1] <<  "\n";
    std::cout<<"Phase Frac0 "<<phaseCompFraction.value[0][0] << " "<< phaseCompFraction.value[0][1]<<"\n";
    std::cout<<"Phase Frac1 "<<phaseCompFraction.value[1][0] << " "<< phaseCompFraction.value[1][1]<<"\n"<<std::flush; 
    
    std::cout<<"Dpres "<<phaseFraction.derivs[0][Deriv::dP] <<" "<<phaseFraction.derivs[1][Deriv::dP]<<"\n"
               << phaseCompFraction.derivs[0][0][Deriv::dP] << " "<<phaseCompFraction.derivs[0][1][Deriv::dP]<<" "
               << phaseCompFraction.derivs[1][0][Deriv::dP] << " "<<phaseCompFraction.derivs[1][1][Deriv::dP]<<"\n";
    std::cout<<"DT "<<phaseFraction.derivs[0][Deriv::dT] <<" "<<phaseFraction.derivs[1][Deriv::dP]<<"\n"
               << phaseCompFraction.derivs[0][0][Deriv::dT] << " "<<phaseCompFraction.derivs[0][1][Deriv::dP]<<" "
               << phaseCompFraction.derivs[1][0][Deriv::dT] << " "<<phaseCompFraction.derivs[1][1][Deriv::dP]<<"\n";
    std::cout<<"DZ "<<phaseFraction.derivs[0][Deriv::dC] <<" "<<phaseFraction.derivs[0][Deriv::dC+1]<<" "<<
                        phaseFraction.derivs[1][Deriv::dC]<<" "<<phaseFraction.derivs[1][Deriv::dC+1]<<"\n"
               << phaseCompFraction.derivs[0][0][Deriv::dC] << " "<<" "<<phaseCompFraction.derivs[0][0][Deriv::dC+1]<<" "<<
                  phaseCompFraction.derivs[0][1][Deriv::dC]<<" "<<phaseCompFraction.derivs[0][1][Deriv::dC+1]<<" "
               << phaseCompFraction.derivs[1][0][Deriv::dC] << " "<<phaseCompFraction.derivs[1][0][Deriv::dC+1] 
               <<" "<<phaseCompFraction.derivs[1][1][Deriv::dC]<<" "<<phaseCompFraction.derivs[1][1][Deriv::dC+1]<<"\n";  
  */
  
  
  }
  
  
}
#if 1
template< typename T >
void 
CO2SolubilityDelftUpdate::doflash(real64 const & pressure, real64 const & temperature, T const & compFraction,
                                  std::vector<real64> & phaseFraction,
                                  std::vector<std::vector<real64> > & phaseCompFraction, int &phase_state) const
{
  std::vector<std::string> comp{"CO2", "H2O"};
  // std::vector<std::string> comp{"CO2", "H2O"};
  std::vector<std::string> ph{"Aq", "V"};

 
  std::vector<double> z{compFraction[m_CO2Index], compFraction[m_waterIndex]};
  //normalize
  double eps = 1e-8;
  double sum = 0;
  if (z[0] < eps)
    z[0]=eps;
  if (z[1] < eps)
    z[1] = eps;
  sum = z[0]+z[1];
  z[0] /= sum;
  z[1] /= sum;
  std::unordered_map<std::string, std::string> eos = {{"V", "PR"}, {"L", "PR"}, {"Aq", "AQ1"}};
	// define which EoS for each phase
	// Aq: 0) activity model (eos_aq1.cpp) 1) CSMGem activity model (eos_aq2.cpp)
	// V: 0) PR-EoS (eos_pr.cpp) 2) SRK-EoS (eos_srk.cpp)
  NPhaseSSI flash(comp, eos, ph);
  std::vector<std::string> phases;
 
  //convert to bars and Kelvins
  double pbar = pressure*1e-5;
  double tempK = temperature + 273.15; 

  //run flash
  phases = flash.runFlash(pbar, tempK, z); 
	
  //get results
  std::vector<double> &V = flash.getV();
	std::vector<double> &x = flash.getx();  
  if (phases.size() == 1 && phases[0] == "V" )
  {
    phase_state = 0;
    phaseFraction[m_phaseLiquidIndex] = 0.;
    phaseFraction[m_phaseGasIndex] = 1.0;
    for( localIndex ic = 0; ic < 2; ++ic )
      phaseCompFraction[m_phaseGasIndex][ic] = compFraction[ic];
    phaseCompFraction[m_phaseLiquidIndex][m_CO2Index] = .0;
    phaseCompFraction[m_phaseLiquidIndex][m_waterIndex] = 1.0;
/*
    if (phaseFraction.value[0] < 1-1e-2)
    {
      std::cout<<"PURE GAS \n";
      std::cout<<"INPUT P, T, z: "<<" "<<pbar<<" "<<tempK<<" "<<" " <<z[0]<< " "<< z[1]<<"\n"<<std::flush; 
      int i;
      std::cout<<"GEOSX Flash\n";
      std::cout<<phaseFraction.value[0]<<" "<<phaseFraction.value[1]<<"\n";
      std::cout<<"Phase Frac0 "<<phaseCompFraction.value[0][0] << " "<< phaseCompFraction.value[0][1]<<"\n";
      std::cout<<"Phase Frac1 "<<phaseCompFraction.value[1][0] << " "<< phaseCompFraction.value[1][1]<<"\n"<<std::flush;      
      std::cin>>i;
    }
    */
  }
  else if (phases.size() == 1 && phases[0] == "Aq" )
  {
    phase_state = 1;
    phaseFraction[m_phaseLiquidIndex] = 1.0;
    phaseFraction[m_phaseGasIndex] = 0.0;
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction[m_phaseLiquidIndex][ic] = compFraction[ic];
    }
    phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction[m_phaseGasIndex][m_waterIndex] = .0;
    
    /*
    if (phaseFraction.value[1] < 1-1e-2)
    {
      std::cout<<"PURE WATER \n";
      std::cout<<"INPUT P, T, z: "<<" "<<pbar<<" "<<tempK<<" "<<" " <<z[0]<< " "<< z[1]<<"\n"<<std::flush;     
      std::cout<<"GEOSX Flash\n";
      std::cout<<phaseFraction.value[0]<<" "<<phaseFraction.value[1]<<"\n";
      std::cout<<"Phase Frac0 "<<phaseCompFraction.value[0][0] << " "<< phaseCompFraction.value[0][1]<<"\n";
      std::cout<<"Phase Frac1 "<<phaseCompFraction.value[1][0] << " "<< phaseCompFraction.value[1][1]<<"\n"<<std::flush;
      int i;
      std::cin>>i;
    }*/
  }
  else
  { 
    phase_state = 2;
    int idx_vapor_delft = 1;
    int idx_aq_delft = 0;    
    phaseFraction[m_phaseLiquidIndex] = V[idx_aq_delft];
    phaseFraction[m_phaseGasIndex] = V[idx_vapor_delft];
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      phaseCompFraction[m_phaseGasIndex][ic] = x[2*idx_vapor_delft + ic];
      phaseCompFraction[m_phaseLiquidIndex][ic] = x[2*idx_aq_delft + ic];
    }
    /*
    {
    std::cout<<"TWO PHASEINSIDE\n";
    std::cout<<"INPUT P, T, z: "<<" "<<pbar<<" "<<tempK<<" "<<" " <<z[0]<< " "<< z[1]<<"\n"<<std::flush; 
    std::cout<<"Phases "<<V.size() << " "<<phases.size()<<" "<<phases[0]<<"\n";
	  std::cout << "V: " << V[0] << " " << V[1] <<  "\n";
	  std::cout << "x0: " << x[0] << " " << x[1] <<  "\n";
	  std::cout << "x1: " << x[2] << " " << x[3] << "\n";
    int i;
    std::cin>>i;
    } 
    */ 
  }  
}
#endif

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITYDELFT_HPP_
