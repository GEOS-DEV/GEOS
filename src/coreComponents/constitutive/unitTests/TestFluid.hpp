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
 * @file TestFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace testing
{

struct Fluid
{
  static constexpr integer NC =  11;

  static constexpr integer H2O =  0;
  static constexpr integer CO2 =  1;
  static constexpr integer N2  =  2;
  static constexpr integer H2S =  3;
  static constexpr integer C1  =  4;
  static constexpr integer C2  =  5;
  static constexpr integer C3  =  6;
  static constexpr integer C4  =  7;
  static constexpr integer C5  =  8;
  static constexpr integer C8  =  9;
  static constexpr integer C10 = 10;

  static constexpr integer Pc = 0;
  static constexpr integer Tc = 1;
  static constexpr integer Ac = 2;
  static constexpr integer Mw = 3;
  static constexpr integer Vs = 4;

  static std::array< real64, 55 > data;
};

template< int NC >
using Feed = std::array< real64, NC >;

template< int NC >
class TestFluid
{
public:
  ~TestFluid() = default;

  static std::unique_ptr< TestFluid< NC > > create( std::array< integer, NC > const & components )
  {
    std::unique_ptr< TestFluid< NC > > testFluid( new TestFluid() );
    createArray( testFluid->criticalPressure, components, Fluid::Pc, Fluid::data );
    createArray( testFluid->criticalTemperature, components, Fluid::Tc, Fluid::data );
    createArray( testFluid->acentricFactor, components, Fluid::Ac, Fluid::data );
    createArray( testFluid->molecularWeight, components, Fluid::Mw, Fluid::data );
    createArray( testFluid->volumeShift, components, Fluid::Vs, Fluid::data );
    return testFluid;
  }

  arrayView1d< real64 const > const getCriticalPressure() const { return criticalPressure.toViewConst(); }
  arrayView1d< real64 const > const getCriticalTemperature() const { return criticalTemperature.toViewConst(); }
  arrayView1d< real64 const > const getCriticalVolume() const { return criticalVolume.toViewConst(); }
  arrayView1d< real64 const > const getAcentricFactor() const { return acentricFactor.toViewConst(); }
  arrayView1d< real64 const > const getMolecularWeight() const { return molecularWeight.toViewConst(); }
  arrayView1d< real64 const > const getVolumeShift() const { return volumeShift.toViewConst(); }

private:
  TestFluid() = default;

  array1d< real64 > criticalPressure;
  array1d< real64 > criticalTemperature;
  array1d< real64 > criticalVolume;
  array1d< real64 > acentricFactor;
  array1d< real64 > molecularWeight;
  array1d< real64 > volumeShift;

private:
  template< typename ARRAY, typename LIST, typename DATAARRAY >
  static void createArray( ARRAY & array, LIST const & indices, integer const row, DATAARRAY const & data )
  {
    for( auto const i : indices )
    {
      array.emplace_back( data[Fluid::NC *row + i] );
    }
  }
public:
  template< typename ARRAY, typename LIST >
  static void createArray( ARRAY & array, LIST const & data )
  {
    for( auto const value : data )
    {
      array.emplace_back( value );
    }
  }
};

std::array< real64, 55 > Fluid::data = {
  // -- Pc
  2.2050e+07, 7.3750e+06, 3.4000e+06, 8.9630e+06, 1.2960e+06, 4.8721e+06,
  4.2481e+06, 3.6400e+06, 4.5990e+06, 2.5300e+06, 1.4600e+06,
  // -- Tc
  6.4700e+02, 3.0410e+02, 1.2620e+02, 3.7353e+02, 3.3150e+01, 3.0532e+02,
  3.6983e+02, 4.0785e+02, 1.9060e+02, 6.2200e+02, 7.8200e+02,
  // -- Ac
  3.4400e-01, 2.3900e-01, 4.0000e-02, 9.4200e-02, -2.1900e-01, 9.9500e-02,
  1.5230e-01, 1.8440e-01, 1.1400e-02, 4.4300e-01, 8.1600e-01,
  // -- Mw
  1.8015e+01, 4.4010e+01, 2.8013e+01, 3.4100e+01, 1.6043e+01, 3.0070e+01,
  4.4097e+01, 5.8124e+01, 7.2151e+01, 1.1423e+02, 1.4228e+02,
  // -- Vs
  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
};

}// testing

}// geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
