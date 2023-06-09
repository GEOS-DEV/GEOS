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

struct FluidCatalogue
{
  static std::map< std::string, real64 const > const Pc;
  static std::map< std::string, real64 const > const Tc;
  static std::map< std::string, real64 const > const Ac;
  static std::map< std::string, real64 const > const Mw;
};

template< int NC >
using Feed = std::array< real64, NC >;

template< int NC >
class TestFluid
{
public:
  ~TestFluid() = default;

  static std::unique_ptr< TestFluid< NC > > create( std::array< std::string, NC > const & componentNames )
  {
    std::unique_ptr< TestFluid< NC > > testFluid( new TestFluid());
    Feed< NC > data;
    std::transform( componentNames.begin(), componentNames.end(), data.begin(), []( auto const & name ){ return 1.0e5 * FluidCatalogue::Pc.at( name ); } );
    createArray( testFluid->criticalPressure, data );
    std::transform( componentNames.begin(), componentNames.end(), data.begin(), []( auto const & name ){ return FluidCatalogue::Tc.at( name ); } );
    createArray( testFluid->criticalTemperature, data );
    std::transform( componentNames.begin(), componentNames.end(), data.begin(), []( auto const & name ){ return FluidCatalogue::Ac.at( name ); } );
    createArray( testFluid->acentricFactor, data );
    std::transform( componentNames.begin(), componentNames.end(), data.begin(), []( auto const & name ){ return FluidCatalogue::Mw.at( name ); } );
    createArray( testFluid->molecularWeight, data );
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

std::map< std::string, real64 const > const FluidCatalogue::Pc = {
  {"H2O", 2.2050e+02}, {"CO2", 7.3750e+01}, {"N2", 3.4000e+01}, {"H2S", 8.9630e+01},
  {"C1", 1.2960e+01}, {"C2", 4.8721e+01}, {"C3", 4.2481e+01}, {"C4", 3.6400e+01},
  {"C5", 4.5990e+01}, {"C8", 2.5300e+01}, {"C10", 1.4600e+01}
};
std::map< std::string, real64 const > const FluidCatalogue::Tc = {
  {"H2O", 6.4700e+02}, {"CO2", 3.0410e+02}, {"N2", 1.2620e+02}, {"H2S", 3.7353e+02},
  {"C1", 3.3150e+01}, {"C2", 3.0532e+02}, {"C3", 3.6983e+02}, {"C4", 4.0785e+02},
  {"C5", 1.9060e+02}, {"C8", 6.2200e+02}, {"C10", 7.8200e+02}
};
std::map< std::string, real64 const > const FluidCatalogue::Ac = {
  {"H2O", 3.4400e-01}, {"CO2", 2.3900e-01}, {"N2", 4.0000e-02}, {"H2S", 9.4200e-02},
  {"C1", -2.1900e-01}, {"C2", 9.9500e-02}, {"C3", 1.5230e-01}, {"C4", 1.8440e-01},
  {"C5", 1.1400e-02}, {"C8", 4.4300e-01}, {"C10", 8.1600e-01}
};
std::map< std::string, real64 const > const FluidCatalogue::Mw = {
  {"H2O", 1.8015e+01}, {"CO2", 4.4010e+01}, {"N2", 2.8013e+01}, {"H2S", 3.4100e+01},
  {"C1", 1.6043e+01}, {"C2", 3.0070e+01}, {"C3", 4.4097e+01}, {"C4", 5.8124e+01},
  {"C5", 7.2151e+01}, {"C8", 1.1423e+02}, {"C10", 1.4228e+02}
};

}// testing

}// geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
