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
 * @file TestFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"

#include <unordered_map>

namespace geos
{

namespace testing
{

struct Fluid
{
  /* UNCRUSTIFY-OFF */
  static constexpr integer H2O    =  0; // water
  static constexpr integer CO2    =  1; // carbon dioxide
  static constexpr integer N2     =  2; // nitrogen
  static constexpr integer H2S    =  3; // hydrogen sulfide
  static constexpr integer H2     =  4; // hydrogen
  static constexpr integer CH4    =  5; // methane
  static constexpr integer C2H6   =  6; // ethane
  static constexpr integer C3H8   =  7; // propane
  static constexpr integer C4H10  =  8; // butane
  static constexpr integer C5H12  =  9; // pentane
  static constexpr integer C8H18  = 10; // octane
  static constexpr integer C10H22 = 11; // decane
  static constexpr integer NACL   = 12; // sodium chloride
  static constexpr integer KCL    = 13; // potassium chloride

  static constexpr integer NC    = 14; // number of components

  static constexpr integer Mw =  0; // Molecular weight
  static constexpr integer Pc =  1; // Critical Pressure
  static constexpr integer Tc =  2; // Critical Temperature
  static constexpr integer Vc =  3; // Critical Volume
  static constexpr integer Ac =  4; // Acentric Factor
  static constexpr integer Pr =  5; // Parachor

  static constexpr integer NP = 6; // Number of properties

  static constexpr std::array<real64,84> data = {
  //  Mw            Pc            Tc            Vc            Ac            Pr              Name
      1.80153e-02,  2.20640e+07,  6.47096e+02,  5.59480e-05,  3.44300e-01,  9.36563e-06, // H2O    (water)
      4.40095e-02,  7.37730e+06,  3.04128e+02,  9.41185e-05,  2.23940e-01,  7.37268e-06, // CO2    (carbon dioxide)
      2.80134e-02,  3.39580e+06,  1.26192e+02,  8.94142e-05,  3.72000e-02,  0.00000e+00, // N2     (nitrogen)
      3.40809e-02,  9.00000e+06,  3.73100e+02,  9.81354e-05,  1.00500e-01,  1.43639e-05, // H2S    (hydrogen sulfide)
      2.01588e-03,  1.29640e+06,  3.31450e+01,  6.44828e-05, -2.19000e-01,  0.00000e+00, // H2     (hydrogen)
      1.60425e-02,  4.59920e+06,  1.90564e+02,  9.86278e-05,  1.14200e-02,  0.00000e+00, // CH4    (methane)
      3.00690e-02,  4.87220e+06,  3.05322e+02,  1.45839e-04,  9.95000e-02,  1.12855e-05, // C2H6   (ethane)
      4.40956e-02,  4.25120e+06,  3.69890e+02,  2.00000e-04,  1.52100e-01,  2.60222e-05, // C3H8   (propane)
      5.81222e-02,  3.79600e+06,  4.25125e+02,  2.54922e-04,  2.01000e-01,  3.36663e-05, // C4H10  (butane)
      7.21488e-02,  3.36750e+06,  4.69700e+02,  3.11526e-04,  2.51000e-01,  4.11451e-05, // C5H12  (pentane)
      1.14229e-01,  2.48359e+06,  5.68740e+02,  4.92368e-04,  3.98000e-01,  6.27709e-05, // C8H18  (octane)
      1.42282e-01,  2.10300e+06,  6.17700e+02,  6.09756e-04,  4.88400e-01,  7.72027e-05, // C10H22 (decane)
      5.84428e-02,  3.55000e+07,  3.40000e+03,  2.66000e-04,  1.89400e-01,  2.89938e-05, // NACL   (sodium chloride)
      7.45513e-02,  1.80000e+07,  3.47000e+03,  6.25000e-04, -1.21200e-01,  3.05992e-05, // KCL    (potassium chloride)
  };

  static stdUnorderedMap<integer,string> const componentNames;
  /* UNCRUSTIFY-ON */
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
    for( integer const ic : components )
    {
      testFluid->componentNames.emplace_back( Fluid::componentNames.at( ic ) );
    }
    createArray( testFluid->criticalPressure, components, Fluid::Pc, Fluid::data );
    createArray( testFluid->criticalTemperature, components, Fluid::Tc, Fluid::data );
    createArray( testFluid->criticalVolume, components, Fluid::Vc, Fluid::data );
    createArray( testFluid->acentricFactor, components, Fluid::Ac, Fluid::data );
    createArray( testFluid->molecularWeight, components, Fluid::Mw, Fluid::data );
    testFluid->volumeShift.resize( NC );
    testFluid->volumeShift.zero();
    testFluid->binaryCoeff.resize( NC, NC );
    testFluid->binaryCoeff.zero();
    return testFluid;
  }

  template< typename LIST >
  void setBinaryCoefficients( LIST const & bics )
  {
    auto bic = bics.begin();
    for( int i = 0; i < NC; ++i )
    {
      for( int j = 0; j < i; ++j )
      {
        binaryCoeff( i, j ) = *bic++;
        binaryCoeff( j, i ) = binaryCoeff( i, j );
      }
      binaryCoeff( i, i ) = 0.0;
    }
  }

  constitutive::compositional::ComponentProperties & getComponentProperties()
  {
    if( !m_component_properties )
    {
      m_component_properties = std::make_unique< constitutive::compositional::ComponentProperties >(
        componentNames,
        molecularWeight );
      createArray( m_component_properties->m_componentCriticalPressure, criticalPressure );
      createArray( m_component_properties->m_componentCriticalTemperature, criticalTemperature );
      createArray( m_component_properties->m_componentAcentricFactor, acentricFactor );
      createArray( m_component_properties->m_componentVolumeShift, volumeShift );
      m_component_properties->m_componentBinaryCoeff.resize( NC, NC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        for( integer jc = 0; jc < NC; ++jc )
        {
          m_component_properties->m_componentBinaryCoeff( ic, jc ) = binaryCoeff( ic, jc );
        }
      }
    }
    return *m_component_properties;
  }

  constitutive::compositional::ComponentProperties::KernelWrapper createKernelWrapper()
  {
    return getComponentProperties().createKernelWrapper();
  }

private:
  TestFluid() = default;

public:
  string_array componentNames;
  array1d< real64 > criticalPressure;
  array1d< real64 > criticalTemperature;
  array1d< real64 > criticalVolume;
  array1d< real64 > acentricFactor;
  array1d< real64 > molecularWeight;
  array1d< real64 > volumeShift;
  array2d< real64 > binaryCoeff;

private:
  std::unique_ptr< constitutive::compositional::ComponentProperties > m_component_properties{};

private:
  template< typename ARRAY, typename LIST, typename DATAARRAY >
  static void createArray( ARRAY & array, LIST const & indices, integer const row, DATAARRAY const & data )
  {
    for( auto const i : indices )
    {
      array.emplace_back( data[i*Fluid::NP + row] );
    }
  }
public:
  template< typename ARRAY, typename LIST >
  static void createArray( ARRAY & array, LIST const & data )
  {
    for( auto const & value : data )
    {
      array.emplace_back( value );
    }
  }
  template< typename ARRAY, typename LIST >
  static void populateArray( ARRAY & array, LIST const & data )
  {
    integer i = 0;
    for( auto const value : data )
    {
      array[i++] = value;
    }
  }
};

/* UNCRUSTIFY-OFF */
stdUnorderedMap<integer,string> const Fluid::componentNames = {
  { H2O,    "H2O"    }, // water
  { CO2,    "CO2"    }, // carbon dioxide
  { N2,     "N2"     }, // nitrogen
  { H2S,    "H2S"    }, // hydrogen sulfide
  { H2,     "H2"     }, // hydrogen
  { CH4,    "CH4"    }, // methane
  { C2H6,   "C2H6"   }, // ethane
  { C3H8,   "C3H8"   }, // propane
  { C4H10,  "C4H10"  }, // butane
  { C5H12,  "C5H12"  }, // pentane
  { C8H18,  "C8H18"  }, // octane
  { C10H22, "C10H22" }, // decane
  { NACL,   "NACL"   }, // sodium chloride
  { KCL,    "KCL"    }, // potassium chloride
};
/* UNCRUSTIFY-ON */

}// testing
}// geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUID_HPP_
