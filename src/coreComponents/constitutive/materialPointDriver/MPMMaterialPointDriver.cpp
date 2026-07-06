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
 * @file MPMMaterialPointDriver.cpp
 *
 * Standalone material-point driver for exercising MPM constitutive models without
 * constructing an MPM grid.  This executable is built only when
 * GEOS_ENABLE_MPM_MATERIAL_POINT_DRIVER is enabled, so the production MPM solver
 * and its default build are unchanged.
 */

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TypeDispatch.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/ContinuumBase.hpp"
#include "constitutive/solid/SolidUtilities.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "conduit.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <type_traits>
#include <utility>
#include <vector>

namespace geos
{
namespace constitutive
{
namespace materialPointDriver
{

using dataRepository::Group;
using dataRepository::Wrapper;
using dataRepository::WrapperBase;

constexpr localIndex k = 0;
constexpr localIndex q = 0;
constexpr int VOIGT_SIZE = 6;
constexpr int DIM = 3;

using Vec3 = std::array< real64, 3 >;
using Vec6 = std::array< real64, 6 >;
using Mat3 = std::array< std::array< real64, 3 >, 3 >;

enum class ControlMode
{
  StrainIncrement,
  StrainRate,
  TrueStrainRate,
  Stress
};

enum class MaterialFrameUpdate
{
  Auto,
  Fixed,
  Rotation,
  Fiber,
  Normal,
  Graphite,
  MpmCofactor
};

enum class EnergyMode
{
  Off,
  StressPower,
  Material
};

enum class TemperatureMode
{
  Prescribed,
  Isothermal,
  Adiabatic,
  FromMaterial
};

enum class StressControlFailurePolicy
{
  Error,
  Stop,
  Continue
};

enum class StressControlAlgorithm
{
  Newton,
  RegularizedNewton,
  Servo,
  Hybrid
};

enum class StressControlDiagnosticsLevel
{
  Off,
  Step,
  Iteration,
  Full
};

struct InitialField
{
  std::string name;
  std::vector< real64 > values;
};

struct Options
{
  std::string materialXmlPath;
  std::string materialName;
  std::string pathCsvPath;
  std::string outputCsvPath = "material_point.csv";

  std::array< ControlMode, 6 > controlModes = {
    ControlMode::StrainIncrement,
    ControlMode::StrainIncrement,
    ControlMode::StrainIncrement,
    ControlMode::StrainIncrement,
    ControlMode::StrainIncrement,
    ControlMode::StrainIncrement };

  Vec6 initialStress = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  std::vector< InitialField > initialFields;

  std::vector< real64 > materialDirectionValues = { 1.0, 0.0, 0.0,
                                                    0.0, 1.0, 0.0,
                                                    0.0, 0.0, 1.0 };
  Vec3 tangentHint = { 1.0, 0.0, 0.0 };
  MaterialFrameUpdate materialFrameUpdate = MaterialFrameUpdate::Auto;

  real64 initialTemperature = 300.0;
  real64 initialTemperatureRate = 0.0;
  real64 initialSpecificInternalEnergy = 0.0;
  real64 initialDensity = -1.0;
  real64 initialLengthScale = 1.0;
  real64 initialStrengthScale = 1.0;
  real64 heatCapacity = 1.0;
  real64 retentionFactor = 1.0;

  EnergyMode energyMode = EnergyMode::StressPower;
  TemperatureMode temperatureMode = TemperatureMode::Isothermal;

  real64 stressTolerance = 1.0e-8;
  real64 finiteDifferenceEpsilon = 1.0e-8;
  int maxNewtonIterations = 25;
  int maxLineSearchIterations = 12;
  int maxStressBracketIterations = 32;
  int maxStressBisectionIterations = 64;
  real64 stressBracketInitialScale = 0.0;
  real64 stressBracketMaxStrain = 5.0e-2;
  real64 stressBracketGrowth = 2.0;
  StressControlAlgorithm stressControlAlgorithm = StressControlAlgorithm::Hybrid;
  real64 stressControlRegularization = 1.0e-12;
  real64 stressControlMaxStrainCorrection = 5.0e-2;
  real64 stressControlServoCompliance = 1.0e-2;
  real64 stressControlServoRelaxation = 0.5;
  real64 stressControlServoDerivativeFloor = 1.0e-8;
  int stressControlServoIterations = 12;
  int stressControlPatternIterations = 0;
  real64 stressControlPatternInitialStep = 0.0;
  real64 stressControlPatternMinStep = 1.0e-12;
  real64 stressControlPatternShrink = 0.5;
  real64 stressControlPatternGrowth = 1.25;
  Vec6 stressControlMinStrain = { -std::numeric_limits< real64 >::infinity(),
                                  -std::numeric_limits< real64 >::infinity(),
                                  -std::numeric_limits< real64 >::infinity(),
                                  -std::numeric_limits< real64 >::infinity(),
                                  -std::numeric_limits< real64 >::infinity(),
                                  -std::numeric_limits< real64 >::infinity() };
  Vec6 stressControlMaxStrain = {  std::numeric_limits< real64 >::infinity(),
                                   std::numeric_limits< real64 >::infinity(),
                                   std::numeric_limits< real64 >::infinity(),
                                   std::numeric_limits< real64 >::infinity(),
                                   std::numeric_limits< real64 >::infinity(),
                                   std::numeric_limits< real64 >::infinity() };
  StressControlFailurePolicy stressControlFailurePolicy = StressControlFailurePolicy::Continue;
  std::string stressControlDiagnosticsPath;
  StressControlDiagnosticsLevel stressControlDiagnosticsLevel = StressControlDiagnosticsLevel::Off;

  localIndex constantSteps = 0;
  real64 constantDt = 0.0;
  Vec6 constantValues = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

struct PathStep
{
  real64 dt = 0.0;
  Vec6 values = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
};

struct DriverState
{
  real64 time = 0.0;
  Mat3 F = {};
  Mat3 materialFrameReference = {};
  Mat3 materialFrame = {};
  Vec6 totalStrain = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  Vec6 previousStrainIncrement = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  Vec6 stress = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  real64 referenceDensity = 1.0;
  real64 density = 1.0;
  real64 jacobian = 1.0;
  real64 lengthScale = 1.0;
  real64 strengthScale = 1.0;
  real64 temperature = 300.0;
  real64 temperatureRate = 0.0;
  real64 specificInternalEnergy = 0.0;
  real64 accumulatedStressPower = 0.0;
  real64 lastStressPower = 0.0;
};

struct TrialResult
{
  DriverState state;
  Mat3 L = {};
  Vec6 strainIncrement = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  int newtonIterations = 0;
  real64 residualNorm = 0.0;
  bool converged = true;
};

std::string trim( std::string const & input )
{
  auto begin = input.begin();
  while( begin != input.end() && std::isspace( static_cast< unsigned char >( *begin ) ) )
  {
    ++begin;
  }
  auto end = input.end();
  do
  {
    if( begin == end )
    {
      break;
    }
    --end;
  }
  while( std::isspace( static_cast< unsigned char >( *end ) ) );
  return std::string( begin, end + 1 );
}

std::string toLower( std::string input )
{
  std::transform( input.begin(), input.end(), input.begin(), []( unsigned char c )
  {
    return static_cast< char >( std::tolower( c ) );
  } );
  return input;
}

std::vector< std::string > split( std::string const & input, char const delimiter )
{
  std::vector< std::string > result;
  std::stringstream ss( input );
  std::string item;
  while( std::getline( ss, item, delimiter ) )
  {
    result.push_back( trim( item ) );
  }
  return result;
}

std::vector< real64 > parseRealList( std::string const & input )
{
  std::vector< real64 > values;
  for( std::string const & token : split( input, ',' ) )
  {
    if( token.empty() )
    {
      continue;
    }
    values.push_back( std::stod( token ) );
  }
  return values;
}

Vec6 parseVec6( std::string const & input )
{
  std::vector< real64 > values = parseRealList( input );
  if( values.size() != 6 )
  {
    throw std::runtime_error( "Expected six comma-separated values, got " + std::to_string( values.size() ) );
  }
  Vec6 out = {};
  for( int i = 0; i < 6; ++i )
  {
    out[i] = values[static_cast< std::size_t >( i )];
  }
  return out;
}

Vec3 parseVec3( std::string const & input )
{
  std::vector< real64 > values = parseRealList( input );
  if( values.size() != 3 )
  {
    throw std::runtime_error( "Expected three comma-separated values, got " + std::to_string( values.size() ) );
  }
  return Vec3{ values[0], values[1], values[2] };
}

std::string voigtComponentName( int const component )
{
  static std::array< char const *, 6 > const names = { "xx", "yy", "zz", "yz", "xz", "xy" };
  if( component >= 0 && component < 6 )
  {
    return names[static_cast< std::size_t >( component )];
  }
  return "unknown";
}

int parseVoigtComponentIndex( std::string const & input )
{
  std::string const name = toLower( trim( input ) );
  static std::array< char const *, 6 > const names = { "xx", "yy", "zz", "yz", "xz", "xy" };
  for( int c = 0; c < 6; ++c )
  {
    if( name == names[static_cast< std::size_t >( c )] )
    {
      return c;
    }
  }
  throw std::runtime_error( "Unknown Voigt component name in strain bound: " + input );
}

real64 parseStrainBoundValue( std::string const & token, bool const isLowerBound )
{
  std::string const value = toLower( trim( token ) );
  if( value.empty() || value == "none" || value == "off" || value == "nan" || value == "na" )
  {
    return isLowerBound ? -std::numeric_limits< real64 >::infinity()
                        :  std::numeric_limits< real64 >::infinity();
  }
  if( value == "inf" || value == "+inf" || value == "infinity" || value == "+infinity" )
  {
    return std::numeric_limits< real64 >::infinity();
  }
  if( value == "-inf" || value == "-infinity" )
  {
    return -std::numeric_limits< real64 >::infinity();
  }
  return std::stod( value );
}

Vec6 parseStressControlStrainBound( std::string const & input, bool const isLowerBound )
{
  real64 const disabled = isLowerBound ? -std::numeric_limits< real64 >::infinity()
                                       :  std::numeric_limits< real64 >::infinity();
  Vec6 bound = { disabled, disabled, disabled, disabled, disabled, disabled };
  std::vector< std::string > const tokens = split( input, ',' );

  bool const assignments = input.find( '=' ) != std::string::npos;
  if( assignments )
  {
    for( std::string const & token : tokens )
    {
      if( token.empty() )
      {
        continue;
      }
      std::size_t const pos = token.find( '=' );
      if( pos == std::string::npos )
      {
        throw std::runtime_error( "Stress-control strain bounds must be all component assignments or all six values" );
      }
      int const component = parseVoigtComponentIndex( token.substr( 0, pos ) );
      bound[static_cast< std::size_t >( component )] = parseStrainBoundValue( token.substr( pos + 1 ), isLowerBound );
    }
    return bound;
  }

  if( tokens.size() != 6 )
  {
    throw std::runtime_error( "Stress-control strain bound expects six comma-separated values or component assignments" );
  }
  for( int c = 0; c < 6; ++c )
  {
    bound[static_cast< std::size_t >( c )] = parseStrainBoundValue( tokens[static_cast< std::size_t >( c )], isLowerBound );
  }
  return bound;
}

std::string readFile( std::string const & path )
{
  std::ifstream input( path );
  if( !input )
  {
    throw std::runtime_error( "Could not open file: " + path );
  }
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

void writeUsage( std::ostream & os )
{
  os << "Usage: geos_mpm_material_point_driver --material-xml material.xml --material-name name "
     << "--path load_path.csv --output result.csv [options]\n\n"
     << "Load path CSV columns are: dt,xx,yy,zz,yz,xz,xy.  The six component values "
     << "are interpreted according to --control.\n\n"
     << "Key options:\n"
     << "  --control strain,stress,trueStrainRate,stress,strain,strain\n"
     << "  --initial-stress sx,sy,sz,syz,sxz,sxy\n"
     << "  --material-direction n0,n1,n2 or nine row-wise 3x3 values\n"
     << "  --material-tangent-hint x,y,z\n"
     << "  --material-direction-update auto|fixed|rotation|fiber|normal|graphite|mpmCofactor\n"
     << "  --temperature-mode prescribed|isothermal|adiabatic|fromMaterial\n"
     << "  --energy-mode off|stressPower|material\n"
     << "  --stress-control-failure-policy error|stop|continue\n"
     << "  --stress-control-algorithm newton|regularizedNewton|servo|hybrid\n"
     << "  --stress-control-regularization value\n"
     << "  --stress-control-max-strain-correction value\n"
     << "  --stress-control-servo-compliance value\n"
     << "  --stress-control-servo-relaxation value\n"
     << "  --stress-control-servo-derivative-floor value\n"
     << "  --stress-control-servo-iterations n\n"
     << "  --stress-control-pattern-iterations n\n"
     << "  --stress-control-pattern-initial-step value (0 means automatic)\n"
     << "  --stress-control-pattern-min-step value\n"
     << "  --stress-control-pattern-shrink value\n"
     << "  --stress-control-pattern-growth value\n"
     << "  --stress-control-min-strain xx=0,yy=0 or six comma-separated values\n"
     << "  --stress-control-max-strain xx=0.5,yy=0.5 or six comma-separated values\n"
     << "  --stress-control-diagnostics path/to/trace.csv\n"
     << "  --stress-control-diagnostics-level off|step|iteration|full\n"
     << "  --stress-bracket-initial-scale value (0 means automatic)\n"
     << "  --stress-bracket-max-strain value\n"
     << "  --stress-bracket-growth value\n"
     << "  --max-stress-bracket-iterations n\n"
     << "  --max-stress-bisection-iterations n\n"
     << "  --initial-field name=v0[,v1,...] (may be repeated)\n";
}

ControlMode parseControlMode( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "strain" || mode == "strainincrement" || mode == "increment" )
  {
    return ControlMode::StrainIncrement;
  }
  if( mode == "strainrate" || mode == "engineeringstrainrate" )
  {
    return ControlMode::StrainRate;
  }
  if( mode == "truestrainrate" || mode == "logstrainrate" )
  {
    return ControlMode::TrueStrainRate;
  }
  if( mode == "stress" || mode == "stresscontrol" )
  {
    return ControlMode::Stress;
  }
  throw std::runtime_error( "Unknown control mode: " + token );
}

MaterialFrameUpdate parseMaterialFrameUpdate( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "auto" ) return MaterialFrameUpdate::Auto;
  if( mode == "fixed" ) return MaterialFrameUpdate::Fixed;
  if( mode == "rotation" || mode == "rotate" ) return MaterialFrameUpdate::Rotation;
  if( mode == "fiber" || mode == "fibre" ) return MaterialFrameUpdate::Fiber;
  if( mode == "normal" ) return MaterialFrameUpdate::Normal;
  if( mode == "graphite" ) return MaterialFrameUpdate::Graphite;
  if( mode == "mpmcofactor" || mode == "cofactor" ) return MaterialFrameUpdate::MpmCofactor;
  throw std::runtime_error( "Unknown material direction update mode: " + token );
}

EnergyMode parseEnergyMode( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "off" || mode == "none" ) return EnergyMode::Off;
  if( mode == "stresspower" || mode == "stress_power" ) return EnergyMode::StressPower;
  if( mode == "material" || mode == "frommaterial" ) return EnergyMode::Material;
  throw std::runtime_error( "Unknown energy mode: " + token );
}

TemperatureMode parseTemperatureMode( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "prescribed" ) return TemperatureMode::Prescribed;
  if( mode == "isothermal" ) return TemperatureMode::Isothermal;
  if( mode == "adiabatic" ) return TemperatureMode::Adiabatic;
  if( mode == "frommaterial" || mode == "material" ) return TemperatureMode::FromMaterial;
  throw std::runtime_error( "Unknown temperature mode: " + token );
}

StressControlFailurePolicy parseStressControlFailurePolicy( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "error" || mode == "abort" || mode == "fail" ) return StressControlFailurePolicy::Error;
  if( mode == "stop" || mode == "partial" || mode == "partialoutput" ) return StressControlFailurePolicy::Stop;
  if( mode == "continue" || mode == "keepgoing" || mode == "keep_going" ) return StressControlFailurePolicy::Continue;
  throw std::runtime_error( "Unknown stress-control failure policy: " + token );
}

StressControlAlgorithm parseStressControlAlgorithm( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "newton" || mode == "exact" ) return StressControlAlgorithm::Newton;
  if( mode == "regularizednewton" || mode == "regularized" || mode == "lm" ||
      mode == "levenberg" || mode == "levenbergmarquardt" || mode == "levenberg_marquardt" )
  {
    return StressControlAlgorithm::RegularizedNewton;
  }
  if( mode == "servo" || mode == "barostat" || mode == "controller" || mode == "relaxed" ) return StressControlAlgorithm::Servo;
  if( mode == "hybrid" || mode == "robust" || mode == "newtonservo" || mode == "newton_servo" ||
      mode == "newtonlmservo" || mode == "newton_lm_servo" )
  {
    return StressControlAlgorithm::Hybrid;
  }
  throw std::runtime_error( "Unknown stress-control algorithm: " + token );
}

StressControlDiagnosticsLevel parseStressControlDiagnosticsLevel( std::string const & token )
{
  std::string const mode = toLower( token );
  if( mode == "off" || mode == "none" || mode == "false" || mode == "0" ) return StressControlDiagnosticsLevel::Off;
  if( mode == "step" || mode == "summary" ) return StressControlDiagnosticsLevel::Step;
  if( mode == "iteration" || mode == "iterations" || mode == "iter" ) return StressControlDiagnosticsLevel::Iteration;
  if( mode == "full" || mode == "debug" || mode == "trace" || mode == "all" ) return StressControlDiagnosticsLevel::Full;
  throw std::runtime_error( "Unknown stress-control diagnostics level: " + token );
}

bool stressControlUsesRegularization( StressControlAlgorithm const algorithm )
{
  return algorithm == StressControlAlgorithm::RegularizedNewton || algorithm == StressControlAlgorithm::Hybrid;
}

bool stressControlUsesServo( StressControlAlgorithm const algorithm )
{
  return algorithm == StressControlAlgorithm::Servo || algorithm == StressControlAlgorithm::Hybrid;
}

bool permissiveStressControlFailurePolicy( StressControlFailurePolicy const policy )
{
  return policy == StressControlFailurePolicy::Stop || policy == StressControlFailurePolicy::Continue;
}

bool stressControlDiagnosticsAtLeast( StressControlDiagnosticsLevel const level,
                                      StressControlDiagnosticsLevel const requested )
{
  return static_cast< int >( level ) >= static_cast< int >( requested );
}

bool hasStressControl( Options const & options )
{
  for( ControlMode const mode : options.controlModes )
  {
    if( mode == ControlMode::Stress )
    {
      return true;
    }
  }
  return false;
}

Options parseCommandLine( int argc, char ** argv )
{
  Options options;
  for( int i = 1; i < argc; ++i )
  {
    std::string const key = argv[i];
    if( key == "--help" || key == "-h" )
    {
      writeUsage( std::cout );
      std::exit( 0 );
    }

    auto requireValue = [&]( std::string const & optionName ) -> std::string
    {
      if( i + 1 >= argc )
      {
        throw std::runtime_error( "Missing value after " + optionName );
      }
      ++i;
      return std::string( argv[i] );
    };

    if( key == "--material-xml" ) options.materialXmlPath = requireValue( key );
    else if( key == "--material-name" ) options.materialName = requireValue( key );
    else if( key == "--path" ) options.pathCsvPath = requireValue( key );
    else if( key == "--output" ) options.outputCsvPath = requireValue( key );
    else if( key == "--control" )
    {
      std::vector< std::string > tokens = split( requireValue( key ), ',' );
      if( tokens.size() != 6 )
      {
        throw std::runtime_error( "--control requires six comma-separated control modes" );
      }
      for( int c = 0; c < 6; ++c )
      {
        options.controlModes[c] = parseControlMode( tokens[static_cast< std::size_t >( c )] );
      }
    }
    else if( key == "--initial-stress" ) options.initialStress = parseVec6( requireValue( key ) );
    else if( key == "--material-direction" ) options.materialDirectionValues = parseRealList( requireValue( key ) );
    else if( key == "--material-tangent-hint" ) options.tangentHint = parseVec3( requireValue( key ) );
    else if( key == "--material-direction-update" ) options.materialFrameUpdate = parseMaterialFrameUpdate( requireValue( key ) );
    else if( key == "--initial-temperature" ) options.initialTemperature = std::stod( requireValue( key ) );
    else if( key == "--initial-temperature-rate" ) options.initialTemperatureRate = std::stod( requireValue( key ) );
    else if( key == "--initial-specific-internal-energy" ) options.initialSpecificInternalEnergy = std::stod( requireValue( key ) );
    else if( key == "--initial-density" ) options.initialDensity = std::stod( requireValue( key ) );
    else if( key == "--initial-length-scale" ) options.initialLengthScale = std::stod( requireValue( key ) );
    else if( key == "--initial-strength-scale" ) options.initialStrengthScale = std::stod( requireValue( key ) );
    else if( key == "--heat-capacity" ) options.heatCapacity = std::stod( requireValue( key ) );
    else if( key == "--retention-factor" ) options.retentionFactor = std::stod( requireValue( key ) );
    else if( key == "--energy-mode" ) options.energyMode = parseEnergyMode( requireValue( key ) );
    else if( key == "--temperature-mode" ) options.temperatureMode = parseTemperatureMode( requireValue( key ) );
    else if( key == "--stress-tolerance" ) options.stressTolerance = std::stod( requireValue( key ) );
    else if( key == "--fd-epsilon" ) options.finiteDifferenceEpsilon = std::stod( requireValue( key ) );
    else if( key == "--max-newton-iterations" ) options.maxNewtonIterations = std::stoi( requireValue( key ) );
    else if( key == "--max-line-search-iterations" ) options.maxLineSearchIterations = std::stoi( requireValue( key ) );
    else if( key == "--max-stress-bracket-iterations" || key == "--max-stress-control-bracket-iterations" ) options.maxStressBracketIterations = std::stoi( requireValue( key ) );
    else if( key == "--max-stress-bisection-iterations" || key == "--max-stress-control-bisection-iterations" ) options.maxStressBisectionIterations = std::stoi( requireValue( key ) );
    else if( key == "--stress-bracket-initial-scale" || key == "--stress-control-bracket-initial-scale" ) options.stressBracketInitialScale = std::stod( requireValue( key ) );
    else if( key == "--stress-bracket-max-strain" || key == "--stress-control-bracket-max-strain" ) options.stressBracketMaxStrain = std::stod( requireValue( key ) );
    else if( key == "--stress-bracket-growth" || key == "--stress-control-bracket-growth" ) options.stressBracketGrowth = std::stod( requireValue( key ) );
    else if( key == "--stress-control-algorithm" || key == "--stress-algorithm" || key == "--stress-control-method" ) options.stressControlAlgorithm = parseStressControlAlgorithm( requireValue( key ) );
    else if( key == "--stress-control-regularization" || key == "--stress-regularization" || key == "--stress-lm-regularization" ) options.stressControlRegularization = std::stod( requireValue( key ) );
    else if( key == "--stress-control-max-strain-correction" || key == "--stress-max-strain-correction" ) options.stressControlMaxStrainCorrection = std::stod( requireValue( key ) );
    else if( key == "--stress-control-servo-compliance" || key == "--stress-servo-compliance" || key == "--barostat-compliance" ) options.stressControlServoCompliance = std::stod( requireValue( key ) );
    else if( key == "--stress-control-servo-relaxation" || key == "--stress-servo-relaxation" || key == "--barostat-gain" ) options.stressControlServoRelaxation = std::stod( requireValue( key ) );
    else if( key == "--stress-control-servo-derivative-floor" || key == "--stress-servo-derivative-floor" ) options.stressControlServoDerivativeFloor = std::stod( requireValue( key ) );
    else if( key == "--stress-control-servo-iterations" || key == "--stress-servo-iterations" || key == "--barostat-iterations" ) options.stressControlServoIterations = std::stoi( requireValue( key ) );
    else if( key == "--stress-control-pattern-iterations" || key == "--stress-pattern-iterations" || key == "--max-stress-pattern-iterations" ) options.stressControlPatternIterations = std::stoi( requireValue( key ) );
    else if( key == "--stress-control-pattern-initial-step" || key == "--stress-pattern-initial-step" ) options.stressControlPatternInitialStep = std::stod( requireValue( key ) );
    else if( key == "--stress-control-pattern-min-step" || key == "--stress-pattern-min-step" ) options.stressControlPatternMinStep = std::stod( requireValue( key ) );
    else if( key == "--stress-control-pattern-shrink" || key == "--stress-pattern-shrink" ) options.stressControlPatternShrink = std::stod( requireValue( key ) );
    else if( key == "--stress-control-pattern-growth" || key == "--stress-pattern-growth" ) options.stressControlPatternGrowth = std::stod( requireValue( key ) );
    else if( key == "--stress-control-min-strain" || key == "--stress-min-strain" || key == "--min-stress-control-strain" )
    {
      options.stressControlMinStrain = parseStressControlStrainBound( requireValue( key ), true );
    }
    else if( key == "--stress-control-max-strain" || key == "--stress-max-strain" || key == "--max-stress-control-strain" )
    {
      options.stressControlMaxStrain = parseStressControlStrainBound( requireValue( key ), false );
    }
    else if( key == "--stress-control-diagnostics" || key == "--stress-diagnostics" || key == "--stress-control-trace" )
    {
      options.stressControlDiagnosticsPath = requireValue( key );
      if( options.stressControlDiagnosticsLevel == StressControlDiagnosticsLevel::Off )
      {
        options.stressControlDiagnosticsLevel = StressControlDiagnosticsLevel::Iteration;
      }
    }
    else if( key == "--stress-control-diagnostics-level" || key == "--stress-diagnostics-level" || key == "--stress-control-trace-level" )
    {
      options.stressControlDiagnosticsLevel = parseStressControlDiagnosticsLevel( requireValue( key ) );
    }
    else if( key == "--stress-control-failure-policy" || key == "--stress-failure-policy" )
    {
      options.stressControlFailurePolicy = parseStressControlFailurePolicy( requireValue( key ) );
    }
    else if( key == "--allow-partial-output" )
    {
      options.stressControlFailurePolicy = StressControlFailurePolicy::Stop;
    }
    else if( key == "--continue-on-stress-control-failure" || key == "--continue-on-nonconvergence" )
    {
      options.stressControlFailurePolicy = StressControlFailurePolicy::Continue;
    }
    else if( key == "--steps" ) options.constantSteps = static_cast< localIndex >( std::stoll( requireValue( key ) ) );
    else if( key == "--dt" ) options.constantDt = std::stod( requireValue( key ) );
    else if( key == "--values" ) options.constantValues = parseVec6( requireValue( key ) );
    else if( key == "--initial-field" )
    {
      std::string field = requireValue( key );
      std::size_t const pos = field.find( '=' );
      if( pos == std::string::npos )
      {
        throw std::runtime_error( "--initial-field expects name=value[,value...]" );
      }
      options.initialFields.push_back( InitialField{ field.substr( 0, pos ), parseRealList( field.substr( pos + 1 ) ) } );
    }
    else
    {
      throw std::runtime_error( "Unknown command-line option: " + key );
    }
  }

  if( options.materialXmlPath.empty() )
  {
    throw std::runtime_error( "--material-xml is required" );
  }
  if( options.materialName.empty() )
  {
    throw std::runtime_error( "--material-name is required" );
  }
  if( options.pathCsvPath.empty() && ( options.constantSteps <= 0 || options.constantDt <= 0.0 ) )
  {
    throw std::runtime_error( "Provide either --path or both --steps and --dt" );
  }
  if( options.materialDirectionValues.size() != 3 && options.materialDirectionValues.size() != 9 )
  {
    throw std::runtime_error( "--material-direction expects either 3 values or 9 row-wise 3x3 values" );
  }
  if( options.heatCapacity <= 0.0 )
  {
    throw std::runtime_error( "--heat-capacity must be positive" );
  }
  if( options.stressControlRegularization < 0.0 )
  {
    throw std::runtime_error( "--stress-control-regularization must be non-negative" );
  }
  if( options.stressControlMaxStrainCorrection < 0.0 )
  {
    throw std::runtime_error( "--stress-control-max-strain-correction must be non-negative" );
  }
  if( options.stressControlServoCompliance < 0.0 )
  {
    throw std::runtime_error( "--stress-control-servo-compliance must be non-negative" );
  }
  if( options.stressControlServoRelaxation < 0.0 )
  {
    throw std::runtime_error( "--stress-control-servo-relaxation must be non-negative" );
  }
  if( options.stressControlServoDerivativeFloor < 0.0 )
  {
    throw std::runtime_error( "--stress-control-servo-derivative-floor must be non-negative" );
  }
  if( options.stressControlServoIterations < 0 )
  {
    throw std::runtime_error( "--stress-control-servo-iterations must be non-negative" );
  }
  if( options.stressControlPatternIterations < 0 )
  {
    throw std::runtime_error( "--stress-control-pattern-iterations must be non-negative" );
  }
  if( options.stressControlPatternInitialStep < 0.0 )
  {
    throw std::runtime_error( "--stress-control-pattern-initial-step must be non-negative" );
  }
  if( options.stressControlPatternMinStep < 0.0 )
  {
    throw std::runtime_error( "--stress-control-pattern-min-step must be non-negative" );
  }
  if( options.stressControlPatternShrink <= 0.0 || options.stressControlPatternShrink >= 1.0 )
  {
    throw std::runtime_error( "--stress-control-pattern-shrink must be greater than 0 and less than 1" );
  }
  if( options.stressControlPatternGrowth < 1.0 )
  {
    throw std::runtime_error( "--stress-control-pattern-growth must be at least 1" );
  }
  for( int c = 0; c < 6; ++c )
  {
    if( options.stressControlMinStrain[static_cast< std::size_t >( c )] >
        options.stressControlMaxStrain[static_cast< std::size_t >( c )] )
    {
      throw std::runtime_error( "Stress-control strain lower bound exceeds upper bound for component " +
                                voigtComponentName( c ) );
    }
  }
  return options;
}

std::vector< PathStep > readPathCsv( Options const & options )
{
  if( options.pathCsvPath.empty() )
  {
    std::vector< PathStep > path( static_cast< std::size_t >( options.constantSteps ) );
    for( PathStep & step : path )
    {
      step.dt = options.constantDt;
      step.values = options.constantValues;
    }
    return path;
  }

  std::ifstream input( options.pathCsvPath );
  if( !input )
  {
    throw std::runtime_error( "Could not open load path CSV: " + options.pathCsvPath );
  }

  std::string headerLine;
  if( !std::getline( input, headerLine ) )
  {
    throw std::runtime_error( "Load path CSV is empty: " + options.pathCsvPath );
  }
  std::vector< std::string > headers = split( headerLine, ',' );
  std::unordered_map< std::string, std::size_t > column;
  for( std::size_t i = 0; i < headers.size(); ++i )
  {
    column[toLower( headers[i] )] = i;
  }

  std::array< std::string, 6 > const componentNames = { "xx", "yy", "zz", "yz", "xz", "xy" };
  if( column.count( "dt" ) == 0 )
  {
    throw std::runtime_error( "Load path CSV must contain a dt column" );
  }
  for( std::string const & name : componentNames )
  {
    if( column.count( name ) == 0 )
    {
      throw std::runtime_error( "Load path CSV must contain column " + name );
    }
  }

  std::vector< PathStep > path;
  std::string line;
  localIndex lineNumber = 1;
  while( std::getline( input, line ) )
  {
    ++lineNumber;
    if( trim( line ).empty() )
    {
      continue;
    }
    std::vector< std::string > tokens = split( line, ',' );
    if( tokens.size() < headers.size() )
    {
      throw std::runtime_error( "Too few columns in load path CSV at line " + std::to_string( lineNumber ) );
    }
    PathStep step;
    step.dt = std::stod( tokens[column["dt"]] );
    for( int c = 0; c < 6; ++c )
    {
      step.values[c] = std::stod( tokens[column[componentNames[c]]] );
    }
    path.push_back( step );
  }
  return path;
}

Mat3 identityMatrix()
{
  Mat3 A = {};
  for( int i = 0; i < 3; ++i )
  {
    A[i][i] = 1.0;
  }
  return A;
}

real64 dot( Vec3 const & a, Vec3 const & b )
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 cross( Vec3 const & a, Vec3 const & b )
{
  return Vec3{ a[1] * b[2] - a[2] * b[1],
               a[2] * b[0] - a[0] * b[2],
               a[0] * b[1] - a[1] * b[0] };
}

real64 norm( Vec3 const & a )
{
  return std::sqrt( dot( a, a ) );
}

Vec3 normalize( Vec3 v )
{
  real64 const n = norm( v );
  if( n <= 1.0e-30 )
  {
    throw std::runtime_error( "Encountered zero-magnitude vector while constructing material frame" );
  }
  for( real64 & x : v )
  {
    x /= n;
  }
  return v;
}

Vec3 matrixVectorProduct( Mat3 const & A, Vec3 const & x )
{
  Vec3 y = { 0.0, 0.0, 0.0 };
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      y[i] += A[i][j] * x[j];
    }
  }
  return y;
}

Mat3 matrixProduct( Mat3 const & A, Mat3 const & B )
{
  Mat3 C = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      for( int l = 0; l < 3; ++l )
      {
        C[i][j] += A[i][l] * B[l][j];
      }
    }
  }
  return C;
}

Mat3 addScaledIdentity( Mat3 A, real64 const scale )
{
  for( int i = 0; i < 3; ++i )
  {
    A[i][i] += scale;
  }
  return A;
}

Mat3 cofactor( Mat3 const & F )
{
  real64 f[3][3] = {};
  real64 c[3][3] = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      f[i][j] = F[i][j];
    }
  }
  LvArray::tensorOps::cofactor< 3 >( c, f );
  Mat3 out = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      out[i][j] = c[i][j];
    }
  }
  return out;
}

real64 determinant( Mat3 const & F )
{
  real64 f[3][3] = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      f[i][j] = F[i][j];
    }
  }
  return LvArray::tensorOps::determinant< 3 >( f );
}

Mat3 polarRotation( Mat3 const & F )
{
  real64 f[3][3] = {};
  real64 r[3][3] = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      f[i][j] = F[i][j];
    }
  }
  LvArray::tensorOps::polarDecomposition< 3 >( r, f );
  Mat3 out = {};
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      out[i][j] = r[i][j];
    }
  }
  return out;
}

Mat3 buildMaterialFrame( std::vector< real64 > const & values, Vec3 const & tangentHint )
{
  Mat3 frame = {};
  if( values.size() == 9 )
  {
    for( int i = 0; i < 3; ++i )
    {
      Vec3 row{ values[static_cast< std::size_t >( 3 * i )],
                values[static_cast< std::size_t >( 3 * i + 1 )],
                values[static_cast< std::size_t >( 3 * i + 2 )] };
      if( i == 0 )
      {
        row = normalize( row );
      }
      else
      {
        for( int previous = 0; previous < i; ++previous )
        {
          Vec3 prev{ frame[previous][0], frame[previous][1], frame[previous][2] };
          real64 const projection = dot( row, prev );
          for( int c = 0; c < 3; ++c )
          {
            row[c] -= projection * prev[c];
          }
        }
        row = normalize( row );
      }
      frame[i] = row;
    }
    Vec3 row2 = cross( frame[0], frame[1] );
    frame[2] = normalize( row2 );
    return frame;
  }

  Vec3 normal{ values[0], values[1], values[2] };
  frame[0] = normalize( normal );

  Vec3 tangent = tangentHint;
  real64 projection = dot( tangent, frame[0] );
  for( int c = 0; c < 3; ++c )
  {
    tangent[c] -= projection * frame[0][c];
  }
  if( norm( tangent ) <= 1.0e-12 )
  {
    Vec3 fallback = std::abs( frame[0][0] ) < 0.9 ? Vec3{ 1.0, 0.0, 0.0 } : Vec3{ 0.0, 1.0, 0.0 };
    projection = dot( fallback, frame[0] );
    for( int c = 0; c < 3; ++c )
    {
      fallback[c] -= projection * frame[0][c];
    }
    tangent = fallback;
  }
  frame[1] = normalize( tangent );
  frame[2] = normalize( cross( frame[0], frame[1] ) );
  return frame;
}

MaterialFrameUpdate resolveMaterialFrameUpdate( MaterialFrameUpdate requested, ContinuumBase const & model )
{
  if( requested != MaterialFrameUpdate::Auto )
  {
    return requested;
  }
  if( model.getCatalogName() == "Graphite" )
  {
    return MaterialFrameUpdate::Graphite;
  }
  if( model.hasWrapper( "isFiber" ) )
  {
    return MaterialFrameUpdate::Fiber;
  }
  if( model.hasWrapper( "materialDirection" ) )
  {
    return MaterialFrameUpdate::Rotation;
  }
  return MaterialFrameUpdate::Fixed;
}

Mat3 updateMaterialFrame( Mat3 const & referenceFrame,
                          Mat3 const & F,
                          Mat3 const & R,
                          MaterialFrameUpdate mode )
{
  Mat3 out = referenceFrame;
  Mat3 const cofF = cofactor( F );

  if( mode == MaterialFrameUpdate::Fixed )
  {
    return referenceFrame;
  }

  if( mode == MaterialFrameUpdate::Rotation )
  {
    for( int i = 0; i < 3; ++i )
    {
      out[i] = normalize( matrixVectorProduct( R, referenceFrame[i] ) );
    }
    return out;
  }

  if( mode == MaterialFrameUpdate::Fiber )
  {
    for( int i = 0; i < 3; ++i )
    {
      out[i] = normalize( matrixVectorProduct( F, referenceFrame[i] ) );
    }
    return out;
  }

  if( mode == MaterialFrameUpdate::Normal )
  {
    out[0] = normalize( matrixVectorProduct( cofF, referenceFrame[0] ) );
    Vec3 tangent = matrixVectorProduct( F, referenceFrame[1] );
    real64 const projection = dot( tangent, out[0] );
    for( int c = 0; c < 3; ++c )
    {
      tangent[c] -= projection * out[0][c];
    }
    out[1] = normalize( tangent );
    out[2] = normalize( cross( out[0], out[1] ) );
    return out;
  }

  if( mode == MaterialFrameUpdate::Graphite )
  {
    out[0] = normalize( matrixVectorProduct( cofF, referenceFrame[0] ) );
    Vec3 tangent = matrixVectorProduct( F, referenceFrame[1] );
    real64 const projection = dot( tangent, out[0] );
    for( int c = 0; c < 3; ++c )
    {
      tangent[c] -= projection * out[0][c];
    }
    out[1] = normalize( tangent );
    out[2] = normalize( cross( out[0], out[1] ) );
    return out;
  }

  if( mode == MaterialFrameUpdate::MpmCofactor )
  {
    for( int i = 0; i < 3; ++i )
    {
      out[i] = normalize( matrixVectorProduct( cofF, referenceFrame[i] ) );
    }
    return out;
  }

  throw std::runtime_error( "Unhandled material frame update mode" );
}

Mat3 velocityGradientFromStrainIncrement( Vec6 const & strainIncrement, real64 const dt )
{
  Mat3 L = {};
  L[0][0] = strainIncrement[0] / dt;
  L[1][1] = strainIncrement[1] / dt;
  L[2][2] = strainIncrement[2] / dt;
  L[1][2] = 0.5 * strainIncrement[3] / dt;
  L[2][1] = 0.5 * strainIncrement[3] / dt;
  L[0][2] = 0.5 * strainIncrement[4] / dt;
  L[2][0] = 0.5 * strainIncrement[4] / dt;
  L[0][1] = 0.5 * strainIncrement[5] / dt;
  L[1][0] = 0.5 * strainIncrement[5] / dt;
  return L;
}

Vec6 strainIncrementFromVelocityGradient( Mat3 const & L, real64 const dt )
{
  Vec6 strainIncrement = {};
  strainIncrement[0] = L[0][0] * dt;
  strainIncrement[1] = L[1][1] * dt;
  strainIncrement[2] = L[2][2] * dt;
  strainIncrement[3] = ( L[1][2] + L[2][1] ) * dt;
  strainIncrement[4] = ( L[0][2] + L[2][0] ) * dt;
  strainIncrement[5] = ( L[0][1] + L[1][0] ) * dt;
  return strainIncrement;
}

Mat3 integrateDeformationGradient( Mat3 const & F, Mat3 const & L, real64 const dt )
{
  Mat3 A = identityMatrix();
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      A[i][j] += dt * L[i][j];
    }
  }
  return matrixProduct( A, F );
}

void buildTrialKinematics( Options const & options,
                           PathStep const & step,
                           std::vector< real64 > const & unknowns,
                           Vec6 & strainIncrement,
                           Mat3 & L )
{
  strainIncrement = {};
  int unknownIndex = 0;
  for( int c = 0; c < 6; ++c )
  {
    switch( options.controlModes[c] )
    {
      case ControlMode::Stress:
        strainIncrement[c] = unknowns[static_cast< std::size_t >( unknownIndex++ )];
        break;
      case ControlMode::StrainIncrement:
        strainIncrement[c] = step.values[c];
        break;
      case ControlMode::StrainRate:
      case ControlMode::TrueStrainRate:
        strainIncrement[c] = step.values[c] * step.dt;
        break;
    }
  }
  L = velocityGradientFromStrainIncrement( strainIncrement, step.dt );

  // For the diagonal true-strain-rate controls, the velocity-gradient diagonal is the log strain rate.
  for( int c = 0; c < 3; ++c )
  {
    if( options.controlModes[c] == ControlMode::TrueStrainRate )
    {
      L[c][c] = step.values[c];
      strainIncrement[c] = step.values[c] * step.dt;
    }
  }
}

bool finiteLowerStrainBound( real64 const value )
{
  return std::isfinite( value );
}

bool finiteUpperStrainBound( real64 const value )
{
  return std::isfinite( value );
}

void enforceStressControlStrainBounds( Options const & options,
                                       DriverState const & stateOld,
                                       std::vector< real64 > & unknowns )
{
  int unknownIndex = 0;
  for( int c = 0; c < 6; ++c )
  {
    if( options.controlModes[c] != ControlMode::Stress )
    {
      continue;
    }
    if( static_cast< std::size_t >( unknownIndex ) >= unknowns.size() )
    {
      return;
    }

    real64 const currentTotalStrain = stateOld.totalStrain[static_cast< std::size_t >( c )];
    real64 lowerIncrement = -std::numeric_limits< real64 >::infinity();
    real64 upperIncrement =  std::numeric_limits< real64 >::infinity();
    real64 const lower = options.stressControlMinStrain[static_cast< std::size_t >( c )];
    real64 const upper = options.stressControlMaxStrain[static_cast< std::size_t >( c )];
    if( finiteLowerStrainBound( lower ) )
    {
      lowerIncrement = lower - currentTotalStrain;
    }
    if( finiteUpperStrainBound( upper ) )
    {
      upperIncrement = upper - currentTotalStrain;
    }
    if( lowerIncrement > upperIncrement )
    {
      lowerIncrement = upperIncrement;
    }

    real64 & unknown = unknowns[static_cast< std::size_t >( unknownIndex )];
    unknown = std::min( std::max( unknown, lowerIncrement ), upperIncrement );
    ++unknownIndex;
  }
}

std::vector< real64 > boundedStressControlUnknowns( Options const & options,
                                                    DriverState const & stateOld,
                                                    std::vector< real64 > const & unknowns )
{
  std::vector< real64 > bounded = unknowns;
  enforceStressControlStrainBounds( options, stateOld, bounded );
  return bounded;
}

void validateStressControlStrainBounds( Options const & options, DriverState const & state )
{
  real64 const tolerance = 256.0 * std::numeric_limits< real64 >::epsilon();
  for( int c = 0; c < 6; ++c )
  {
    if( options.controlModes[c] != ControlMode::Stress )
    {
      continue;
    }
    real64 const value = state.totalStrain[static_cast< std::size_t >( c )];
    real64 const lower = options.stressControlMinStrain[static_cast< std::size_t >( c )];
    real64 const upper = options.stressControlMaxStrain[static_cast< std::size_t >( c )];
    if( finiteLowerStrainBound( lower ) && value + tolerance < lower )
    {
      std::ostringstream message;
      message << "Stress-controlled component " << voigtComponentName( c )
              << " has total strain " << value << " below minimum " << lower;
      throw std::runtime_error( message.str() );
    }
    if( finiteUpperStrainBound( upper ) && value - tolerance > upper )
    {
      std::ostringstream message;
      message << "Stress-controlled component " << voigtComponentName( c )
              << " has total strain " << value << " above maximum " << upper;
      throw std::runtime_error( message.str() );
    }
  }
}

real64 stressPower( Vec6 const & stress, Mat3 const & L )
{
  return stress[0] * L[0][0]
       + stress[1] * L[1][1]
       + stress[2] * L[2][2]
       + stress[3] * ( L[1][2] + L[2][1] )
       + stress[4] * ( L[0][2] + L[2][0] )
       + stress[5] * ( L[0][1] + L[1][0] );
}

std::vector< real64 > flatten( Mat3 const & A )
{
  std::vector< real64 > values;
  values.reserve( 9 );
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      values.push_back( A[i][j] );
    }
  }
  return values;
}

std::vector< real64 > flatten( Vec6 const & v )
{
  return std::vector< real64 >( v.begin(), v.end() );
}

struct RepositorySnapshot
{
  std::vector< std::function< void() > > restoreFunctions;

  void capture( Group & group )
  {
    restoreFunctions.clear();
    group.forWrappers( [&]( WrapperBase & wrapper )
    {
      if( wrapper.numArrayDims() == 0 )
      {
        return;
      }
      types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
      {
        using ArrayType = camp::first< decltype( tupleOfTypes ) >;
        if constexpr ( ArrayType::NDIM > 0 )
        {
          ArrayType const & array = Wrapper< ArrayType >::cast( wrapper ).reference();
          auto snapshot = std::make_shared< ArrayType >( array );
          WrapperBase * wrapperPtr = &wrapper;
          restoreFunctions.push_back( [snapshot, wrapperPtr]()
          {
            Wrapper< ArrayType >::cast( *wrapperPtr ).reference() = *snapshot;
          } );
        }
      }, wrapper );
    } );
  }

  void restore() const
  {
    for( auto const & restoreFunction : restoreFunctions )
    {
      restoreFunction();
    }
  }
};

template< typename ArrayType >
void setArrayFirstEntry( ArrayType & array,
                         std::vector< real64 > const & values,
                         std::string const & fieldName )
{
  using T = typename ArrayType::ValueType;
  static_assert( ArrayType::NDIM > 0, "setArrayFirstEntry expects an array" );
  auto view = array.toView();
  if( view.size( 0 ) < 1 )
  {
    throw std::runtime_error( "Field " + fieldName + " has zero parent size" );
  }

  localIndex entries = 1;
  for( int dim = 1; dim < ArrayType::NDIM; ++dim )
  {
    entries *= view.size( dim );
  }
  if( values.size() != 1 && values.size() != static_cast< std::size_t >( entries ) )
  {
    throw std::runtime_error( "Field " + fieldName + " expected 1 or " + std::to_string( entries ) + " values" );
  }

  auto valueAt = [&]( localIndex const flatIndex ) -> T
  {
    std::size_t const valueIndex = values.size() == 1 ? 0 : static_cast< std::size_t >( flatIndex );
    return static_cast< T >( values[valueIndex] );
  };

  if constexpr ( ArrayType::NDIM == 1 )
  {
    view[0] = valueAt( 0 );
  }
  else if constexpr ( ArrayType::NDIM == 2 )
  {
    localIndex flat = 0;
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      view[0][i] = valueAt( flat++ );
    }
  }
  else if constexpr ( ArrayType::NDIM == 3 )
  {
    localIndex flat = 0;
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      for( localIndex j = 0; j < view.size( 2 ); ++j )
      {
        view[0][i][j] = valueAt( flat++ );
      }
    }
  }
  else if constexpr ( ArrayType::NDIM == 4 )
  {
    localIndex flat = 0;
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      for( localIndex j = 0; j < view.size( 2 ); ++j )
      {
        for( localIndex l = 0; l < view.size( 3 ); ++l )
        {
          view[0][i][j][l] = valueAt( flat++ );
        }
      }
    }
  }
  else
  {
    throw std::runtime_error( "Unsupported array rank for field " + fieldName );
  }
}

void setFieldFirstEntryIfPresent( ContinuumBase & model,
                                  std::string const & fieldName,
                                  std::vector< real64 > const & values,
                                  bool const requireField = false )
{
  if( !model.hasWrapper( fieldName ) )
  {
    if( requireField )
    {
      throw std::runtime_error( "Material " + model.getName() + " does not have field " + fieldName );
    }
    return;
  }

  WrapperBase & wrapper = model.getWrapperBase( fieldName );
  if( wrapper.numArrayDims() == 0 )
  {
    if( requireField )
    {
      throw std::runtime_error( "Field " + fieldName + " is scalar; this driver only initializes per-point arrays" );
    }
    return;
  }

  bool handled = false;
  types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    using T = typename ArrayType::ValueType;
    if constexpr ( ArrayType::NDIM > 0 && std::is_arithmetic< T >::value )
    {
      setArrayFirstEntry( Wrapper< ArrayType >::cast( wrapper ).reference(), values, fieldName );
      handled = true;
    }
  }, wrapper );

  if( requireField && !handled )
  {
    throw std::runtime_error( "Could not initialize field " + fieldName );
  }
}

template< typename ArrayType >
std::vector< real64 > getArrayFirstEntry( ArrayType const & array )
{
  auto view = array.toViewConst();
  std::vector< real64 > values;
  if( view.size( 0 ) < 1 )
  {
    return values;
  }
  if constexpr ( ArrayType::NDIM == 1 )
  {
    values.push_back( static_cast< real64 >( view[0] ) );
  }
  else if constexpr ( ArrayType::NDIM == 2 )
  {
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      values.push_back( static_cast< real64 >( view[0][i] ) );
    }
  }
  else if constexpr ( ArrayType::NDIM == 3 )
  {
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      for( localIndex j = 0; j < view.size( 2 ); ++j )
      {
        values.push_back( static_cast< real64 >( view[0][i][j] ) );
      }
    }
  }
  else if constexpr ( ArrayType::NDIM == 4 )
  {
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      for( localIndex j = 0; j < view.size( 2 ); ++j )
      {
        for( localIndex l = 0; l < view.size( 3 ); ++l )
        {
          values.push_back( static_cast< real64 >( view[0][i][j][l] ) );
        }
      }
    }
  }
  return values;
}

std::vector< real64 > getFieldFirstEntryIfPresent( ContinuumBase const & model, std::string const & fieldName )
{
  if( !model.hasWrapper( fieldName ) )
  {
    return {};
  }
  WrapperBase const & wrapper = model.getWrapperBase( fieldName );
  if( wrapper.numArrayDims() == 0 )
  {
    return {};
  }

  std::vector< real64 > values;
  types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    using T = typename ArrayType::ValueType;
    if constexpr ( ArrayType::NDIM > 0 && std::is_arithmetic< T >::value )
    {
      values = getArrayFirstEntry( Wrapper< ArrayType >::cast( wrapper ).reference() );
    }
  }, wrapper );
  return values;
}

void setMaterialDirectionFieldIfPresent( ContinuumBase & model, Mat3 const & materialFrame )
{
  if( !model.hasWrapper( "materialDirection" ) )
  {
    return;
  }

  WrapperBase & wrapper = model.getWrapperBase( "materialDirection" );
  if( wrapper.numArrayDims() == 2 )
  {
    setFieldFirstEntryIfPresent( model, "materialDirection", std::vector< real64 >{ materialFrame[0][0],
                                                                                     materialFrame[0][1],
                                                                                     materialFrame[0][2] },
                                  true );
  }
  else if( wrapper.numArrayDims() == 3 )
  {
    setFieldFirstEntryIfPresent( model, "materialDirection", flatten( materialFrame ), true );
  }
  else
  {
    throw std::runtime_error( "Unsupported materialDirection rank " + std::to_string( wrapper.numArrayDims() ) );
  }
}

void pushDependencies( ContinuumBase & model, DriverState const & state, Mat3 const & L )
{
  setMaterialDirectionFieldIfPresent( model, state.materialFrame );

  setFieldFirstEntryIfPresent( model, "deformationGradient", flatten( state.F ) );
  setFieldFirstEntryIfPresent( model, "velocityGradient", flatten( L ) );
  setFieldFirstEntryIfPresent( model, "temperature", std::vector< real64 >{ state.temperature } );
  setFieldFirstEntryIfPresent( model, "temperatureRate", std::vector< real64 >{ state.temperatureRate } );
  setFieldFirstEntryIfPresent( model, "density", std::vector< real64 >{ state.density } );
  setFieldFirstEntryIfPresent( model, "jacobian", std::vector< real64 >{ state.jacobian } );
  setFieldFirstEntryIfPresent( model, "lengthScale", std::vector< real64 >{ state.lengthScale } );
  setFieldFirstEntryIfPresent( model, "strengthScale", std::vector< real64 >{ state.strengthScale } );
  setFieldFirstEntryIfPresent( model, "specificInternalEnergy", std::vector< real64 >{ state.specificInternalEnergy } );
  setFieldFirstEntryIfPresent( model, "internalEnergy", std::vector< real64 >{ state.specificInternalEnergy } );
}

void callConstitutiveUpdate( ContinuumBase & model,
                             DriverState const & stateOld,
                             TrialResult & result,
                             real64 const dt,
                             bool const commit )
{
  real64 fNew[3][3] = {};
  real64 rotBeginning[3][3] = {};
  real64 rotEnd[3][3] = {};
  real64 ddt[6] = {};
  real64 stress[6] = {};

  Mat3 const Rold = polarRotation( stateOld.F );
  Mat3 const Rnew = polarRotation( result.state.F );
  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      fNew[i][j] = result.state.F[i][j];
      rotBeginning[i][j] = Rold[i][j];
      rotEnd[i][j] = Rnew[i][j];
    }
  }
  for( int i = 0; i < 6; ++i )
  {
    ddt[i] = result.strainIncrement[i];
  }

  bool executed = false;
  bool const hyperelasticUpdate = model.getCatalogName() == "HyperelasticMMS" ||
                                  model.getCatalogName() == "Hyperelastic" ||
                                  model.getCatalogName() == "Chiumenti";

  ConstitutivePassThruMPM< ContinuumBase >::execute( model, [&]( auto & castedConstitutiveModel )
  {
    auto constitutiveWrapper = castedConstitutiveModel.createKernelUpdates();
    if( hyperelasticUpdate )
    {
      real64 FminusI[3][3] = {};
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          FminusI[i][j] = fNew[i][j];
        }
      }
      LvArray::tensorOps::addIdentity< 3 >( FminusI, -1.0 );
      constitutiveWrapper.hyperUpdate( k, q, FminusI, stress );
    }
    else
    {
      SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                              k,
                                              q,
                                              dt,
                                              ddt,
                                              rotBeginning,
                                              rotEnd,
                                              stress );
    }

    if( commit )
    {
      constitutiveWrapper.saveConvergedState( k, q );
    }
    executed = true;
  } );

  if( !executed )
  {
    throw std::runtime_error( "The material was not handled by ConstitutivePassThruMPM" );
  }

  for( int i = 0; i < 6; ++i )
  {
    result.state.stress[i] = stress[i];
  }
}

TrialResult evaluateTrial( ContinuumBase & model,
                           DriverState const & stateOld,
                           RepositorySnapshot const & snapshot,
                           Options const & options,
                           PathStep const & step,
                           std::vector< real64 > const & unknowns,
                           MaterialFrameUpdate const materialFrameUpdate,
                           bool const commit )
{
  snapshot.restore();

  TrialResult result;
  result.state = stateOld;
  std::vector< real64 > const boundedUnknowns = boundedStressControlUnknowns( options, stateOld, unknowns );
  buildTrialKinematics( options, step, boundedUnknowns, result.strainIncrement, result.L );

  result.state.F = integrateDeformationGradient( stateOld.F, result.L, step.dt );
  result.state.jacobian = determinant( result.state.F );
  if( result.state.jacobian <= 0.0 )
  {
    throw std::runtime_error( "Trial deformation gradient has non-positive determinant" );
  }
  result.state.density = result.state.referenceDensity / result.state.jacobian;
  result.state.materialFrame = updateMaterialFrame( stateOld.materialFrameReference,
                                                    result.state.F,
                                                    polarRotation( result.state.F ),
                                                    materialFrameUpdate );
  for( int c = 0; c < 6; ++c )
  {
    result.state.totalStrain[c] += result.strainIncrement[c];
    result.state.previousStrainIncrement[c] = result.strainIncrement[c];
  }
  validateStressControlStrainBounds( options, result.state );

  // The constitutive model sees the beginning-of-step temperature by default.
  // Adiabatic temperature increments are applied after the accepted stress update.
  if( options.temperatureMode == TemperatureMode::Isothermal ||
      options.temperatureMode == TemperatureMode::Prescribed )
  {
    result.state.temperature = options.initialTemperature;
    result.state.temperatureRate = 0.0;
  }

  pushDependencies( model, result.state, result.L );
  callConstitutiveUpdate( model, stateOld, result, step.dt, commit );

  if( commit )
  {
    if( options.energyMode == EnergyMode::StressPower )
    {
      real64 const oldPower = stressPower( stateOld.stress, result.L );
      real64 const newPower = stressPower( result.state.stress, result.L );
      real64 const workIncrement = 0.5 * step.dt * ( oldPower + newPower ) / result.state.density;
      result.state.lastStressPower = newPower;
      result.state.accumulatedStressPower += workIncrement;
      result.state.specificInternalEnergy += options.retentionFactor * workIncrement;

      if( options.temperatureMode == TemperatureMode::Adiabatic )
      {
        real64 const oldTemperature = result.state.temperature;
        result.state.temperature += options.retentionFactor * workIncrement / options.heatCapacity;
        result.state.temperatureRate = ( result.state.temperature - oldTemperature ) / step.dt;
        setFieldFirstEntryIfPresent( model, "temperature", std::vector< real64 >{ result.state.temperature } );
        setFieldFirstEntryIfPresent( model, "temperatureRate", std::vector< real64 >{ result.state.temperatureRate } );
      }
      setFieldFirstEntryIfPresent( model, "specificInternalEnergy", std::vector< real64 >{ result.state.specificInternalEnergy } );
      setFieldFirstEntryIfPresent( model, "internalEnergy", std::vector< real64 >{ result.state.specificInternalEnergy } );
    }
    else if( options.energyMode == EnergyMode::Material )
    {
      std::vector< real64 > values = getFieldFirstEntryIfPresent( model, "specificInternalEnergy" );
      if( values.empty() )
      {
        values = getFieldFirstEntryIfPresent( model, "internalEnergy" );
      }
      if( !values.empty() )
      {
        result.state.specificInternalEnergy = values[0];
      }
      if( options.temperatureMode == TemperatureMode::FromMaterial )
      {
        std::vector< real64 > temperature = getFieldFirstEntryIfPresent( model, "temperature" );
        if( !temperature.empty() )
        {
          result.state.temperature = temperature[0];
        }
      }
    }
  }

  return result;
}

std::vector< real64 > stressResidual( Options const & options, PathStep const & step, Vec6 const & stress )
{
  std::vector< real64 > residual;
  for( int c = 0; c < 6; ++c )
  {
    if( options.controlModes[c] == ControlMode::Stress )
    {
      residual.push_back( stress[c] - step.values[c] );
    }
  }
  return residual;
}

real64 l2Norm( std::vector< real64 > const & values )
{
  real64 sum = 0.0;
  for( real64 const value : values )
  {
    sum += value * value;
  }
  return std::sqrt( sum );
}

std::string csvEscape( std::string const & text )
{
  bool needsQuotes = text.empty();
  for( char const c : text )
  {
    needsQuotes = needsQuotes || c == ',' || c == '"' || c == '\n' || c == '\r';
  }
  if( !needsQuotes )
  {
    return text;
  }

  std::string escaped;
  escaped.reserve( text.size() + 2 );
  escaped.push_back( '"' );
  for( char const c : text )
  {
    if( c == '"' )
    {
      escaped.push_back( '"' );
      escaped.push_back( '"' );
    }
    else if( c == '\n' || c == '\r' )
    {
      escaped.push_back( ' ' );
    }
    else
    {
      escaped.push_back( c );
    }
  }
  escaped.push_back( '"' );
  return escaped;
}

void appendDiagnosticValue( std::ofstream & output, real64 const value )
{
  output << ',' << std::setprecision( 17 ) << value;
}

void appendDiagnosticString( std::ofstream & output, std::string const & value )
{
  output << ',' << csvEscape( value );
}

Vec6 nanVec6()
{
  real64 const nan = std::numeric_limits< real64 >::quiet_NaN();
  return Vec6{ nan, nan, nan, nan, nan, nan };
}

Mat3 nanMat3()
{
  real64 const nan = std::numeric_limits< real64 >::quiet_NaN();
  return Mat3{ Vec3{ nan, nan, nan }, Vec3{ nan, nan, nan }, Vec3{ nan, nan, nan } };
}

real64 pressureFromStress( Vec6 const & stress )
{
  return -( stress[0] + stress[1] + stress[2] ) / 3.0;
}

real64 firstFieldValueOrNaN( ContinuumBase const & model, std::string const & fieldName )
{
  std::vector< real64 > const values = getFieldFirstEntryIfPresent( model, fieldName );
  if( values.empty() )
  {
    return std::numeric_limits< real64 >::quiet_NaN();
  }
  return values[0];
}

std::vector< real64 > fieldValuesOrNaN( ContinuumBase const & model,
                                        std::string const & fieldName,
                                        std::size_t const count )
{
  std::vector< real64 > values = getFieldFirstEntryIfPresent( model, fieldName );
  values.resize( count, std::numeric_limits< real64 >::quiet_NaN() );
  return values;
}

struct StressControlDiagnostics
{
  StressControlDiagnosticsLevel level = StressControlDiagnosticsLevel::Off;
  std::ofstream output;
  localIndex recordCount = 0;

  bool active() const
  {
    return output.is_open() && level != StressControlDiagnosticsLevel::Off;
  }

  bool shouldRecord( StressControlDiagnosticsLevel const requested ) const
  {
    return active() && stressControlDiagnosticsAtLeast( level, requested );
  }

  void open( Options const & options )
  {
    level = options.stressControlDiagnosticsLevel;
    if( options.stressControlDiagnosticsPath.empty() || level == StressControlDiagnosticsLevel::Off )
    {
      return;
    }
    output.open( options.stressControlDiagnosticsPath );
    if( !output )
    {
      throw std::runtime_error( "Could not open stress-control diagnostic CSV: " + options.stressControlDiagnosticsPath );
    }
    writeHeader();
  }

  void writeHeader()
  {
    std::array< std::string, 6 > const voigt = { "xx", "yy", "zz", "yz", "xz", "xy" };
    output << "record,step,timeOld,timeTrial,dt,stage,iteration,subiteration,accepted,converged,residualNorm,message";
    for( std::string const & name : voigt ) output << ",startStress_" << name;
    output << ",startPressure,trialPressure";
    for( std::string const & name : voigt ) output << ",unknown_" << name;
    for( std::string const & name : voigt ) output << ",strainIncrement_" << name;
    for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) output << ",L_" << i << j;
    for( std::string const & name : voigt ) output << ",stress_" << name;
    for( std::string const & name : voigt ) output << ",target_" << name;
    for( std::string const & name : voigt ) output << ",residual_" << name;
    for( std::string const & name : voigt ) output << ",totalStrain_" << name;
    for( std::string const & name : voigt ) output << ",strainMin_" << name;
    for( std::string const & name : voigt ) output << ",strainMax_" << name;
    for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) output << ",F_" << i << j;
    for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) output << ",materialDirection_" << i << j;
    output << ",density,jacobian,temperature,effectiveBulkModulus,effectiveShearModulus,damage,basalPlaneDamage,comminutionDamage";
    for( std::string const & name : voigt ) output << ",plasticStrain_" << name;
    output << '\n';
  }

  void record( ContinuumBase const & model,
               localIndex const stepIndex,
               DriverState const & stateOld,
               Options const & options,
               PathStep const & step,
               std::string const & stage,
               int const iteration,
               int const subiteration,
               std::vector< real64 > const & unknowns,
               TrialResult const * const trial,
               bool const accepted,
               std::string const & message,
               StressControlDiagnosticsLevel const requested = StressControlDiagnosticsLevel::Iteration )
  {
    if( !shouldRecord( requested ) )
    {
      return;
    }

    std::vector< real64 > const boundedUnknowns = boundedStressControlUnknowns( options, stateOld, unknowns );
    Vec6 unknownByComponent = nanVec6();
    int unknownIndex = 0;
    for( int c = 0; c < 6; ++c )
    {
      if( options.controlModes[c] == ControlMode::Stress )
      {
        if( static_cast< std::size_t >( unknownIndex ) < boundedUnknowns.size() )
        {
          unknownByComponent[c] = boundedUnknowns[static_cast< std::size_t >( unknownIndex )];
        }
        ++unknownIndex;
      }
    }

    Vec6 strainIncrement = nanVec6();
    Mat3 L = nanMat3();
    try
    {
      buildTrialKinematics( options, step, boundedUnknowns, strainIncrement, L );
    }
    catch( std::exception const & )
    {
      strainIncrement = nanVec6();
      L = nanMat3();
    }

    Vec6 stress = nanVec6();
    Vec6 target = nanVec6();
    Vec6 residualByComponent = nanVec6();
    Vec6 totalStrain = nanVec6();
    Mat3 F = nanMat3();
    Mat3 materialFrame = nanMat3();
    real64 residualNorm = std::numeric_limits< real64 >::quiet_NaN();
    real64 trialPressure = std::numeric_limits< real64 >::quiet_NaN();
    int converged = -1;
    real64 density = std::numeric_limits< real64 >::quiet_NaN();
    real64 jacobian = std::numeric_limits< real64 >::quiet_NaN();
    real64 temperature = std::numeric_limits< real64 >::quiet_NaN();

    if( trial != nullptr )
    {
      stress = trial->state.stress;
      trialPressure = pressureFromStress( stress );
      totalStrain = trial->state.totalStrain;
      F = trial->state.F;
      materialFrame = trial->state.materialFrame;
      residualNorm = l2Norm( stressResidual( options, step, stress ) );
      converged = residualNorm <= options.stressTolerance ? 1 : 0;
      density = trial->state.density;
      jacobian = trial->state.jacobian;
      temperature = trial->state.temperature;
    }

    for( int c = 0; c < 6; ++c )
    {
      if( options.controlModes[c] == ControlMode::Stress )
      {
        target[c] = step.values[c];
        if( trial != nullptr )
        {
          residualByComponent[c] = stress[c] - target[c];
        }
      }
    }

    output << recordCount++;
    appendDiagnosticValue( output, static_cast< real64 >( stepIndex ) );
    appendDiagnosticValue( output, stateOld.time );
    appendDiagnosticValue( output, stateOld.time + step.dt );
    appendDiagnosticValue( output, step.dt );
    appendDiagnosticString( output, stage );
    output << ',' << iteration << ',' << subiteration << ',' << ( accepted ? 1 : 0 ) << ',' << converged;
    appendDiagnosticValue( output, residualNorm );
    appendDiagnosticString( output, message );
    for( real64 const value : stateOld.stress ) appendDiagnosticValue( output, value );
    appendDiagnosticValue( output, pressureFromStress( stateOld.stress ) );
    appendDiagnosticValue( output, trialPressure );
    for( real64 const value : unknownByComponent ) appendDiagnosticValue( output, value );
    for( real64 const value : strainIncrement ) appendDiagnosticValue( output, value );
    for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) appendDiagnosticValue( output, L[i][j] );
    for( real64 const value : stress ) appendDiagnosticValue( output, value );
    for( real64 const value : target ) appendDiagnosticValue( output, value );
    for( real64 const value : residualByComponent ) appendDiagnosticValue( output, value );
    for( real64 const value : totalStrain ) appendDiagnosticValue( output, value );
    for( real64 const value : options.stressControlMinStrain ) appendDiagnosticValue( output, value );
    for( real64 const value : options.stressControlMaxStrain ) appendDiagnosticValue( output, value );
    for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) appendDiagnosticValue( output, F[i][j] );
    for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) appendDiagnosticValue( output, materialFrame[i][j] );
    appendDiagnosticValue( output, density );
    appendDiagnosticValue( output, jacobian );
    appendDiagnosticValue( output, temperature );
    appendDiagnosticValue( output, firstFieldValueOrNaN( model, "effectiveBulkModulus" ) );
    appendDiagnosticValue( output, firstFieldValueOrNaN( model, "effectiveShearModulus" ) );
    appendDiagnosticValue( output, firstFieldValueOrNaN( model, "damage" ) );
    appendDiagnosticValue( output, firstFieldValueOrNaN( model, "basalPlaneDamage" ) );
    appendDiagnosticValue( output, firstFieldValueOrNaN( model, "comminutionDamage" ) );
    std::vector< real64 > const plasticStrain = fieldValuesOrNaN( model, "plasticStrain", 6 );
    for( real64 const value : plasticStrain ) appendDiagnosticValue( output, value );
    output << '\n';
    output.flush();
  }
};

void recordStressControlDiagnostic( StressControlDiagnostics * const diagnostics,
                                    ContinuumBase const & model,
                                    localIndex const stepIndex,
                                    DriverState const & stateOld,
                                    Options const & options,
                                    PathStep const & step,
                                    std::string const & stage,
                                    int const iteration,
                                    int const subiteration,
                                    std::vector< real64 > const & unknowns,
                                    TrialResult const * const trial,
                                    bool const accepted,
                                    std::string const & message,
                                    StressControlDiagnosticsLevel const requested = StressControlDiagnosticsLevel::Iteration )
{
  if( diagnostics != nullptr )
  {
    diagnostics->record( model,
                         stepIndex,
                         stateOld,
                         options,
                         step,
                         stage,
                         iteration,
                         subiteration,
                         unknowns,
                         trial,
                         accepted,
                         message,
                         requested );
  }
}

std::vector< real64 > solveLinearSystem( std::vector< std::vector< real64 > > A,
                                         std::vector< real64 > b )
{
  std::size_t const n = b.size();
  for( std::size_t i = 0; i < n; ++i )
  {
    std::size_t pivot = i;
    real64 pivotAbs = std::abs( A[i][i] );
    for( std::size_t r = i + 1; r < n; ++r )
    {
      real64 const candidateAbs = std::abs( A[r][i] );
      if( candidateAbs > pivotAbs )
      {
        pivot = r;
        pivotAbs = candidateAbs;
      }
    }
    if( pivotAbs <= 1.0e-30 )
    {
      throw std::runtime_error( "Singular stress-control Jacobian" );
    }
    if( pivot != i )
    {
      std::swap( A[i], A[pivot] );
      std::swap( b[i], b[pivot] );
    }
    real64 const diagonal = A[i][i];
    for( std::size_t c = i; c < n; ++c )
    {
      A[i][c] /= diagonal;
    }
    b[i] /= diagonal;
    for( std::size_t r = 0; r < n; ++r )
    {
      if( r == i )
      {
        continue;
      }
      real64 const factor = A[r][i];
      for( std::size_t c = i; c < n; ++c )
      {
        A[r][c] -= factor * A[i][c];
      }
      b[r] -= factor * b[i];
    }
  }
  return b;
}

bool isFinite( real64 const value )
{
  return std::isfinite( value );
}

bool isFinite( std::vector< real64 > const & values )
{
  for( real64 const value : values )
  {
    if( !isFinite( value ) )
    {
      return false;
    }
  }
  return true;
}

std::vector< real64 > solveRegularizedLeastSquares( std::vector< std::vector< real64 > > const & jacobian,
                                                    std::vector< real64 > const & rhs,
                                                    real64 const requestedRegularization )
{
  std::size_t const n = rhs.size();
  std::vector< std::vector< real64 > > normal( n, std::vector< real64 >( n, 0.0 ) );
  std::vector< real64 > normalRhs( n, 0.0 );

  for( std::size_t row = 0; row < n; ++row )
  {
    for( std::size_t i = 0; i < n; ++i )
    {
      normalRhs[i] += jacobian[row][i] * rhs[row];
      for( std::size_t j = 0; j < n; ++j )
      {
        normal[i][j] += jacobian[row][i] * jacobian[row][j];
      }
    }
  }

  real64 diagonalScale = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    diagonalScale = std::max( diagonalScale, std::abs( normal[i][i] ) );
  }
  real64 const regularization = requestedRegularization > 0.0 ?
                                requestedRegularization :
                                1.0e-12 * std::max( real64( 1.0 ), diagonalScale );
  for( std::size_t i = 0; i < n; ++i )
  {
    normal[i][i] += regularization;
  }
  return solveLinearSystem( normal, normalRhs );
}

void limitStressControlStepDirection( Options const & options, std::vector< real64 > & stepDirection )
{
  real64 const limit = options.stressControlMaxStrainCorrection;
  if( limit <= 0.0 )
  {
    return;
  }
  for( real64 & value : stepDirection )
  {
    value = std::max( -limit, std::min( limit, value ) );
  }
}

bool updateBestStressControlCandidate( ContinuumBase & model,
                                       DriverState const & stateOld,
                                       RepositorySnapshot const & snapshot,
                                       Options const & options,
                                       PathStep const & step,
                                       MaterialFrameUpdate const materialFrameUpdate,
                                       std::vector< real64 > const & candidateUnknowns,
                                       std::vector< real64 > & bestUnknowns,
                                       real64 & bestResidualNorm,
                                       bool & haveBest,
                                       real64 & candidateResidualNorm,
                                       StressControlDiagnostics * const diagnostics,
                                       localIndex const stepIndex,
                                       std::string const & stage,
                                       int const iteration,
                                       int const subiteration,
                                       StressControlDiagnosticsLevel const diagnosticsLevel )
{
  candidateResidualNorm = std::numeric_limits< real64 >::infinity();
  std::vector< real64 > const boundedUnknowns = boundedStressControlUnknowns( options, stateOld, candidateUnknowns );
  try
  {
    TrialResult trial = evaluateTrial( model, stateOld, snapshot, options, step, boundedUnknowns, materialFrameUpdate, false );
    recordStressControlDiagnostic( diagnostics,
                                   model,
                                   stepIndex,
                                   stateOld,
                                   options,
                                   step,
                                   stage,
                                   iteration,
                                   subiteration,
                                   boundedUnknowns,
                                   &trial,
                                   false,
                                   std::string(),
                                   diagnosticsLevel );
    std::vector< real64 > const residual = stressResidual( options, step, trial.state.stress );
    candidateResidualNorm = l2Norm( residual );
    if( isFinite( candidateResidualNorm ) && candidateResidualNorm < bestResidualNorm )
    {
      bestResidualNorm = candidateResidualNorm;
      bestUnknowns = boundedUnknowns;
      haveBest = true;
    }
    return isFinite( candidateResidualNorm ) && candidateResidualNorm <= options.stressTolerance;
  }
  catch( std::exception const & error )
  {
    recordStressControlDiagnostic( diagnostics,
                                   model,
                                   stepIndex,
                                   stateOld,
                                   options,
                                   step,
                                   stage,
                                   iteration,
                                   subiteration,
                                   boundedUnknowns,
                                   nullptr,
                                   false,
                                   error.what(),
                                   diagnosticsLevel );
    return false;
  }

}

real64 characteristicPrescribedStrainIncrement( Options const & options, PathStep const & step );

std::vector< std::vector< real64 > > stressControlPatternDirections( std::size_t const numUnknowns )
{
  std::vector< std::vector< real64 > > directions;
  if( numUnknowns == 0 )
  {
    return directions;
  }

  for( std::size_t j = 0; j < numUnknowns; ++j )
  {
    std::vector< real64 > positive( numUnknowns, 0.0 );
    positive[j] = 1.0;
    directions.push_back( positive );

    std::vector< real64 > negative( numUnknowns, 0.0 );
    negative[j] = -1.0;
    directions.push_back( negative );
  }

  if( numUnknowns <= 4 )
  {
    std::size_t const numSignPatterns = std::size_t( 1 ) << numUnknowns;
    real64 const scale = 1.0 / std::sqrt( static_cast< real64 >( numUnknowns ) );
    for( std::size_t mask = 0; mask < numSignPatterns; ++mask )
    {
      std::vector< real64 > direction( numUnknowns, 0.0 );
      for( std::size_t j = 0; j < numUnknowns; ++j )
      {
        bool const positive = ( mask & ( std::size_t( 1 ) << j ) ) != 0;
        direction[j] = positive ? scale : -scale;
      }
      directions.push_back( direction );
    }
  }
  else
  {
    real64 const scale = 1.0 / std::sqrt( static_cast< real64 >( numUnknowns ) );
    std::vector< real64 > positive( numUnknowns, scale );
    std::vector< real64 > negative( numUnknowns, -scale );
    directions.push_back( positive );
    directions.push_back( negative );
  }

  return directions;
}

bool patternStressControlFallback( ContinuumBase & model,
                                   DriverState const & stateOld,
                                   RepositorySnapshot const & snapshot,
                                   Options const & options,
                                   PathStep const & step,
                                   MaterialFrameUpdate const materialFrameUpdate,
                                   StressControlDiagnostics * const diagnostics,
                                   localIndex const stepIndex,
                                   std::vector< real64 > & bestUnknowns,
                                   real64 & bestResidualNorm,
                                   bool & haveBest,
                                   std::string & stopReason )
{
  if( bestUnknowns.empty() || options.stressControlPatternIterations <= 0 )
  {
    return false;
  }

  std::vector< real64 > unknowns = haveBest ? bestUnknowns :
                                   std::vector< real64 >( bestUnknowns.size(), 0.0 );
  enforceStressControlStrainBounds( options, stateOld, unknowns );
  std::vector< std::vector< real64 > > const directions = stressControlPatternDirections( unknowns.size() );
  if( directions.empty() )
  {
    return false;
  }

  real64 stepSize = options.stressControlPatternInitialStep > 0.0 ?
                    options.stressControlPatternInitialStep :
                    characteristicPrescribedStrainIncrement( options, step );
  stepSize = std::max( stepSize, 10.0 * options.finiteDifferenceEpsilon );
  real64 const minStepSize = std::max( options.stressControlPatternMinStep,
                                       10.0 * std::numeric_limits< real64 >::epsilon() );
  real64 const maxStepSize = std::max( options.stressControlMaxStrainCorrection,
                                       stepSize );
  bool improved = false;

  for( int iteration = 0; iteration < options.stressControlPatternIterations && stepSize >= minStepSize; ++iteration )
  {
    real64 currentNorm = std::numeric_limits< real64 >::infinity();
    bool const currentExact = updateBestStressControlCandidate( model,
                                                                stateOld,
                                                                snapshot,
                                                                options,
                                                                step,
                                                                materialFrameUpdate,
                                                                unknowns,
                                                                bestUnknowns,
                                                                bestResidualNorm,
                                                                haveBest,
                                                                currentNorm,
                                                                diagnostics,
                                                                stepIndex,
                                                                "pattern_current",
                                                                iteration,
                                                                0,
                                                                StressControlDiagnosticsLevel::Iteration );
    if( currentExact )
    {
      stopReason = "used stress-control pattern-search fallback";
      return true;
    }
    if( !isFinite( currentNorm ) )
    {
      currentNorm = bestResidualNorm;
    }

    real64 bestLocalNorm = currentNorm;
    std::vector< real64 > bestLocalUnknowns = unknowns;
    for( std::size_t directionIndex = 0; directionIndex < directions.size(); ++directionIndex )
    {
      std::vector< real64 > correction = directions[directionIndex];
      for( real64 & value : correction )
      {
        value *= stepSize;
      }
      limitStressControlStepDirection( options, correction );

      bool nonzeroCorrection = false;
      for( real64 const value : correction )
      {
        if( std::abs( value ) > 10.0 * std::numeric_limits< real64 >::epsilon() )
        {
          nonzeroCorrection = true;
          break;
        }
      }
      if( !nonzeroCorrection )
      {
        continue;
      }

      std::vector< real64 > candidate = unknowns;
      for( std::size_t j = 0; j < candidate.size(); ++j )
      {
        candidate[j] += correction[j];
      }
      enforceStressControlStrainBounds( options, stateOld, candidate );

      real64 candidateNorm = std::numeric_limits< real64 >::infinity();
      bool const exact = updateBestStressControlCandidate( model,
                                                           stateOld,
                                                           snapshot,
                                                           options,
                                                           step,
                                                           materialFrameUpdate,
                                                           candidate,
                                                           bestUnknowns,
                                                           bestResidualNorm,
                                                           haveBest,
                                                           candidateNorm,
                                                           diagnostics,
                                                           stepIndex,
                                                           "pattern_candidate",
                                                           iteration,
                                                           static_cast< int >( directionIndex ),
                                                           StressControlDiagnosticsLevel::Iteration );
      if( exact )
      {
        stopReason = "used stress-control pattern-search fallback";
        return true;
      }
      if( isFinite( candidateNorm ) && candidateNorm < bestLocalNorm )
      {
        bestLocalNorm = candidateNorm;
        bestLocalUnknowns = candidate;
      }
    }

    if( bestLocalNorm < currentNorm )
    {
      unknowns = bestLocalUnknowns;
      enforceStressControlStrainBounds( options, stateOld, unknowns );
      improved = true;
      stepSize = std::min( maxStepSize, stepSize * options.stressControlPatternGrowth );
    }
    else
    {
      stepSize *= options.stressControlPatternShrink;
    }
  }

  if( improved )
  {
    stopReason = "used best available stress-control pattern-search fallback trial";
  }
  return improved;
}

bool servoStressControlFallback( ContinuumBase & model,
                                 DriverState const & stateOld,
                                 RepositorySnapshot const & snapshot,
                                 Options const & options,
                                 PathStep const & step,
                                 MaterialFrameUpdate const materialFrameUpdate,
                                 StressControlDiagnostics * const diagnostics,
                                 localIndex const stepIndex,
                                 std::vector< real64 > & bestUnknowns,
                                 real64 & bestResidualNorm,
                                 bool & haveBest,
                                 std::string & stopReason )
{
  if( bestUnknowns.empty() || options.stressControlServoIterations <= 0 )
  {
    return false;
  }

  std::vector< real64 > unknowns = haveBest ? bestUnknowns :
                                   std::vector< real64 >( bestUnknowns.size(), 0.0 );
  enforceStressControlStrainBounds( options, stateOld, unknowns );
  bool improved = false;
  real64 const relaxation = std::min( real64( 1.0 ), std::max( real64( 0.0 ), options.stressControlServoRelaxation ) );
  real64 const compliance = std::max( options.stressControlServoCompliance, real64( 0.0 ) );
  real64 const derivativeFloor = std::max( options.stressControlServoDerivativeFloor,
                                           std::numeric_limits< real64 >::epsilon() );

  for( int iteration = 0; iteration < options.stressControlServoIterations; ++iteration )
  {
    try
    {
      TrialResult trial = evaluateTrial( model, stateOld, snapshot, options, step, unknowns, materialFrameUpdate, false );
      recordStressControlDiagnostic( diagnostics,
                                     model,
                                     stepIndex,
                                     stateOld,
                                     options,
                                     step,
                                     "servo_current",
                                     iteration,
                                     0,
                                     unknowns,
                                     &trial,
                                     false,
                                     std::string() );
      std::vector< real64 > const residual = stressResidual( options, step, trial.state.stress );
      real64 const residualNorm = l2Norm( residual );
      if( isFinite( residualNorm ) && residualNorm < bestResidualNorm )
      {
        bestResidualNorm = residualNorm;
        bestUnknowns = unknowns;
        haveBest = true;
        improved = true;
      }
      if( isFinite( residualNorm ) && residualNorm <= options.stressTolerance )
      {
        stopReason = "used stress-control servo fallback";
        return true;
      }
      if( !isFinite( residualNorm ) || !isFinite( residual ) )
      {
        stopReason = "stress-control servo fallback encountered a non-finite residual";
        return improved;
      }

      std::vector< real64 > correction( unknowns.size(), 0.0 );
      for( std::size_t j = 0; j < unknowns.size(); ++j )
      {
        real64 const dx = std::max( options.finiteDifferenceEpsilon,
                                    options.finiteDifferenceEpsilon * std::max( real64( 1.0 ), std::abs( unknowns[j] ) ) );
        real64 derivative = 0.0;
        bool haveDerivative = false;
        try
        {
          std::vector< real64 > perturbed = unknowns;
          perturbed[j] += dx;
          enforceStressControlStrainBounds( options, stateOld, perturbed );
          TrialResult perturbedResult = evaluateTrial( model,
                                                       stateOld,
                                                       snapshot,
                                                       options,
                                                       step,
                                                       perturbed,
                                                       materialFrameUpdate,
                                                       false );
          recordStressControlDiagnostic( diagnostics,
                                         model,
                                         stepIndex,
                                         stateOld,
                                         options,
                                         step,
                                         "servo_fd",
                                         iteration,
                                         static_cast< int >( j ),
                                         perturbed,
                                         &perturbedResult,
                                         false,
                                         std::string(),
                                         StressControlDiagnosticsLevel::Full );
          std::vector< real64 > const perturbedResidual = stressResidual( options, step, perturbedResult.state.stress );
          if( perturbedResidual.size() == residual.size() && isFinite( perturbedResidual ) )
          {
            derivative = ( perturbedResidual[j] - residual[j] ) / dx;
            haveDerivative = std::abs( derivative ) > derivativeFloor;
          }
        }
        catch( std::exception const & )
        {
          haveDerivative = false;
        }

        if( haveDerivative )
        {
          correction[j] = -relaxation * residual[j] / derivative;
        }
        else if( compliance > real64( 0.0 ) )
        {
          correction[j] = -relaxation * compliance * residual[j];
        }
      }
      limitStressControlStepDirection( options, correction );

      bool nonzeroCorrection = false;
      for( real64 const value : correction )
      {
        if( std::abs( value ) > 10.0 * std::numeric_limits< real64 >::epsilon() )
        {
          nonzeroCorrection = true;
          break;
        }
      }
      if( !nonzeroCorrection )
      {
        stopReason = "stress-control servo correction vanished";
        return improved;
      }

      real64 bestLocalNorm = residualNorm;
      std::vector< real64 > bestLocalUnknowns = unknowns;
      std::array< real64, 5 > const servoLineFactors = { 1.0, 0.5, 0.25, -0.5, -1.0 };
      for( std::size_t lineIndex = 0; lineIndex < servoLineFactors.size(); ++lineIndex )
      {
        real64 const sign = servoLineFactors[lineIndex];
        std::vector< real64 > candidate = unknowns;
        for( std::size_t j = 0; j < unknowns.size(); ++j )
        {
          candidate[j] += sign * correction[j];
        }
        real64 candidateNorm = std::numeric_limits< real64 >::infinity();
        bool const exact = updateBestStressControlCandidate( model,
                                                             stateOld,
                                                             snapshot,
                                                             options,
                                                             step,
                                                             materialFrameUpdate,
                                                             candidate,
                                                             bestUnknowns,
                                                             bestResidualNorm,
                                                             haveBest,
                                                             candidateNorm,
                                                             diagnostics,
                                                             stepIndex,
                                                             "servo_candidate",
                                                             iteration,
                                                             static_cast< int >( lineIndex ),
                                                             StressControlDiagnosticsLevel::Iteration );
        if( exact )
        {
          stopReason = "used stress-control servo fallback";
          return true;
        }
        if( isFinite( candidateNorm ) && candidateNorm < bestLocalNorm )
        {
          bestLocalNorm = candidateNorm;
          bestLocalUnknowns = candidate;
          improved = true;
        }
      }

      if( bestLocalNorm < residualNorm )
      {
        unknowns = bestLocalUnknowns;
        enforceStressControlStrainBounds( options, stateOld, unknowns );
      }
      else
      {
        break;
      }
    }
    catch( std::exception const & error )
    {
      stopReason = std::string( "stress-control servo fallback trial failed: " ) + error.what();
      return improved;
    }
  }

  if( improved )
  {
    stopReason = "used best available stress-control servo fallback trial";
  }
  return improved;
}

real64 characteristicPrescribedStrainIncrement( Options const & options, PathStep const & step )
{
  Vec6 strainIncrement = {};
  Mat3 L = {};
  int numUnknowns = 0;
  for( ControlMode const mode : options.controlModes )
  {
    if( mode == ControlMode::Stress )
    {
      ++numUnknowns;
    }
  }
  buildTrialKinematics( options, step, std::vector< real64 >( static_cast< std::size_t >( numUnknowns ), 0.0 ), strainIncrement, L );

  real64 scale = 0.0;
  for( int c = 0; c < 6; ++c )
  {
    if( options.controlModes[c] != ControlMode::Stress )
    {
      scale = std::max( scale, std::abs( strainIncrement[c] ) );
    }
  }
  return std::max( scale, 10.0 * options.finiteDifferenceEpsilon );
}

bool residualsBracketRoot( real64 const a, real64 const b )
{
  return ( a <= 0.0 && b >= 0.0 ) || ( a >= 0.0 && b <= 0.0 );
}

bool scalarStressControlFallback( ContinuumBase & model,
                                  DriverState const & stateOld,
                                  RepositorySnapshot const & snapshot,
                                  Options const & options,
                                  PathStep const & step,
                                  MaterialFrameUpdate const materialFrameUpdate,
                                  StressControlDiagnostics * const diagnostics,
                                  localIndex const stepIndex,
                                  std::vector< real64 > & bestUnknowns,
                                  real64 & bestResidualNorm,
                                  bool & haveBest,
                                  std::string & stopReason )
{
  if( bestUnknowns.size() != 1 )
  {
    return false;
  }

  struct Sample
  {
    real64 x;
    real64 r;
  };

  std::vector< Sample > samples;
  bool exact = false;

  auto evaluate = [&]( real64 const x,
                       real64 & residualValue,
                       std::string const & stage,
                       int const iteration,
                       int const subiteration ) -> bool
  {
    std::vector< real64 > trialUnknowns{ x };
    enforceStressControlStrainBounds( options, stateOld, trialUnknowns );
    real64 const boundedX = trialUnknowns[0];
    try
    {
      TrialResult trial = evaluateTrial( model, stateOld, snapshot, options, step, trialUnknowns, materialFrameUpdate, false );
      recordStressControlDiagnostic( diagnostics,
                                     model,
                                     stepIndex,
                                     stateOld,
                                     options,
                                     step,
                                     stage,
                                     iteration,
                                     subiteration,
                                     trialUnknowns,
                                     &trial,
                                     false,
                                     std::string() );
      std::vector< real64 > const residual = stressResidual( options, step, trial.state.stress );
      if( residual.size() != 1 || !isFinite( residual ) )
      {
        return false;
      }
      residualValue = residual[0];
      real64 const residualNorm = std::abs( residualValue );
      if( isFinite( residualNorm ) && residualNorm < bestResidualNorm )
      {
        bestResidualNorm = residualNorm;
        bestUnknowns[0] = boundedX;
        haveBest = true;
      }
      if( residualNorm <= options.stressTolerance )
      {
        exact = true;
      }
      samples.push_back( Sample{ boundedX, residualValue } );
      return true;
    }
    catch( std::exception const & error )
    {
      recordStressControlDiagnostic( diagnostics,
                                     model,
                                     stepIndex,
                                     stateOld,
                                     options,
                                     step,
                                     stage,
                                     iteration,
                                     subiteration,
                                     trialUnknowns,
                                     nullptr,
                                     false,
                                     error.what() );
      stopReason = std::string( "scalar stress-control fallback trial failed: " ) + error.what();
      return false;
    }
  };

  auto bisectBracket = [&]( Sample lo, Sample hi ) -> bool
  {
    if( std::abs( lo.r ) <= options.stressTolerance )
    {
      bestUnknowns[0] = lo.x;
      bestResidualNorm = std::abs( lo.r );
      haveBest = true;
      return true;
    }
    if( std::abs( hi.r ) <= options.stressTolerance )
    {
      bestUnknowns[0] = hi.x;
      bestResidualNorm = std::abs( hi.r );
      haveBest = true;
      return true;
    }

    for( int iteration = 0; iteration < options.maxStressBisectionIterations; ++iteration )
    {
      real64 const xMid = 0.5 * ( lo.x + hi.x );
      real64 rMid = 0.0;
      if( !evaluate( xMid, rMid, "scalar_bisection", iteration, 0 ) )
      {
        break;
      }
      real64 const boundedMid = samples.empty() ? xMid : samples.back().x;
      if( std::abs( rMid ) <= options.stressTolerance )
      {
        bestUnknowns[0] = boundedMid;
        bestResidualNorm = std::abs( rMid );
        haveBest = true;
        return true;
      }
      if( residualsBracketRoot( lo.r, rMid ) )
      {
        hi = Sample{ boundedMid, rMid };
      }
      else
      {
        lo = Sample{ boundedMid, rMid };
      }
    }
    return false;
  };

  real64 const x0 = haveBest ? bestUnknowns[0] : 0.0;
  real64 r0 = 0.0;
  evaluate( x0, r0, "scalar_initial", 0, 0 );
  if( exact )
  {
    return true;
  }

  real64 scale = options.stressBracketInitialScale > 0.0 ?
                 options.stressBracketInitialScale :
                 characteristicPrescribedStrainIncrement( options, step );
  scale = std::max( scale, 10.0 * options.finiteDifferenceEpsilon );
  real64 const maxScale = std::max( scale, options.stressBracketMaxStrain );
  real64 const growth = options.stressBracketGrowth > 1.0 ? options.stressBracketGrowth : 2.0;

  for( int iteration = 0; iteration < options.maxStressBracketIterations && scale <= maxScale; ++iteration )
  {
    for( real64 const sign : { -1.0, 1.0 } )
    {
      real64 const x = x0 + sign * scale;
      real64 r = 0.0;
      if( !evaluate( x, r, "scalar_bracket", iteration, sign < 0.0 ? 0 : 1 ) )
      {
        continue;
      }
      if( exact )
      {
        return true;
      }
      Sample const current{ samples.empty() ? x : samples.back().x, r };
      for( Sample const & previous : samples )
      {
        real64 const duplicateTolerance = std::numeric_limits< real64 >::epsilon() *
                                          std::max( real64( 1.0 ),
                                                    std::max( std::abs( previous.x ),
                                                              std::abs( current.x ) ) );
        if( std::abs( previous.x - current.x ) <= duplicateTolerance )
        {
          continue;
        }
        if( residualsBracketRoot( previous.r, current.r ) )
        {
          stopReason = "used scalar bracketed stress-control fallback";
          if( bisectBracket( previous, current ) )
          {
            return true;
          }
          return false;
        }
      }
    }
    scale *= growth;
  }

  if( haveBest )
  {
    stopReason = "scalar stress-control fallback could not bracket a root; using best available trial";
  }
  return false;
}

TrialResult solveStep( ContinuumBase & model,
                       DriverState const & stateOld,
                       Options const & options,
                       PathStep const & step,
                       MaterialFrameUpdate const materialFrameUpdate,
                       StressControlDiagnostics * const diagnostics,
                       localIndex const stepIndex )
{
  RepositorySnapshot snapshot;
  snapshot.capture( model );

  int numUnknowns = 0;
  for( ControlMode const mode : options.controlModes )
  {
    if( mode == ControlMode::Stress )
    {
      ++numUnknowns;
    }
  }

  if( numUnknowns == 0 )
  {
    TrialResult result = evaluateTrial( model, stateOld, snapshot, options, step, {}, materialFrameUpdate, true );
    result.newtonIterations = 0;
    result.residualNorm = 0.0;
    result.converged = true;
    return result;
  }

  std::vector< real64 > unknowns( static_cast< std::size_t >( numUnknowns ), 0.0 );
  for( int c = 0, unknownIndex = 0; c < 6; ++c )
  {
    if( options.controlModes[c] == ControlMode::Stress )
    {
      unknowns[static_cast< std::size_t >( unknownIndex++ )] = stateOld.previousStrainIncrement[c];
    }
  }
  enforceStressControlStrainBounds( options, stateOld, unknowns );
  std::vector< real64 > bestUnknowns = unknowns;
  real64 bestResidualNorm = std::numeric_limits< real64 >::infinity();
  bool haveBest = false;
  std::string stopReason;
  int iteration = 0;

  if( options.stressControlAlgorithm != StressControlAlgorithm::Servo )
  {
    for( iteration = 0; iteration < options.maxNewtonIterations; ++iteration )
    {
    enforceStressControlStrainBounds( options, stateOld, unknowns );
    TrialResult current = evaluateTrial( model, stateOld, snapshot, options, step, unknowns, materialFrameUpdate, false );
    recordStressControlDiagnostic( diagnostics,
                                   model,
                                   stepIndex,
                                   stateOld,
                                   options,
                                   step,
                                   "newton_current",
                                   iteration,
                                   0,
                                   unknowns,
                                   &current,
                                   false,
                                   std::string() );
    std::vector< real64 > residual = stressResidual( options, step, current.state.stress );
    real64 const residualNorm = l2Norm( residual );
    if( isFinite( residualNorm ) && residualNorm < bestResidualNorm )
    {
      bestResidualNorm = residualNorm;
      bestUnknowns = unknowns;
      haveBest = true;
    }
    if( isFinite( residualNorm ) && residualNorm <= options.stressTolerance )
    {
      break;
    }
    if( !isFinite( residualNorm ) || !isFinite( residual ) )
    {
      stopReason = "non-finite stress-control residual";
      if( permissiveStressControlFailurePolicy( options.stressControlFailurePolicy ) )
      {
        break;
      }
      break;
    }

    std::vector< std::vector< real64 > > jacobian( static_cast< std::size_t >( numUnknowns ),
                                                   std::vector< real64 >( static_cast< std::size_t >( numUnknowns ), 0.0 ) );
    for( int j = 0; j < numUnknowns; ++j )
    {
      std::vector< real64 > perturbed = unknowns;
      real64 const dx = std::max( options.finiteDifferenceEpsilon,
                                  options.finiteDifferenceEpsilon * std::max( 1.0, std::abs( perturbed[j] ) ) );
      real64 const baseUnknown = perturbed[static_cast< std::size_t >( j )];
      perturbed[static_cast< std::size_t >( j )] += dx;
      enforceStressControlStrainBounds( options, stateOld, perturbed );
      real64 effectiveDx = perturbed[static_cast< std::size_t >( j )] - baseUnknown;
      if( std::abs( effectiveDx ) <= 16.0 * std::numeric_limits< real64 >::epsilon() * std::max( real64( 1.0 ), std::abs( baseUnknown ) ) )
      {
        perturbed = unknowns;
        perturbed[static_cast< std::size_t >( j )] -= dx;
        enforceStressControlStrainBounds( options, stateOld, perturbed );
        effectiveDx = perturbed[static_cast< std::size_t >( j )] - baseUnknown;
      }
      TrialResult perturbedResult = evaluateTrial( model, stateOld, snapshot, options, step, perturbed, materialFrameUpdate, false );
      recordStressControlDiagnostic( diagnostics,
                                     model,
                                     stepIndex,
                                     stateOld,
                                     options,
                                     step,
                                     "newton_fd",
                                     iteration,
                                     j,
                                     perturbed,
                                     &perturbedResult,
                                     false,
                                     std::string(),
                                     StressControlDiagnosticsLevel::Full );
      std::vector< real64 > perturbedResidual = stressResidual( options, step, perturbedResult.state.stress );
      for( int i = 0; i < numUnknowns; ++i )
      {
        jacobian[static_cast< std::size_t >( i )][static_cast< std::size_t >( j )] =
          ( std::abs( effectiveDx ) > 0.0 ?
            ( perturbedResidual[static_cast< std::size_t >( i )] - residual[static_cast< std::size_t >( i )] ) / effectiveDx :
            real64( 0.0 ) );
      }
    }

    std::vector< real64 > rhs = residual;
    for( real64 & value : rhs )
    {
      value *= -1.0;
    }

    std::vector< real64 > stepDirection;
    try
    {
      if( stressControlUsesRegularization( options.stressControlAlgorithm ) )
      {
        stepDirection = solveRegularizedLeastSquares( jacobian, rhs, options.stressControlRegularization );
      }
      else
      {
        stepDirection = solveLinearSystem( jacobian, rhs );
      }
    }
    catch( std::exception const & error )
    {
      stopReason = error.what();
      if( stressControlUsesRegularization( options.stressControlAlgorithm ) )
      {
        try
        {
          stepDirection = solveRegularizedLeastSquares( jacobian, rhs, 1.0e-8 );
          stopReason = "used regularized stress-control Jacobian after linear solve failure";
        }
        catch( std::exception const & regularizedError )
        {
          stopReason = regularizedError.what();
          if( permissiveStressControlFailurePolicy( options.stressControlFailurePolicy ) )
          {
            break;
          }
          break;
        }
      }
      else
      {
        if( permissiveStressControlFailurePolicy( options.stressControlFailurePolicy ) )
        {
          break;
        }
        break;
      }
    }
    limitStressControlStepDirection( options, stepDirection );

    real64 acceptedNorm = residualNorm;
    std::vector< real64 > acceptedUnknowns = unknowns;
    bool accepted = false;
    real64 damping = 1.0;
    for( int lineSearch = 0; lineSearch < options.maxLineSearchIterations; ++lineSearch )
    {
      std::vector< real64 > trialUnknowns = unknowns;
      for( int j = 0; j < numUnknowns; ++j )
      {
        trialUnknowns[static_cast< std::size_t >( j )] += damping * stepDirection[static_cast< std::size_t >( j )];
      }
      enforceStressControlStrainBounds( options, stateOld, trialUnknowns );
      TrialResult trial = evaluateTrial( model, stateOld, snapshot, options, step, trialUnknowns, materialFrameUpdate, false );
      recordStressControlDiagnostic( diagnostics,
                                     model,
                                     stepIndex,
                                     stateOld,
                                     options,
                                     step,
                                     "line_search",
                                     iteration,
                                     lineSearch,
                                     trialUnknowns,
                                     &trial,
                                     false,
                                     std::string() );
      std::vector< real64 > trialResidual = stressResidual( options, step, trial.state.stress );
      real64 const trialNorm = l2Norm( trialResidual );
      if( isFinite( trialNorm ) && trialNorm < bestResidualNorm )
      {
        bestResidualNorm = trialNorm;
        bestUnknowns = trialUnknowns;
        haveBest = true;
      }
      if( isFinite( trialNorm ) && trialNorm < acceptedNorm )
      {
        acceptedNorm = trialNorm;
        acceptedUnknowns = trialUnknowns;
        accepted = true;
        break;
      }
      damping *= 0.5;
    }
    if( !accepted )
    {
      stopReason = "line search failed to reduce the stress-control residual";
      if( permissiveStressControlFailurePolicy( options.stressControlFailurePolicy ) )
      {
        break;
      }
      break;
    }
    unknowns = acceptedUnknowns;
    enforceStressControlStrainBounds( options, stateOld, unknowns );
    }
  }
  else
  {
    stopReason = "stress-control servo algorithm selected";
  }

  if( numUnknowns == 1 &&
      options.stressControlAlgorithm != StressControlAlgorithm::Servo &&
      ( !haveBest || bestResidualNorm > options.stressTolerance ) )
  {
    scalarStressControlFallback( model,
                                 stateOld,
                                 snapshot,
                                 options,
                                 step,
                                 materialFrameUpdate,
                                 diagnostics,
                                 stepIndex,
                                 bestUnknowns,
                                 bestResidualNorm,
                                 haveBest,
                                 stopReason );
  }

  if( options.stressControlPatternIterations > 0 &&
      ( !haveBest || bestResidualNorm > options.stressTolerance ) )
  {
    patternStressControlFallback( model,
                                  stateOld,
                                  snapshot,
                                  options,
                                  step,
                                  materialFrameUpdate,
                                  diagnostics,
                                  stepIndex,
                                  bestUnknowns,
                                  bestResidualNorm,
                                  haveBest,
                                  stopReason );
  }

  if( stressControlUsesServo( options.stressControlAlgorithm ) &&
      ( !haveBest || bestResidualNorm > options.stressTolerance ) )
  {
    servoStressControlFallback( model,
                                stateOld,
                                snapshot,
                                options,
                                step,
                                materialFrameUpdate,
                                diagnostics,
                                stepIndex,
                                bestUnknowns,
                                bestResidualNorm,
                                haveBest,
                                stopReason );
  }

  std::vector< real64 > finalUnknowns = haveBest ? bestUnknowns : unknowns;
  enforceStressControlStrainBounds( options, stateOld, finalUnknowns );
  TrialResult finalTrial = evaluateTrial( model, stateOld, snapshot, options, step, finalUnknowns, materialFrameUpdate, false );
  recordStressControlDiagnostic( diagnostics,
                                 model,
                                 stepIndex,
                                 stateOld,
                                 options,
                                 step,
                                 "final_best",
                                 iteration,
                                 0,
                                 finalUnknowns,
                                 &finalTrial,
                                 false,
                                 stopReason,
                                 StressControlDiagnosticsLevel::Step );
  finalTrial.newtonIterations = std::min( iteration + 1, options.maxNewtonIterations );
  finalTrial.residualNorm = l2Norm( stressResidual( options, step, finalTrial.state.stress ) );
  finalTrial.converged = finalTrial.residualNorm <= options.stressTolerance;

  if( !finalTrial.converged && options.stressControlFailurePolicy == StressControlFailurePolicy::Error )
  {
    std::ostringstream message;
    message << "stress control did not converge after " << finalTrial.newtonIterations
            << " iterations; residual norm = " << finalTrial.residualNorm;
    if( !stopReason.empty() )
    {
      message << "; reason: " << stopReason;
    }
    throw std::runtime_error( message.str() );
  }

  TrialResult accepted = evaluateTrial( model, stateOld, snapshot, options, step, finalUnknowns, materialFrameUpdate, true );
  recordStressControlDiagnostic( diagnostics,
                                 model,
                                 stepIndex,
                                 stateOld,
                                 options,
                                 step,
                                 "accepted",
                                 iteration,
                                 0,
                                 finalUnknowns,
                                 &accepted,
                                 true,
                                 stopReason,
                                 StressControlDiagnosticsLevel::Step );
  accepted.newtonIterations = finalTrial.newtonIterations;
  accepted.residualNorm = finalTrial.residualNorm;
  accepted.converged = finalTrial.converged;
  if( !accepted.converged )
  {
    std::cerr << "Warning: stress control did not converge after " << accepted.newtonIterations
              << " iterations; residual norm = " << accepted.residualNorm;
    if( !stopReason.empty() )
    {
      std::cerr << "; reason: " << stopReason;
    }
    std::cerr << "\n";
  }
  return accepted;
}

void appendCsvValue( std::ofstream & output, real64 const value )
{
  output << ',' << std::setprecision( 17 ) << value;
}

void writeOutputHeader( std::ofstream & output, ContinuumBase const & model )
{
  output << "step,time,dt";
  std::array< std::string, 6 > const voigt = { "xx", "yy", "zz", "yz", "xz", "xy" };
  for( std::string const & name : voigt ) output << ",eps_" << name;
  for( std::string const & name : voigt ) output << ",stressControlStrainMin_" << name;
  for( std::string const & name : voigt ) output << ",stressControlStrainMax_" << name;
  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) output << ",F_" << i << j;
  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) output << ",L_" << i << j;
  for( std::string const & name : voigt ) output << ",stress_" << name;
  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) output << ",materialDirection_" << i << j;
  output << ",specificInternalEnergy,accumulatedStressPower,stressPower,temperature,temperatureRate,density,jacobian,lengthScale,strengthScale,newtonIterations,stressResidualNorm,converged";

  std::vector< std::string > optionalFields = { "damage",
                                                "basalPlaneDamage",
                                                "comminutionDamage",
                                                "plasticWork",
                                                "basalPlanePlasticWork",
                                                "plasticStrain",
                                                "relaxation",
                                                "effectiveBulkModulus",
                                                "effectiveShearModulus",
                                                "crackTipStressConcentration",
                                                "distanceToCrackTip" };
  for( std::string const & field : optionalFields )
  {
    if( model.hasWrapper( field ) )
    {
      std::vector< real64 > values = getFieldFirstEntryIfPresent( model, field );
      if( values.size() <= 1 )
      {
        output << ',' << field;
      }
      else
      {
        for( std::size_t i = 0; i < values.size(); ++i )
        {
          output << ',' << field << '_' << i;
        }
      }
    }
  }
  output << '\n';
}

void writeOutputRow( std::ofstream & output,
                     ContinuumBase const & model,
                     Options const & options,
                     localIndex const stepIndex,
                     PathStep const & step,
                     TrialResult const & result )
{
  output << stepIndex;
  appendCsvValue( output, result.state.time );
  appendCsvValue( output, step.dt );
  for( real64 const value : result.state.totalStrain ) appendCsvValue( output, value );
  for( real64 const value : options.stressControlMinStrain ) appendCsvValue( output, value );
  for( real64 const value : options.stressControlMaxStrain ) appendCsvValue( output, value );
  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) appendCsvValue( output, result.state.F[i][j] );
  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) appendCsvValue( output, result.L[i][j] );
  for( real64 const value : result.state.stress ) appendCsvValue( output, value );
  for( int i = 0; i < 3; ++i ) for( int j = 0; j < 3; ++j ) appendCsvValue( output, result.state.materialFrame[i][j] );
  appendCsvValue( output, result.state.specificInternalEnergy );
  appendCsvValue( output, result.state.accumulatedStressPower );
  appendCsvValue( output, result.state.lastStressPower );
  appendCsvValue( output, result.state.temperature );
  appendCsvValue( output, result.state.temperatureRate );
  appendCsvValue( output, result.state.density );
  appendCsvValue( output, result.state.jacobian );
  appendCsvValue( output, result.state.lengthScale );
  appendCsvValue( output, result.state.strengthScale );
  output << ',' << result.newtonIterations;
  appendCsvValue( output, result.residualNorm );
  output << ',' << ( result.converged ? 1 : 0 );

  std::vector< std::string > optionalFields = { "damage",
                                                "basalPlaneDamage",
                                                "comminutionDamage",
                                                "plasticWork",
                                                "basalPlanePlasticWork",
                                                "plasticStrain",
                                                "relaxation",
                                                "effectiveBulkModulus",
                                                "effectiveShearModulus",
                                                "crackTipStressConcentration",
                                                "distanceToCrackTip" };
  for( std::string const & field : optionalFields )
  {
    if( model.hasWrapper( field ) )
    {
      std::vector< real64 > values = getFieldFirstEntryIfPresent( model, field );
      if( values.empty() )
      {
        appendCsvValue( output, std::numeric_limits< real64 >::quiet_NaN() );
      }
      else
      {
        for( real64 const value : values ) appendCsvValue( output, value );
      }
    }
  }
  output << '\n';
}

ContinuumBase & createMaterial( Options const & options,
                                Group & /*rootGroup*/,
                                ConstitutiveManager & constitutiveManager,
                                Group & pointGroup )
{
  std::string const materialXml = readFile( options.materialXmlPath );
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( materialXml.c_str() );
  if( !xmlResult )
  {
    std::ostringstream message;
    message << "Material XML parse error: " << xmlResult.description() << " at offset " << xmlResult.offset;
    throw std::runtime_error( message.str() );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.getChild( "Constitutive" );
  constitutiveManager.processInputFileRecursive( xmlDocument, xmlConstitutiveNode );
  constitutiveManager.postInputInitializationRecursive();

  pointGroup.resize( 1 );
  ContinuumBase & model = constitutiveManager.getConstitutiveRelation< ContinuumBase >( options.materialName );
  model.allocateConstitutiveData( pointGroup, 1 );
  model.resize( 1 );
  return model;
}

DriverState initializeState( ContinuumBase & model, Options const & options )
{
  DriverState state;
  state.F = identityMatrix();
  state.materialFrameReference = buildMaterialFrame( options.materialDirectionValues, options.tangentHint );
  state.materialFrame = state.materialFrameReference;
  state.stress = options.initialStress;
  state.temperature = options.initialTemperature;
  state.temperatureRate = options.initialTemperatureRate;
  state.specificInternalEnergy = options.initialSpecificInternalEnergy;
  state.lengthScale = options.initialLengthScale;
  state.strengthScale = options.initialStrengthScale;

  setFieldFirstEntryIfPresent( model, "stress", flatten( state.stress ), true );
  setFieldFirstEntryIfPresent( model, "oldStress", flatten( state.stress ), true );
  setMaterialDirectionFieldIfPresent( model, state.materialFrame );
  setFieldFirstEntryIfPresent( model, "temperature", std::vector< real64 >{ state.temperature } );
  setFieldFirstEntryIfPresent( model, "temperatureRate", std::vector< real64 >{ state.temperatureRate } );
  setFieldFirstEntryIfPresent( model, "jacobian", std::vector< real64 >{ 1.0 } );
  setFieldFirstEntryIfPresent( model, "lengthScale", std::vector< real64 >{ state.lengthScale } );
  setFieldFirstEntryIfPresent( model, "strengthScale", std::vector< real64 >{ state.strengthScale } );
  setFieldFirstEntryIfPresent( model, "specificInternalEnergy", std::vector< real64 >{ state.specificInternalEnergy } );
  setFieldFirstEntryIfPresent( model, "internalEnergy", std::vector< real64 >{ state.specificInternalEnergy } );

  for( InitialField const & field : options.initialFields )
  {
    setFieldFirstEntryIfPresent( model, field.name, field.values, true );
  }

  if( options.initialDensity > 0.0 )
  {
    state.referenceDensity = options.initialDensity;
  }
  else
  {
    std::vector< real64 > densityValues = getFieldFirstEntryIfPresent( model, "density" );
    if( densityValues.empty() || densityValues[0] <= 0.0 )
    {
      throw std::runtime_error( "Material density was not initialized; pass --initial-density" );
    }
    state.referenceDensity = densityValues[0];
  }
  state.density = state.referenceDensity;
  setFieldFirstEntryIfPresent( model, "density", std::vector< real64 >{ state.density } );

  return state;
}

int run( int argc, char ** argv )
{
  Options const options = parseCommandLine( argc, argv );
  std::vector< PathStep > const path = readPathCsv( options );

  conduit::Node node;
  Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );
  Group pointGroup( "materialPoint", &rootGroup );
  ContinuumBase & model = createMaterial( options, rootGroup, constitutiveManager, pointGroup );
  DriverState state = initializeState( model, options );
  MaterialFrameUpdate const frameUpdate = resolveMaterialFrameUpdate( options.materialFrameUpdate, model );

  std::ofstream output( options.outputCsvPath );
  if( !output )
  {
    throw std::runtime_error( "Could not open output CSV: " + options.outputCsvPath );
  }
  writeOutputHeader( output, model );

  StressControlDiagnostics stressControlDiagnostics;
  stressControlDiagnostics.open( options );

  for( localIndex stepIndex = 0; stepIndex < static_cast< localIndex >( path.size() ); ++stepIndex )
  {
    PathStep const & step = path[static_cast< std::size_t >( stepIndex )];
    if( step.dt <= 0.0 )
    {
      throw std::runtime_error( "Non-positive dt at load step " + std::to_string( stepIndex ) );
    }

    TrialResult result;
    try
    {
      result = solveStep( model,
                          state,
                          options,
                          step,
                          frameUpdate,
                          &stressControlDiagnostics,
                          stepIndex + 1 );
    }
    catch( std::exception const & error )
    {
      if( permissiveStressControlFailurePolicy( options.stressControlFailurePolicy ) && hasStressControl( options ) )
      {
        std::cerr << "Warning: stopping material-point driver with partial output at load step "
                  << stepIndex + 1 << " after stress-control failure: " << error.what() << "\n";
        output.flush();
        return 0;
      }
      throw;
    }

    result.state.time = state.time + step.dt;
    writeOutputRow( output, model, options, stepIndex + 1, step, result );
    output.flush();
    state = result.state;

    if( !result.converged &&
        options.stressControlFailurePolicy == StressControlFailurePolicy::Stop &&
        hasStressControl( options ) )
    {
      std::cerr << "Warning: stopping material-point driver with partial output at load step "
                << stepIndex + 1 << " because stress control did not converge; residual norm = "
                << result.residualNorm << "\n";
      return 0;
    }
  }

  return 0;
}

} // namespace materialPointDriver
} // namespace constitutive
} // namespace geos

int main( int argc, char ** argv )
{
  try
  {
    return geos::constitutive::materialPointDriver::run( argc, argv );
  }
  catch( std::exception const & error )
  {
    std::cerr << "geos_mpm_material_point_driver error: " << error.what() << '\n';
    return 1;
  }
}
