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
  buildTrialKinematics( options, step, unknowns, result.strainIncrement, result.L );

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
  }

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

TrialResult solveStep( ContinuumBase & model,
                       DriverState const & stateOld,
                       Options const & options,
                       PathStep const & step,
                       MaterialFrameUpdate const materialFrameUpdate )
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
  TrialResult current;
  real64 residualNorm = std::numeric_limits< real64 >::infinity();
  int iteration = 0;

  for( iteration = 0; iteration < options.maxNewtonIterations; ++iteration )
  {
    current = evaluateTrial( model, stateOld, snapshot, options, step, unknowns, materialFrameUpdate, false );
    std::vector< real64 > residual = stressResidual( options, step, current.state.stress );
    residualNorm = l2Norm( residual );
    if( residualNorm <= options.stressTolerance )
    {
      break;
    }

    std::vector< std::vector< real64 > > jacobian( static_cast< std::size_t >( numUnknowns ),
                                                   std::vector< real64 >( static_cast< std::size_t >( numUnknowns ), 0.0 ) );
    for( int j = 0; j < numUnknowns; ++j )
    {
      std::vector< real64 > perturbed = unknowns;
      real64 const dx = std::max( options.finiteDifferenceEpsilon,
                                  options.finiteDifferenceEpsilon * std::max( 1.0, std::abs( perturbed[j] ) ) );
      perturbed[j] += dx;
      TrialResult perturbedResult = evaluateTrial( model, stateOld, snapshot, options, step, perturbed, materialFrameUpdate, false );
      std::vector< real64 > perturbedResidual = stressResidual( options, step, perturbedResult.state.stress );
      for( int i = 0; i < numUnknowns; ++i )
      {
        jacobian[static_cast< std::size_t >( i )][static_cast< std::size_t >( j )] =
          ( perturbedResidual[static_cast< std::size_t >( i )] - residual[static_cast< std::size_t >( i )] ) / dx;
      }
    }

    std::vector< real64 > rhs = residual;
    for( real64 & value : rhs )
    {
      value *= -1.0;
    }
    std::vector< real64 > stepDirection = solveLinearSystem( jacobian, rhs );

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
      TrialResult trial = evaluateTrial( model, stateOld, snapshot, options, step, trialUnknowns, materialFrameUpdate, false );
      std::vector< real64 > trialResidual = stressResidual( options, step, trial.state.stress );
      real64 const trialNorm = l2Norm( trialResidual );
      if( trialNorm < acceptedNorm || lineSearch == options.maxLineSearchIterations - 1 )
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
      throw std::runtime_error( "Line search failed in material-point stress control" );
    }
    unknowns = acceptedUnknowns;
  }

  TrialResult accepted = evaluateTrial( model, stateOld, snapshot, options, step, unknowns, materialFrameUpdate, true );
  accepted.newtonIterations = iteration + 1;
  accepted.residualNorm = l2Norm( stressResidual( options, step, accepted.state.stress ) );
  accepted.converged = accepted.residualNorm <= options.stressTolerance;
  if( !accepted.converged )
  {
    std::cerr << "Warning: stress control did not converge after " << options.maxNewtonIterations
              << " iterations; residual norm = " << accepted.residualNorm << "\n";
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
                     localIndex const stepIndex,
                     PathStep const & step,
                     TrialResult const & result )
{
  output << stepIndex;
  appendCsvValue( output, result.state.time );
  appendCsvValue( output, step.dt );
  for( real64 const value : result.state.totalStrain ) appendCsvValue( output, value );
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

  for( localIndex stepIndex = 0; stepIndex < static_cast< localIndex >( path.size() ); ++stepIndex )
  {
    PathStep const & step = path[static_cast< std::size_t >( stepIndex )];
    if( step.dt <= 0.0 )
    {
      throw std::runtime_error( "Non-positive dt at load step " + std::to_string( stepIndex ) );
    }
    TrialResult result = solveStep( model, state, options, step, frameUpdate );
    result.state.time = state.time + step.dt;
    writeOutputRow( output, model, stepIndex + 1, step, result );
    state = result.state;
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
