#include "linearAlgebra/interfaces/hypre/hypredrive.hpp"

#include "common/GeosxConfig.hpp"
#include "common/MpiWrapper.hpp"
#include "common/format/Format.hpp"
#include "common/format/StringUtilities.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"
#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/AugmentedLagrangianContactMechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/HybridSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/Hydrofracture.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ImmiscibleMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/LagrangianContactMechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/LagrangianContactMechanicsBubbleStabilization.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/MultiphasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ReactiveCompositionalMultiphaseOBL.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsConformingFractures.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsConformingFracturesALM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsConformingFracturesALMReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsEmbeddedFractures.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhaseReservoirHybridFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SolidMechanicsEmbeddedFractures.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalCompositionalMultiphaseReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalMultiphasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalSinglePhasePoromechanics.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalSinglePhasePoromechanicsReservoirFVM.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/ThermalSinglePhaseReservoirFVM.hpp"

#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cfenv>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace geos
{

namespace
{

void checkHypredriveCall( uint32_t const errorCode,
                          char const * const call,
                          std::string const & context = {} )
{
  if( errorCode != 0 )
  {
    HYPREDRV_ErrorCodeDescribe( errorCode );
    std::fflush( stderr );
    HYPREDRV_ErrorCodeClear();

    std::string const contextMessage = context.empty()
                                             ? std::string{}
                                             : GEOS_FMT( "\nContext: {}", context );
    GEOS_ERROR( GEOS_FMT( "Error in call to {} (HYPREDRV error code: 0x{:08x}, decimal: {}).{}\n"
                          "HypreDrive error details are written to stderr; use '2>&1' to merge them into stdout. "
                          "For additional library tracing, set HYPREDRV_LOG_LEVEL=3 and "
                          "HYPREDRV_LOG_STREAM=stdout.",
                          call,
                          errorCode,
                          errorCode,
                          contextMessage ) );
  }
}

bool envFlagEnabled( char const * const name )
{
  char const * const value = std::getenv( name );
  if( value == nullptr )
  {
    return false;
  }

  std::string normalized( value );
  std::transform( normalized.begin(), normalized.end(), normalized.begin(),
                  []( unsigned char const c ) { return static_cast< char >( std::tolower( c ) ); } );

  return !( normalized.empty() || normalized == "0" || normalized == "false" || normalized == "off" || normalized == "no" );
}

std::optional< std::string > readTextFileContents( std::string const & path )
{
  std::ifstream input( path );
  if( !input )
  {
    return std::nullopt;
  }

  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

std::string makeInputArgsParseTargetSignature( hypre::hypredrive::InputArgsParseTarget const & target )
{
  return GEOS_FMT( "{}\n{}",
                   target.source == hypre::hypredrive::InputSource::authoritativeFile ? "authoritativeFile" : "generatedFallback",
                   target.argument );
}

std::set< std::string > & loggedInputArgsParseTargets()
{
  static std::set< std::string > loggedTargets;
  return loggedTargets;
}

HYPRE_Matrix toHypreMatrix( HypreMatrix::HYPRE_IJMatrix const matrix )
{
  return reinterpret_cast< HYPRE_Matrix >( matrix );
}

void appendLine( std::ostringstream & stream,
                 integer const indentLevel,
                 std::string const & line )
{
  stream << std::string( static_cast< size_t >( indentLevel ) * 2, ' ' ) << line << '\n';
}

void appendIlUDisableRcm( std::ostringstream & stream,
                          integer const indentLevel )
{
  // hypre_ILUCreate and HYPRE_BoomerAMGCreate default local reordering to RCM (1).
  // GEOS never uses that reordering on ILU.
  appendLine( stream, indentLevel, "reordering: 0" );
}

char const * getSolverName( LinearSolverParameters::SolverType const solverType )
{
  using SolverType = LinearSolverParameters::SolverType;

  switch( solverType )
  {
    case SolverType::cg:
      return "pcg";
    case SolverType::gmres:
      return "gmres";
    case SolverType::fgmres:
      return "fgmres";
    case SolverType::bicgstab:
      return "bicgstab";
    case SolverType::direct:
    case SolverType::richardson:
    case SolverType::preconditioner:
      return nullptr;
  }

  return nullptr;
}

bool isSupportedGeneratedAMGSmoother( LinearSolverParameters::AMG::SmootherType const smootherType )
{
  using SmootherType = LinearSolverParameters::AMG::SmootherType;

  switch( smootherType )
  {
    case SmootherType::default_:
    case SmootherType::jacobi:
    case SmootherType::l1jacobi:
    case SmootherType::fgs:
    case SmootherType::bgs:
    case SmootherType::sgs:
    case SmootherType::l1sgs:
    case SmootherType::chebyshev:
    case SmootherType::iluk:
    case SmootherType::ilut:
      return true;
    case SmootherType::ick:
    case SmootherType::ict:
      return false;
  }

  return false;
}

char const * getGeneratedAMGRelaxationName( LinearSolverParameters::AMG::SmootherType const smootherType )
{
  using SmootherType = LinearSolverParameters::AMG::SmootherType;

  switch( smootherType )
  {
    case SmootherType::default_:
    case SmootherType::iluk:
    case SmootherType::ilut:
    case SmootherType::ick:
    case SmootherType::ict:
      return nullptr;
    case SmootherType::jacobi:
      return "jacobi";
    case SmootherType::l1jacobi:
      return "l1-jacobi";
    case SmootherType::fgs:
      return "forward-hl1gs";
    case SmootherType::bgs:
      return "backward-hl1gs";
    case SmootherType::sgs:
      return "hsgs";
    case SmootherType::l1sgs:
      return "l1-hsgs";
    case SmootherType::chebyshev:
      return "chebyshev";
  }

  return nullptr;
}

char const * getGeneratedAMGCoarseName( LinearSolverParameters::AMG::CoarseType const coarseType )
{
  using CoarseType = LinearSolverParameters::AMG::CoarseType;

  switch( coarseType )
  {
    case CoarseType::default_:
      return nullptr;
    case CoarseType::jacobi:
      return "jacobi";
    case CoarseType::l1jacobi:
      return "l1-jacobi";
    case CoarseType::fgs:
      return "forward-hl1gs";
    case CoarseType::bgs:
      return "backward-hl1gs";
    case CoarseType::sgs:
      return "hsgs";
    case CoarseType::l1sgs:
      return "l1-hsgs";
    case CoarseType::chebyshev:
      return "chebyshev";
    case CoarseType::direct:
      return "ge";
    case CoarseType::gsElimWPivoting:
      return "lu_piv";
    case CoarseType::gsElimWInverse:
      return "lu_inv";
  }

  return nullptr;
}

bool supportsGeneratedPreconditioner( LinearSolverParameters::PreconditionerType const preconditionerType )
{
  using PreconditionerType = LinearSolverParameters::PreconditionerType;

  switch( preconditionerType )
  {
    case PreconditionerType::none:
    case PreconditionerType::amg:
    case PreconditionerType::mgr:
    case PreconditionerType::iluk:
    case PreconditionerType::ilut:
      return true;
    case PreconditionerType::jacobi:
    case PreconditionerType::l1jacobi:
    case PreconditionerType::fgs:
    case PreconditionerType::sgs:
    case PreconditionerType::l1sgs:
    case PreconditionerType::chebyshev:
    case PreconditionerType::ick:
    case PreconditionerType::ict:
    case PreconditionerType::block:
    case PreconditionerType::direct:
    case PreconditionerType::bgs:
    case PreconditionerType::multiscale:
      return false;
  }

  return false;
}

HYPRE_Int getAMGPrintLevel( integer const logLevel )
{
  return ( logLevel == 2 || logLevel >= 4 ) ? 1 : 0;
}

HYPRE_Int getMGRPrintLevel( integer const logLevel )
{
  HYPRE_Int hypreLogLevel = LvArray::math::max( LvArray::integerConversion< HYPRE_Int >( logLevel - 1 ),
                                                LvArray::integerConversion< HYPRE_Int >( 0 ) );
  hypreLogLevel &= ~0x2;
  return hypreLogLevel;
}

char const * getGeneratedMGRFRelaxationName( hypre::MGRFRelaxationType const relaxationType )
{
  using RelaxationType = hypre::MGRFRelaxationType;

  switch( relaxationType )
  {
    case RelaxationType::none:
      return "none";
    case RelaxationType::singleVCycleSmoother:
      return "v(1,0)";
    case RelaxationType::amgVCycle:
      return nullptr;
    case RelaxationType::jacobi:
      return "jacobi";
    case RelaxationType::l1jacobi:
      return "l1-jacobi";
    case RelaxationType::ilu:
      return "ilu";
    case RelaxationType::gsElim:
      return "ge";
    case RelaxationType::gsElimWPivoting:
      return "ge-piv";
    case RelaxationType::gsElimWInverse:
      return "ge-inv";
    case RelaxationType::forwardHybridGaussSeidel:
    case RelaxationType::backwardHybridGaussSeidel:
    case RelaxationType::hybridSymmetricGaussSeidel:
    case RelaxationType::l1hybridSymmetricGaussSeidel:
    case RelaxationType::l1forwardGaussSeidel:
    case RelaxationType::l1backwardGaussSeidel:
      return nullptr;
  }

  return nullptr;
}

char const * getGeneratedMGRGlobalSmootherName( hypre::MGRGlobalSmootherType const smootherType )
{
  using SmootherType = hypre::MGRGlobalSmootherType;

  switch( smootherType )
  {
    case SmootherType::none:
      return "none";
    case SmootherType::blockJacobi:
      return "blk-jacobi";
    case SmootherType::blockGaussSeidel:
      return "blk-gs";
    case SmootherType::jacobi:
      // hypredrive follows hypre's MGR global smoother table, where code 2 is exposed as mixed-gs.
      return "mixed-gs";
    case SmootherType::ilu0:
      return "ilu";
  }

  return nullptr;
}

std::string sanitizeLabelToken( std::string const & name )
{
  std::string sanitized;
  sanitized.reserve( name.size() );

  bool lastWasUnderscore = false;
  for( unsigned char const c : name )
  {
    if( std::isalnum( c ) )
    {
      sanitized.push_back( static_cast< char >( c ) );
      lastWasUnderscore = false;
    }
    else if( !lastWasUnderscore )
    {
      sanitized.push_back( '_' );
      lastWasUnderscore = true;
    }
  }

  while( !sanitized.empty() && sanitized.front() == '_' )
  {
    sanitized.erase( sanitized.begin() );
  }
  while( !sanitized.empty() && sanitized.back() == '_' )
  {
    sanitized.pop_back();
  }

  return sanitized;
}

std::string normalizeLabelToken( std::string const & name )
{
  std::string normalized;
  normalized.reserve( name.size() );
  for( unsigned char const c : name )
  {
    normalized.push_back( static_cast< char >( std::tolower( c ) ) );
  }
  return normalized;
}

bool strategyUsesCompositionalSemanticLabels( LinearSolverParameters::MGR::StrategyType const strategy )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  switch( strategy )
  {
    case StrategyType::compositionalMultiphaseFVM:
    case StrategyType::compositionalMultiphaseHybridFVM:
    case StrategyType::compositionalMultiphaseReservoirFVM:
    case StrategyType::compositionalMultiphaseReservoirHybridFVM:
    case StrategyType::immiscibleMultiphaseFVM:
    case StrategyType::reactiveCompositionalMultiphaseOBL:
    case StrategyType::thermalCompositionalMultiphaseFVM:
    case StrategyType::thermalCompositionalMultiphaseReservoirFVM:
    case StrategyType::multiphasePoromechanics:
    case StrategyType::multiphasePoromechanicsReservoirFVM:
    case StrategyType::thermalMultiphasePoromechanics:
      return true;
    case StrategyType::invalid:
    case StrategyType::singlePhaseReservoirFVM:
    case StrategyType::thermalSinglePhaseReservoirFVM:
    case StrategyType::singlePhaseHybridFVM:
    case StrategyType::singlePhaseReservoirHybridFVM:
    case StrategyType::singlePhasePoromechanics:
    case StrategyType::thermalSinglePhasePoromechanics:
    case StrategyType::hybridSinglePhasePoromechanics:
    case StrategyType::singlePhasePoromechanicsEmbeddedFractures:
    case StrategyType::singlePhasePoromechanicsConformingFractures:
    case StrategyType::singlePhasePoromechanicsConformingFracturesALM:
    case StrategyType::singlePhasePoromechanicsConformingFracturesALMReservoirFVM:
    case StrategyType::singlePhasePoromechanicsReservoirFVM:
    case StrategyType::thermalSinglePhasePoromechanicsReservoirFVM:
    case StrategyType::hydrofracture:
    case StrategyType::lagrangianContactMechanics:
    case StrategyType::augmentedLagrangianContactMechanics:
    case StrategyType::lagrangianContactMechanicsBubbleStab:
    case StrategyType::solidMechanicsEmbeddedFractures:
      return false;
  }

  return false;
}

bool strategyUsesTemperatureLabel( LinearSolverParameters::MGR::StrategyType const strategy )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  switch( strategy )
  {
    case StrategyType::thermalCompositionalMultiphaseFVM:
    case StrategyType::thermalCompositionalMultiphaseReservoirFVM:
    case StrategyType::thermalMultiphasePoromechanics:
      return true;
    default:
      return false;
  }
}

bool strategyUsesCompositionalWellSemanticLabels( LinearSolverParameters::MGR::StrategyType const strategy )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  switch( strategy )
  {
    case StrategyType::compositionalMultiphaseReservoirFVM:
    case StrategyType::compositionalMultiphaseReservoirHybridFVM:
    case StrategyType::thermalCompositionalMultiphaseReservoirFVM:
    case StrategyType::multiphasePoromechanicsReservoirFVM:
      return true;
    default:
      return false;
  }
}

bool strategyUsesWellTemperatureLabel( LinearSolverParameters::MGR::StrategyType const strategy )
{
  return strategy == LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseReservoirFVM;
}

stdVector< string > buildCompositionalFieldLabelNames( LinearSolverParameters::MGR::StrategyType const strategy,
                                                       int const numComponents )
{
  stdVector< string > labelNames;
  labelNames.reserve( std::max( numComponents, 0 ) );

  if( numComponents <= 0 )
  {
    return labelNames;
  }

  bool const hasTemperature = strategyUsesTemperatureLabel( strategy );
  labelNames.emplace_back( "pressure" );

  int const densityCount = hasTemperature
                           ? std::max( numComponents - 2, 0 )
                           : std::max( numComponents - 1, 0 );
  for( int densityIndex = 0; densityIndex < densityCount; ++densityIndex )
  {
    labelNames.emplace_back( GEOS_FMT( "density_{}", densityIndex ) );
  }

  if( hasTemperature && numComponents > 1 )
  {
    labelNames.emplace_back( "temperature" );
  }

  return labelNames;
}

stdVector< string > buildCompositionalWellFieldLabelNames( LinearSolverParameters::MGR::StrategyType const strategy,
                                                           int const numComponents )
{
  stdVector< string > labelNames;
  labelNames.reserve( std::max( numComponents, 0 ) );

  if( numComponents <= 0 )
  {
    return labelNames;
  }

  bool const hasTemperature = strategyUsesWellTemperatureLabel( strategy );
  labelNames.emplace_back( "wellPressure" );

  int const densityCount = hasTemperature
                           ? std::max( numComponents - 3, 0 )
                           : std::max( numComponents - 2, 0 );
  for( int densityIndex = 0; densityIndex < densityCount; ++densityIndex )
  {
    labelNames.emplace_back( GEOS_FMT( "wellDensity_{}", densityIndex ) );
  }

  if( numComponents > 1 )
  {
    labelNames.emplace_back( "wellRate" );
  }

  if( hasTemperature && numComponents > 2 )
  {
    labelNames.emplace_back( "wellTemperature" );
  }

  return labelNames;
}

stdVector< string > buildSinglePhaseFieldLabelNames( int const numComponents )
{
  stdVector< string > labelNames;
  labelNames.reserve( std::max( numComponents, 0 ) );

  if( numComponents <= 0 )
  {
    return labelNames;
  }

  labelNames.emplace_back( "pressure" );

  if( numComponents > 1 )
  {
    labelNames.emplace_back( "temperature" );
  }

  for( int component = LvArray::integerConversion< int >( labelNames.size() ); component < numComponents; ++component )
  {
    labelNames.emplace_back( GEOS_FMT( "singlePhaseVariables_{}", component ) );
  }

  return labelNames;
}

stdVector< string > buildSinglePhaseWellFieldLabelNames( int const numComponents )
{
  stdVector< string > labelNames;
  labelNames.reserve( std::max( numComponents, 0 ) );

  if( numComponents <= 0 )
  {
    return labelNames;
  }

  labelNames.emplace_back( "wellPressure" );

  if( numComponents > 1 )
  {
    labelNames.emplace_back( "wellRate" );
  }

  if( numComponents > 2 )
  {
    labelNames.emplace_back( "wellTemperature" );
  }

  for( int component = LvArray::integerConversion< int >( labelNames.size() ); component < numComponents; ++component )
  {
    labelNames.emplace_back( GEOS_FMT( "singlePhaseWellVars_{}", component ) );
  }

  return labelNames;
}

stdVector< string > buildFieldComponentLabelNames( LinearSolverParameters::MGR::StrategyType const strategy,
                                                   std::string const & baseName,
                                                   int const numComponents )
{
  stdVector< string > labelNames;
  labelNames.reserve( std::max( numComponents, 0 ) );

  if( numComponents <= 0 )
  {
    return labelNames;
  }

  std::string const normalizedBaseName = normalizeLabelToken( baseName );

  if( normalizedBaseName == "singlephasevariables" )
  {
    return buildSinglePhaseFieldLabelNames( numComponents );
  }

  if( normalizedBaseName == "singlephasewellvars" )
  {
    return buildSinglePhaseWellFieldLabelNames( numComponents );
  }

  if( normalizedBaseName == "compositionalvariables" && strategyUsesCompositionalSemanticLabels( strategy ) )
  {
    return buildCompositionalFieldLabelNames( strategy, numComponents );
  }

  if( normalizedBaseName == "compositionalwellvars" && strategyUsesCompositionalWellSemanticLabels( strategy ) )
  {
    return buildCompositionalWellFieldLabelNames( strategy, numComponents );
  }

  for( int component = 0; component < numComponents; ++component )
  {
    labelNames.emplace_back( numComponents == 1
                             ? baseName
                             : GEOS_FMT( "{}_{}", baseName, component ) );
  }

  return labelNames;
}

stdVector< string > buildDofLabelNames( LinearSolverParameters::MGR::StrategyType const strategy,
                                        stdVector< string > const & fieldNames,
                                        arrayView1d< int const > const & numComponentsPerField )
{
  stdVector< string > labels;
  labels.reserve( std::accumulate( numComponentsPerField.begin(), numComponentsPerField.end(), 0 ) );

  std::set< string > usedNames;
  for( localIndex fieldIndex = 0; fieldIndex < numComponentsPerField.size(); ++fieldIndex )
  {
    std::string baseName = static_cast< size_t >( fieldIndex ) < fieldNames.size()
                           ? sanitizeLabelToken( fieldNames[fieldIndex] )
                           : "";
    if( baseName.empty() )
    {
      baseName = GEOS_FMT( "field_{}", fieldIndex );
    }

    for( std::string labelName : buildFieldComponentLabelNames( strategy,
                                                                baseName,
                                                                numComponentsPerField[fieldIndex] ) )
    {
      if( labelName.empty() || std::isdigit( static_cast< unsigned char >( labelName.front() ) ) )
      {
        labelName = GEOS_FMT( "dof_{}", labels.size() );
      }

      if( !usedNames.insert( labelName ).second )
      {
        labelName = GEOS_FMT( "{}_{}", labelName, labels.size() );
        usedNames.insert( labelName );
      }

      labels.emplace_back( std::move( labelName ) );
    }
  }

  return labels;
}

std::string joinLabelNames( stdVector< HYPRE_Int > const & values,
                            stdVector< string > const & labelNames )
{
  std::ostringstream stream;
  for( size_t i = 0; i < values.size(); ++i )
  {
    HYPRE_Int const label = values[i];
    GEOS_ERROR_IF( label < 0 || label >= LvArray::integerConversion< HYPRE_Int >( labelNames.size() ),
                   GEOS_FMT( "Invalid dof label {} while generating hypredrive MGR YAML", label ) );
    if( i > 0 )
    {
      stream << ", ";
    }
    stream << labelNames[label];
  }
  return stream.str();
}

void destroyWrapper( HyprePrecWrapper & wrapper )
{
  if( wrapper.ptr != nullptr && wrapper.destroy != nullptr )
  {
    GEOS_LAI_CHECK_ERROR( wrapper.destroy( wrapper.ptr ) );
    wrapper.ptr = nullptr;
  }
}

std::string buildDofLabelsYaml( stdVector< string > const & labelNames )
{
  if( labelNames.empty() )
  {
    return {};
  }

  std::ostringstream stream;
  appendLine( stream, 0, "linear_system:" );
  appendLine( stream, 1, "dof_labels:" );
  for( size_t label = 0; label < labelNames.size(); ++label )
  {
    appendLine( stream, 2, GEOS_FMT( "{}: {}", labelNames[label], label ) );
  }

  return stream.str();
}

void appendSolverYaml( std::ostringstream & stream,
                       LinearSolverParameters const & params )
{
  char const * const solverName = getSolverName( params.solverType );
  GEOS_ERROR_IF( solverName == nullptr,
                 GEOS_FMT( "Unsupported generated hypredrive solver {}", params.solverType ) );

  appendLine( stream, 0, "solver:" );
  appendLine( stream, 1, GEOS_FMT( "{}:", solverName ) );
  appendLine( stream, 2, GEOS_FMT( "max_iter: {}", params.krylov.maxIterations ) );
  appendLine( stream, 2, GEOS_FMT( "relative_tol: {}", params.krylov.relTolerance ) );
  appendLine( stream, 2, "print_level: 0" );

  if( params.solverType == LinearSolverParameters::SolverType::gmres ||
      params.solverType == LinearSolverParameters::SolverType::fgmres )
  {
    appendLine( stream, 2, GEOS_FMT( "krylov_dim: {}", params.krylov.maxRestart ) );
  }
}

bool buildAMGPreconditionerYaml( LinearSolverParameters const & params,
                                 std::string & yaml )
{
  using AMG = LinearSolverParameters::AMG;

  if( params.amg.nullSpaceType == AMG::NullSpaceType::rigidBodyModes )
  {
    return false;
  }

  if( !isSupportedGeneratedAMGSmoother( params.amg.smootherType ) )
  {
    return false;
  }

  std::ostringstream stream;
  appendLine( stream, 0, "preconditioner:" );
  appendLine( stream, 1, "amg:" );
  appendLine( stream, 2, "tolerance: 0.0" );
  appendLine( stream, 2, GEOS_FMT( "max_iter: {}", params.amg.numCycles ) );
  appendLine( stream, 2, GEOS_FMT( "print_level: {}", getAMGPrintLevel( params.logLevel ) ) );

  appendLine( stream, 2, "interpolation:" );
  appendLine( stream, 3, GEOS_FMT( "prolongation_type: {}", hypre::getAMGInterpolationType( params.amg.interpolationType ) ) );
  appendLine( stream, 3, GEOS_FMT( "max_nnz_row: {}", params.amg.interpolationMaxNonZeros ) );

  appendLine( stream, 2, "coarsening:" );
  appendLine( stream, 3, GEOS_FMT( "type: {}", hypre::getAMGCoarseningType( params.amg.coarseningType ) ) );
  appendLine( stream, 3, GEOS_FMT( "max_levels: {}", params.amg.maxLevels ) );
  appendLine( stream, 3, GEOS_FMT( "max_coarse_size: {}", params.amg.maxCoarseSize ) );
  appendLine( stream, 3, GEOS_FMT( "num_functions: {}", params.amg.numFunctions ) );
  appendLine( stream, 3, GEOS_FMT( "filter_functions: {}", params.amg.separateComponents ) );
  if( params.amg.threshold > 0.0 )
  {
    appendLine( stream, 3, GEOS_FMT( "strong_th: {}", params.amg.threshold ) );
  }

  if( params.amg.aggressiveNumLevels > 0 )
  {
    appendLine( stream, 2, "aggressive:" );
    appendLine( stream, 3, GEOS_FMT( "num_levels: {}", params.amg.aggressiveNumLevels ) );
    appendLine( stream, 3, GEOS_FMT( "num_paths: {}", params.amg.aggressiveNumPaths ) );
    appendLine( stream, 3, GEOS_FMT( "prolongation_type: {}", hypre::getAMGAggressiveInterpolationType( params.amg.aggressiveInterpType ) ) );
    appendLine( stream, 3, GEOS_FMT( "max_nnz_row: {}", params.amg.aggressiveInterpMaxNonZeros ) );
  }

  if( params.amg.smootherType == AMG::SmootherType::iluk ||
      params.amg.smootherType == AMG::SmootherType::ilut )
  {
    appendLine( stream, 2, "smoother:" );
    appendLine( stream, 3, GEOS_FMT( "num_levels: {}", params.amg.maxLevels ) );
    appendLine( stream, 3, GEOS_FMT( "num_sweeps: {}", params.amg.numSweeps ) );
    appendLine( stream, 3, "ilu:" );
    appendLine( stream, 4, GEOS_FMT( "type: {}", hypre::getILUType( params.amg.smootherType ) ) );
    appendLine( stream, 4, GEOS_FMT( "fill_level: {}", params.ifact.fill ) );
    appendIlUDisableRcm( stream, 4 );
    if( params.amg.smootherType == AMG::SmootherType::ilut )
    {
      appendLine( stream, 4, GEOS_FMT( "droptol: {}", params.ifact.threshold ) );
    }
  }
  else
  {
    char const * const relaxType = getGeneratedAMGRelaxationName( params.amg.smootherType );

    appendLine( stream, 2, "relaxation:" );
    appendLine( stream, 3, GEOS_FMT( "weight: {}", params.amg.relaxWeight ) );

    if( relaxType != nullptr )
    {
      appendLine( stream, 3, GEOS_FMT( "down_type: {}", relaxType ) );
      appendLine( stream, 3, GEOS_FMT( "up_type: {}", relaxType ) );
    }

    char const * const coarseType = getGeneratedAMGCoarseName( params.amg.coarseType );
    if( coarseType != nullptr )
    {
      appendLine( stream, 3, GEOS_FMT( "coarse_type: {}", coarseType ) );
    }
    else if( relaxType != nullptr )
    {
      appendLine( stream, 3, GEOS_FMT( "coarse_type: {}", relaxType ) );
    }

    switch( params.amg.preOrPostSmoothing )
    {
      case AMG::PreOrPost::both:
      {
        appendLine( stream, 3, GEOS_FMT( "num_sweeps: {}", params.amg.numSweeps ) );
        break;
      }
      case AMG::PreOrPost::pre:
      {
        appendLine( stream, 3, GEOS_FMT( "down_sweeps: {}", params.amg.numSweeps ) );
        appendLine( stream, 3, "up_sweeps: 0" );
        break;
      }
      case AMG::PreOrPost::post:
      {
        appendLine( stream, 3, "down_sweeps: 0" );
        appendLine( stream, 3, GEOS_FMT( "up_sweeps: {}", params.amg.numSweeps ) );
        break;
      }
    }

    if( params.amg.smootherType == AMG::SmootherType::chebyshev &&
        params.amg.numSweeps > 0 && params.amg.numSweeps < 5 )
    {
      appendLine( stream, 3, "chebyshev:" );
      appendLine( stream, 4, GEOS_FMT( "order: {}", params.amg.numSweeps ) );
    }
  }

  yaml = stream.str();
  return true;
}

bool buildILUPreconditionerYaml( LinearSolverParameters const & params,
                                 std::string & yaml )
{
  std::ostringstream stream;
  appendLine( stream, 0, "preconditioner:" );
  appendLine( stream, 1, "ilu:" );
  appendLine( stream, 2, "tolerance: 0.0" );
  appendLine( stream, 2, "max_iter: 1" );
  appendLine( stream, 2, "print_level: 0" );
  appendLine( stream, 2, GEOS_FMT( "type: {}", hypre::getILUType( params.preconditionerType ) ) );

  if( params.ifact.fill >= 0 )
  {
    appendLine( stream, 2, GEOS_FMT( "fill_level: {}", params.ifact.fill ) );
  }
  if( params.preconditionerType == LinearSolverParameters::PreconditionerType::ilut &&
      params.ifact.threshold >= 0.0 )
  {
    appendLine( stream, 2, GEOS_FMT( "droptol: {}", params.ifact.threshold ) );
  }
  // HYPRE_ILUCreate defaults to RCM reordering. Keep the legacy path's
  // behavior for scalar systems while retaining the mechanics safeguard.
  appendLine( stream, 2, GEOS_FMT( "reordering: {}", params.dofsPerNode > 1 ? 0 : 1 ) );

  yaml = stream.str();
  return true;
}

enum class AMGFlavor
{
  pressure,
  pressureTemperature,
  displacementFiltered,
  displacement,
  almDisplacement
};

struct LevelAMGBlock
{
  LevelAMGBlock() = default;

  LevelAMGBlock( HYPRE_Int const inputLevel,
                 AMGFlavor const inputFlavor )
    : level( inputLevel ),
    flavor( inputFlavor )
  {}

  // GCC 12 emits a false positive from std::vector's memmove fast path when this
  // type stays trivially relocatable under -Wnonnull.
  LevelAMGBlock( LevelAMGBlock const & other )
    : level( other.level ),
    flavor( other.flavor )
  {}

  HYPRE_Int level;
  AMGFlavor flavor;
};

struct MGRSpecialization
{
  AMGFlavor coarseFlavor = AMGFlavor::pressure;
  stdVector< LevelAMGBlock > fRelaxAMGLevels;
  HYPRE_Int pmax = 0;
  HYPRE_Int coarseMinCoarseSize = -1;
  char const * cycle = nullptr;
};

MGRSpecialization getSpecialization( LinearSolverParameters::MGR::StrategyType const strategy )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  switch( strategy )
  {
    case StrategyType::thermalSinglePhaseReservoirFVM:
    case StrategyType::thermalCompositionalMultiphaseFVM:
    case StrategyType::thermalCompositionalMultiphaseReservoirFVM:
    {
      return { AMGFlavor::pressureTemperature, {}, 0, -1 };
    }
    case StrategyType::lagrangianContactMechanics:
    case StrategyType::augmentedLagrangianContactMechanics:
    case StrategyType::lagrangianContactMechanicsBubbleStab:
    case StrategyType::solidMechanicsEmbeddedFractures:
    {
      return { AMGFlavor::displacementFiltered, {}, 0, -1 };
    }
    case StrategyType::singlePhasePoromechanics:
    case StrategyType::thermalSinglePhasePoromechanics:
    case StrategyType::hybridSinglePhasePoromechanics:
    case StrategyType::singlePhasePoromechanicsReservoirFVM:
    case StrategyType::thermalSinglePhasePoromechanicsReservoirFVM:
    case StrategyType::multiphasePoromechanics:
    case StrategyType::multiphasePoromechanicsReservoirFVM:
    case StrategyType::thermalMultiphasePoromechanics:
    case StrategyType::hydrofracture:
    {
      MGRSpecialization specialization;
      specialization.coarseFlavor = ( strategy == StrategyType::thermalSinglePhasePoromechanics ||
                                      strategy == StrategyType::thermalSinglePhasePoromechanicsReservoirFVM ||
                                      strategy == StrategyType::thermalMultiphasePoromechanics )
                                      ? AMGFlavor::pressureTemperature
                                      : AMGFlavor::pressure;
      specialization.fRelaxAMGLevels = { LevelAMGBlock{ 0, AMGFlavor::displacementFiltered } };
      if( strategy == StrategyType::hydrofracture )
      {
        specialization.coarseMinCoarseSize = 1000;
      }
      return specialization;
    }
    case StrategyType::singlePhasePoromechanicsEmbeddedFractures:
    case StrategyType::singlePhasePoromechanicsConformingFractures:
    {
      MGRSpecialization specialization;
      specialization.coarseFlavor = AMGFlavor::pressure;
      specialization.fRelaxAMGLevels = { LevelAMGBlock{ 1, AMGFlavor::displacement } };
      specialization.cycle = "v(1,0)";
      return specialization;
    }
    case StrategyType::singlePhasePoromechanicsConformingFracturesALM:
    {
      MGRSpecialization specialization;
      specialization.coarseFlavor = AMGFlavor::pressure;
      specialization.cycle = "v(1,0)";
      return specialization;
    }
    case StrategyType::singlePhasePoromechanicsConformingFracturesALMReservoirFVM:
    {
      MGRSpecialization specialization;
      specialization.coarseFlavor = AMGFlavor::pressure;
      specialization.fRelaxAMGLevels = { LevelAMGBlock{ 1, AMGFlavor::almDisplacement } };
      specialization.cycle = "v(1,0)";
      return specialization;
    }
    case StrategyType::invalid:
    case StrategyType::singlePhaseReservoirFVM:
    case StrategyType::singlePhaseHybridFVM:
    case StrategyType::singlePhaseReservoirHybridFVM:
    case StrategyType::compositionalMultiphaseFVM:
    case StrategyType::compositionalMultiphaseHybridFVM:
    case StrategyType::compositionalMultiphaseReservoirFVM:
    case StrategyType::compositionalMultiphaseReservoirHybridFVM:
    case StrategyType::immiscibleMultiphaseFVM:
    case StrategyType::reactiveCompositionalMultiphaseOBL:
    {
      return {};
    }
  }

  GEOS_ERROR( GEOS_FMT( "Unsupported MGR strategy {}", strategy ) );
  return {};
}

std::optional< AMGFlavor > getFRelaxAMGFlavor( MGRSpecialization const & specialization,
                                               HYPRE_Int const level )
{
  for( LevelAMGBlock const & block : specialization.fRelaxAMGLevels )
  {
    if( block.level == level )
    {
      return block.flavor;
    }
  }
  return std::nullopt;
}

void appendAMGHeader( std::ostringstream & stream,
                      integer const indentLevel )
{
  appendLine( stream, indentLevel, "amg:" );
  appendLine( stream, indentLevel + 1, "tolerance: 0.0" );
  appendLine( stream, indentLevel + 1, "max_iter: 1" );
  appendLine( stream, indentLevel + 1, "print_level: 0" );
  // HYPRE_BoomerAMGCreate uses Schwarz (6) with zero extra smoother levels.
  // hypredrive's YAML AMG smoother defaults to ILU (5) and always writes that
  // onto the BoomerAMG object, so nested MGR AMG must pin the Create values.
  appendLine( stream, indentLevel + 1, "smoother:" );
  appendLine( stream, indentLevel + 2, "type: schwarz" );
  appendLine( stream, indentLevel + 2, "num_levels: 0" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 1" );
  appendLine( stream, indentLevel + 2, "ilu:" );
  appendIlUDisableRcm( stream, indentLevel + 3 );
}

void appendHypreBoomerAMGCreateMaxCoarseSize( std::ostringstream & stream,
                                              integer const indentLevel )
{
  // HYPRE_BoomerAMGCreate defaults max_coarse_size to 9. hypredrive's YAML AMG
  // defaults it to 64 and always writes that value unless the key is present,
  // so every nested MGR AMG block that would otherwise inherit hypre's Create
  // default must emit 9 explicitly. Leaving it out changes the first linear
  // solve, not just later Newton drift.
  appendLine( stream, indentLevel, "max_coarse_size: 9" );
}

void appendDisplacementAMG( std::ostringstream & stream,
                            integer const indentLevel,
                            integer const separateComponents,
                            bool const filterFunctions,
                            bool const useALMSmoother = false )
{
  appendAMGHeader( stream, indentLevel );
  appendLine( stream, indentLevel + 1, "coarsening:" );
  appendHypreBoomerAMGCreateMaxCoarseSize( stream, indentLevel + 2 );
  appendLine( stream, indentLevel + 2, "max_row_sum: 1.0" );
  appendLine( stream, indentLevel + 2, "strong_th: 0.6" );
  appendLine( stream, indentLevel + 2, "num_functions: 3" );
  appendLine( stream, indentLevel + 2, GEOS_FMT( "filter_functions: {}", filterFunctions ? separateComponents : 0 ) );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "type: 8" );
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "down_type: 16" );
  appendLine( stream, indentLevel + 2, "up_type: 16" );
  appendLine( stream, indentLevel + 2, "coarse_type: 16" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 1" );
#else
  appendLine( stream, indentLevel + 1, "relaxation:" );
  if( useALMSmoother )
  {
    appendLine( stream, indentLevel + 2, "down_type: l1sym-hgs" );
    appendLine( stream, indentLevel + 2, "up_type: l1sym-hgs" );
    appendLine( stream, indentLevel + 2, "coarse_type: ge" );
    appendLine( stream, indentLevel + 2, "num_sweeps: 1" );
    appendLine( stream, indentLevel + 2, "order: 0" );
  }
  else
  {
    appendLine( stream, indentLevel + 2, "order: 1" );
  }
#endif
}

void appendPressureAMG( std::ostringstream & stream,
                        integer const indentLevel,
                        HYPRE_Int const minCoarseSize )
{
  appendAMGHeader( stream, indentLevel );
  appendLine( stream, indentLevel + 1, "aggressive:" );
  appendLine( stream, indentLevel + 2, "num_levels: 1" );
  appendLine( stream, indentLevel + 2, "max_nnz_row: 20" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "prolongation_type: 7" );
#else
  appendLine( stream, indentLevel + 2, "prolongation_type: 4" );
#endif
  appendLine( stream, indentLevel + 1, "coarsening:" );
  if( minCoarseSize >= 0 )
  {
    appendLine( stream, indentLevel + 2, GEOS_FMT( "min_coarse_size: {}", minCoarseSize ) );
  }
  appendHypreBoomerAMGCreateMaxCoarseSize( stream, indentLevel + 2 );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "type: 8" );
  appendLine( stream, indentLevel + 2, "max_row_sum: 1.0" );
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "down_type: 18" );
  appendLine( stream, indentLevel + 2, "up_type: 18" );
  appendLine( stream, indentLevel + 2, "coarse_type: 18" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 2" );
#else
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "order: 1" );
#endif
}

void appendALMDisplacementFineAMG( std::ostringstream & stream,
                                   integer const indentLevel )
{
  appendAMGHeader( stream, indentLevel );
  appendLine( stream, indentLevel + 1, "aggressive:" );
  appendLine( stream, indentLevel + 2, "num_levels: 1" );
  appendLine( stream, indentLevel + 1, "interpolation:" );
  appendLine( stream, indentLevel + 2, "max_nnz_row: 20" );
  appendLine( stream, indentLevel + 1, "coarsening:" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "type: 8" );
#else
  appendLine( stream, indentLevel + 2, "type: falgout" );
#endif
  appendLine( stream, indentLevel + 2, "max_row_sum: 1.0" );
  appendLine( stream, indentLevel + 2, "strong_th: 0.8" );
  appendLine( stream, indentLevel + 2, "num_functions: 3" );
  appendLine( stream, indentLevel + 2, "filter_functions: 1" );
  appendLine( stream, indentLevel + 1, "relaxation:" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "down_type: 16" );
  appendLine( stream, indentLevel + 2, "up_type: 16" );
  appendLine( stream, indentLevel + 2, "coarse_type: 16" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 2" );
#else
  appendLine( stream, indentLevel + 2, "down_type: l1sym-hgs" );
  appendLine( stream, indentLevel + 2, "up_type: l1sym-hgs" );
  appendLine( stream, indentLevel + 2, "coarse_type: ge" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 2" );
  appendLine( stream, indentLevel + 2, "order: 0" );
#endif
}

void appendALMDisplacementBubbleAMG( std::ostringstream & stream,
                                     integer const indentLevel )
{
  appendAMGHeader( stream, indentLevel );
  appendLine( stream, indentLevel + 1, "interpolation:" );
  appendLine( stream, indentLevel + 2, "max_nnz_row: 10" );
  appendLine( stream, indentLevel + 1, "coarsening:" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "type: 8" );
#endif
  appendLine( stream, indentLevel + 2, "max_row_sum: 1.0" );
  appendLine( stream, indentLevel + 2, "strong_th: 0.75" );
  appendLine( stream, indentLevel + 2, "num_functions: 3" );
  appendLine( stream, indentLevel + 2, "filter_functions: 0" );
  appendLine( stream, indentLevel + 1, "relaxation:" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "down_type: 16" );
  appendLine( stream, indentLevel + 2, "up_type: 16" );
  appendLine( stream, indentLevel + 2, "coarse_type: 16" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 1" );
#else
  appendLine( stream, indentLevel + 2, "down_type: l1sym-hgs" );
  appendLine( stream, indentLevel + 2, "up_type: l1sym-hgs" );
  appendLine( stream, indentLevel + 2, "coarse_type: ge" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 1" );
  appendLine( stream, indentLevel + 2, "order: 0" );
#endif
}

void appendALMNestedMGR( std::ostringstream & stream,
                         integer const indentLevel,
                         stdVector< string > const & labelNames )
{
  stdVector< HYPRE_Int > const nodalLabels = { 0, 1, 2 };
  appendLine( stream, indentLevel, "mgr:" );
  appendLine( stream, indentLevel + 1, "tolerance: 0.0" );
  appendLine( stream, indentLevel + 1, "max_iter: 1" );
  appendLine( stream, indentLevel + 1, "print_level: 0" );
  appendLine( stream, indentLevel + 1, "cycle: v(1,0)" );
  appendLine( stream, indentLevel + 1, "num_levels: 2" );
  appendLine( stream, indentLevel + 1, "level:" );
  appendLine( stream, indentLevel + 2, "0:" );
  appendLine( stream, indentLevel + 3,
              GEOS_FMT( "f_dofs: [{}]", joinLabelNames( nodalLabels, labelNames ) ) );
  appendLine( stream, indentLevel + 3, "f_relaxation:" );
  appendALMDisplacementFineAMG( stream, indentLevel + 4 );
  appendLine( stream, indentLevel + 3, "g_relaxation: none" );
  appendLine( stream, indentLevel + 3, "prolongation_type: injection" );
  appendLine( stream, indentLevel + 3, "restriction_type: injection" );
  appendLine( stream, indentLevel + 3, "coarse_level_type: rap" );
  appendLine( stream, indentLevel + 1, "coarsest_level:" );
  appendALMDisplacementBubbleAMG( stream, indentLevel + 2 );
}

void appendPressureTemperatureAMG( std::ostringstream & stream,
                                   integer const indentLevel )
{
  appendAMGHeader( stream, indentLevel );
  appendLine( stream, indentLevel + 1, "aggressive:" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "num_levels: 0" );
#else
  appendLine( stream, indentLevel + 2, "num_levels: 1" );
#endif
  appendLine( stream, indentLevel + 2, "max_nnz_row: 16" );
  appendLine( stream, indentLevel + 1, "coarsening:" );
  appendHypreBoomerAMGCreateMaxCoarseSize( stream, indentLevel + 2 );
  appendLine( stream, indentLevel + 2, "num_functions: 2" );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 2, "type: 8" );
  appendLine( stream, indentLevel + 2, "max_row_sum: 1.0" );
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "down_type: 18" );
  appendLine( stream, indentLevel + 2, "up_type: 18" );
  appendLine( stream, indentLevel + 2, "coarse_type: 18" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 2" );
#else
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "order: 1" );
#endif
}

void appendNamedMGRRelaxation( std::ostringstream & stream,
                               integer const indentLevel,
                               char const * const key,
                               char const * const typeName,
                               HYPRE_Int const numSweeps )
{
  std::string const type( typeName );

  // hypredrive attaches HYPRE_ILUCreate for YAML `ilu`. Match hypre MGR's
  // built-in global smoother type 16 (par_mgr_setup.c): BJ-ILUK, fill 0,
  // max_iter = sweep count, tol 0, no RCM, and hypre_ILUCreate's max_row_nnz
  // of 1000. A scalar `g_relaxation: ilu` would keep hypredrive's ILU YAML
  // defaults instead (including max_row_nnz 200).
  if( type == "ilu" )
  {
    HYPRE_Int const iluMaxIter = ( numSweeps > 0 ) ? numSweeps : 1;
    appendLine( stream, indentLevel, GEOS_FMT( "{}:", key ) );
    appendLine( stream, indentLevel + 1, GEOS_FMT( "num_sweeps: {}", numSweeps ) );
    appendLine( stream, indentLevel + 1, "ilu:" );
    appendLine( stream, indentLevel + 2, "type: 0" );
    appendLine( stream, indentLevel + 2, "fill_level: 0" );
    appendLine( stream, indentLevel + 2, GEOS_FMT( "max_iter: {}", iluMaxIter ) );
    appendLine( stream, indentLevel + 2, "tolerance: 0.0" );
    appendIlUDisableRcm( stream, indentLevel + 2 );
    appendLine( stream, indentLevel + 2, "print_level: 0" );
    appendLine( stream, indentLevel + 2, "max_row_nnz: 1000" );
    return;
  }

  // hypredrive defaults num_sweeps to 1 even when the type is none. Legacy
  // setReduction forces 0 sweeps once F/G-relaxation is disabled. Quote the
  // type so YAML 1.1 does not mangle an unquoted `none` token.
  if( type == "none" )
  {
    appendLine( stream, indentLevel, GEOS_FMT( "{}:", key ) );
    appendLine( stream, indentLevel + 1, "type: \"none\"" );
    appendLine( stream, indentLevel + 1, GEOS_FMT( "num_sweeps: {}", numSweeps ) );
    return;
  }

  // hypredrive defaults num_sweeps to 1. When GEOS uses a different sweep
  // count, emit a mapping so the legacy HYPRE_MGRSetLevel*Iters value is
  // preserved.
  if( numSweeps == 1 )
  {
    appendLine( stream, indentLevel, GEOS_FMT( "{}: {}", key, typeName ) );
  }
  else
  {
    appendLine( stream, indentLevel, GEOS_FMT( "{}:", key ) );
    appendLine( stream, indentLevel + 1, GEOS_FMT( "type: {}", typeName ) );
    appendLine( stream, indentLevel + 1, GEOS_FMT( "num_sweeps: {}", numSweeps ) );
  }
}

void appendAMGByFlavor( std::ostringstream & stream,
                        integer const indentLevel,
                        AMGFlavor const flavor,
                        LinearSolverParameters::MGR const & mgrParams,
                        HYPRE_Int const coarseMinCoarseSize = -1 )
{
  switch( flavor )
  {
    case AMGFlavor::pressure:
    {
      appendPressureAMG( stream, indentLevel, coarseMinCoarseSize );
      break;
    }
    case AMGFlavor::pressureTemperature:
    {
      appendPressureTemperatureAMG( stream, indentLevel );
      break;
    }
    case AMGFlavor::displacementFiltered:
    {
      appendDisplacementAMG( stream, indentLevel, mgrParams.separateComponents, true );
      break;
    }
    case AMGFlavor::displacement:
    {
      appendDisplacementAMG( stream, indentLevel, 0, false );
      break;
    }
    case AMGFlavor::almDisplacement:
    {
      appendDisplacementAMG( stream, indentLevel, mgrParams.separateComponents, true, true );
      break;
    }
  }
}

template< typename STRATEGY >
class MGRStrategyProbe final : public STRATEGY
{
public:
  explicit MGRStrategyProbe( arrayView1d< int const > const & numComponentsPerField )
    : STRATEGY( numComponentsPerField )
  {}

  using STRATEGY::m_coarseGridThreshold;
  using STRATEGY::m_labels;
  using STRATEGY::m_levelCoarseGridMethod;
  using STRATEGY::m_levelFRelaxIters;
  using STRATEGY::m_levelFRelaxType;
  using STRATEGY::m_levelGlobalSmootherIters;
  using STRATEGY::m_levelGlobalSmootherType;
  using STRATEGY::m_levelInterpType;
  using STRATEGY::m_levelRestrictType;
  using STRATEGY::m_numBlocks;

  static constexpr HYPRE_Int numReductionLevels = STRATEGY::numLevels;
};

template< typename STRATEGY >
bool buildStrategyYaml( LinearSolverParameters const & params,
                        stdVector< string > const & labelNames,
                        arrayView1d< int const > const & numComponentsPerField,
                        std::string & preconditionerYaml )
{
  MGRStrategyProbe< STRATEGY > strategy( numComponentsPerField );
  MGRSpecialization const specialization = getSpecialization( params.mgr.strategy );

  HyprePrecWrapper precond;
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRCreate( &precond.ptr ) );

  HypreMGRData mgrData;
  mgrData.pointMarkers.resize( strategy.m_numBlocks );
  std::iota( mgrData.pointMarkers.begin(), mgrData.pointMarkers.end(), 0 );

  strategy.setup( params.mgr, precond, mgrData );

  std::ostringstream stream;
  if( params.mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFracturesALM )
  {
    // Keep this representation in lockstep with the fully coupled ALM
    // strategy: the outer F block is displacement plus bubble displacement,
    // and its F-relaxation is a two-level nested MGR.
    appendLine( stream, 0, "preconditioner:" );
    appendLine( stream, 1, "mgr:" );
    appendLine( stream, 2, "tolerance: 0.0" );
    appendLine( stream, 2, "max_iter: 1" );
    appendLine( stream, 2, GEOS_FMT( "print_level: {}", getMGRPrintLevel( params.logLevel ) ) );
    appendLine( stream, 2, "cycle: v(1,0)" );
    appendLine( stream, 2, "non_c_to_f: 1" );
    appendLine( stream, 2, "nonglk_max_elmts: 1" );
    appendLine( stream, 2, "pmax: 0" );
    appendLine( stream, 2, GEOS_FMT( "coarse_th: {}", strategy.m_coarseGridThreshold ) );
    appendLine( stream, 2, "num_levels: 2" );
    appendLine( stream, 2, "level:" );
    appendLine( stream, 3, "0:" );
    stdVector< HYPRE_Int > const displacementLabels = { 0, 1, 2, 3, 4, 5 };
    appendLine( stream, 4,
                GEOS_FMT( "f_dofs: [{}]", joinLabelNames( displacementLabels, labelNames ) ) );
    appendLine( stream, 4, "f_relaxation:" );
    appendALMNestedMGR( stream, 5, labelNames );
    appendLine( stream, 4, "g_relaxation: none" );
    appendLine( stream, 4, "prolongation_type: 2" );
    appendLine( stream, 4, "restriction_type: injection" );
    appendLine( stream, 4, "coarse_level_type: rap" );
    appendLine( stream, 2, "coarsest_level:" );
    appendPressureAMG( stream, 3, -1 );

    preconditionerYaml = stream.str();

    destroyWrapper( mgrData.coarseSolver );
    destroyWrapper( mgrData.mechSolver );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRDestroy( precond.ptr ) );
    destroyWrapper( mgrData.nestedSolver );
    return true;
  }

  appendLine( stream, 0, "preconditioner:" );
  appendLine( stream, 1, "mgr:" );
  appendLine( stream, 2, "tolerance: 0.0" );
  appendLine( stream, 2, "max_iter: 1" );
  // Do not emit relax_type. HYPRE_MGRCreate defaults it to 0, but hypredrive's
  // YAML validator rejects 0 (allowed names map to 3–18, default Jacobi 7).
  // hypre's MGR setup also rewrites F-relax type 0 to 7 unless interp is 12,
  // so omitting the key matches the effective legacy behavior.
  appendLine( stream, 2, GEOS_FMT( "print_level: {}", getMGRPrintLevel( params.logLevel ) ) );
  if( specialization.cycle != nullptr )
  {
    appendLine( stream, 2, GEOS_FMT( "cycle: {}", specialization.cycle ) );
  }
  appendLine( stream, 2, GEOS_FMT( "non_c_to_f: {}", 1 ) );
  appendLine( stream, 2, GEOS_FMT( "nonglk_max_elmts: {}", 1 ) );
  appendLine( stream, 2, GEOS_FMT( "pmax: {}", specialization.pmax ) );
  appendLine( stream, 2, GEOS_FMT( "coarse_th: {}", strategy.m_coarseGridThreshold ) );
  appendLine( stream, 2, GEOS_FMT( "num_levels: {}", MGRStrategyProbe< STRATEGY >::numReductionLevels + 1 ) );
  appendLine( stream, 2, "level:" );

  stdVector< HYPRE_Int > activeLabels( static_cast< size_t >( strategy.m_numBlocks ) );
  std::iota( activeLabels.begin(), activeLabels.end(), 0 );

  for( HYPRE_Int level = 0; level < MGRStrategyProbe< STRATEGY >::numReductionLevels; ++level )
  {
    stdVector< HYPRE_Int > cLabels( strategy.m_labels[level].begin(), strategy.m_labels[level].end() );
    stdVector< HYPRE_Int > fLabels;

    for( HYPRE_Int const activeLabel : activeLabels )
    {
      if( std::find( cLabels.begin(), cLabels.end(), activeLabel ) == cLabels.end() )
      {
        fLabels.push_back( activeLabel );
      }
    }

    appendLine( stream, 3, GEOS_FMT( "{}:", level ) );
    appendLine( stream, 4, GEOS_FMT( "f_dofs: [{}]", joinLabelNames( fLabels, labelNames ) ) );
    if( strategy.m_levelFRelaxType[level] == hypre::MGRFRelaxationType::amgVCycle )
    {
      std::optional< AMGFlavor > const amgFlavor = getFRelaxAMGFlavor( specialization, level );
      GEOS_ERROR_IF( !amgFlavor.has_value(),
                     GEOS_FMT( "Missing hypredrive AMG block description for MGR strategy {} level {}",
                               params.mgr.strategy, level ) );
      appendLine( stream, 4, "f_relaxation:" );
      appendLine( stream, 5, GEOS_FMT( "num_sweeps: {}", strategy.m_levelFRelaxIters[level] ) );
      appendAMGByFlavor( stream, 5, *amgFlavor, params.mgr );
    }
    else
    {
      char const * const fRelaxationName = getGeneratedMGRFRelaxationName( strategy.m_levelFRelaxType[level] );
      if( fRelaxationName == nullptr )
      {
        return false;
      }
      appendNamedMGRRelaxation( stream, 4, "f_relaxation", fRelaxationName, strategy.m_levelFRelaxIters[level] );
    }

    char const * const gRelaxationName = getGeneratedMGRGlobalSmootherName( strategy.m_levelGlobalSmootherType[level] );
    if( gRelaxationName == nullptr )
    {
      return false;
    }
    appendNamedMGRRelaxation( stream, 4, "g_relaxation", gRelaxationName, strategy.m_levelGlobalSmootherIters[level] );
    appendLine( stream, 4, GEOS_FMT( "prolongation_type: {}", static_cast< HYPRE_Int >( strategy.m_levelInterpType[level] ) ) );
    appendLine( stream, 4, GEOS_FMT( "restriction_type: {}", static_cast< HYPRE_Int >( strategy.m_levelRestrictType[level] ) ) );
    appendLine( stream, 4, GEOS_FMT( "coarse_level_type: {}", static_cast< HYPRE_Int >( strategy.m_levelCoarseGridMethod[level] ) ) );

    activeLabels = std::move( cLabels );
  }

  appendLine( stream, 2, "coarsest_level:" );
  appendAMGByFlavor( stream, 3, specialization.coarseFlavor, params.mgr, specialization.coarseMinCoarseSize );

  preconditionerYaml = stream.str();

  destroyWrapper( mgrData.coarseSolver );
  destroyWrapper( mgrData.mechSolver );
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRDestroy( precond.ptr ) );
  destroyWrapper( mgrData.nestedSolver );

  return true;
}

bool buildMGRPreconditionerYaml( LinearSolverParameters const & params,
                                 stdVector< string > const & fieldNames,
                                 arrayView1d< int const > const & numComponentsPerField,
                                 std::string & linearSystemYaml,
                                 std::string & preconditionerYaml )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  if( numComponentsPerField.empty() )
  {
    return false;
  }

  stdVector< string > const labelNames = buildDofLabelNames( params.mgr.strategy, fieldNames, numComponentsPerField );
  linearSystemYaml = buildDofLabelsYaml( labelNames );

  switch( params.mgr.strategy )
  {
    case StrategyType::singlePhaseReservoirFVM:
      return buildStrategyYaml< hypre::mgr::SinglePhaseReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::thermalSinglePhaseReservoirFVM:
      return buildStrategyYaml< hypre::mgr::ThermalSinglePhaseReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhaseHybridFVM:
      return buildStrategyYaml< hypre::mgr::SinglePhaseHybridFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhaseReservoirHybridFVM:
      return buildStrategyYaml< hypre::mgr::SinglePhaseReservoirHybridFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhasePoromechanics:
      return buildStrategyYaml< hypre::mgr::SinglePhasePoromechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::thermalSinglePhasePoromechanics:
      return buildStrategyYaml< hypre::mgr::ThermalSinglePhasePoromechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::hybridSinglePhasePoromechanics:
      return buildStrategyYaml< hypre::mgr::HybridSinglePhasePoromechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhasePoromechanicsEmbeddedFractures:
      return buildStrategyYaml< hypre::mgr::SinglePhasePoromechanicsEmbeddedFractures >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhasePoromechanicsConformingFractures:
      return buildStrategyYaml< hypre::mgr::SinglePhasePoromechanicsConformingFractures >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhasePoromechanicsConformingFracturesALM:
      return buildStrategyYaml< hypre::mgr::SinglePhasePoromechanicsConformingFracturesALM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhasePoromechanicsConformingFracturesALMReservoirFVM:
      return buildStrategyYaml< hypre::mgr::SinglePhasePoromechanicsConformingFracturesALMReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::singlePhasePoromechanicsReservoirFVM:
      return buildStrategyYaml< hypre::mgr::SinglePhasePoromechanicsReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::thermalSinglePhasePoromechanicsReservoirFVM:
      return buildStrategyYaml< hypre::mgr::ThermalSinglePhasePoromechanicsReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::compositionalMultiphaseFVM:
      return buildStrategyYaml< hypre::mgr::CompositionalMultiphaseFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::compositionalMultiphaseHybridFVM:
      return buildStrategyYaml< hypre::mgr::CompositionalMultiphaseHybridFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::compositionalMultiphaseReservoirFVM:
      return buildStrategyYaml< hypre::mgr::CompositionalMultiphaseReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::compositionalMultiphaseReservoirHybridFVM:
      return buildStrategyYaml< hypre::mgr::CompositionalMultiphaseReservoirHybridFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::immiscibleMultiphaseFVM:
      return buildStrategyYaml< hypre::mgr::ImmiscibleMultiphaseFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::reactiveCompositionalMultiphaseOBL:
      return buildStrategyYaml< hypre::mgr::ReactiveCompositionalMultiphaseOBL >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::thermalCompositionalMultiphaseFVM:
      return buildStrategyYaml< hypre::mgr::ThermalCompositionalMultiphaseFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::thermalCompositionalMultiphaseReservoirFVM:
      return buildStrategyYaml< hypre::mgr::ThermalCompositionalMultiphaseReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::multiphasePoromechanics:
      return buildStrategyYaml< hypre::mgr::MultiphasePoromechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::multiphasePoromechanicsReservoirFVM:
      return buildStrategyYaml< hypre::mgr::MultiphasePoromechanicsReservoirFVM >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::thermalMultiphasePoromechanics:
      return buildStrategyYaml< hypre::mgr::ThermalMultiphasePoromechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::hydrofracture:
      return buildStrategyYaml< hypre::mgr::Hydrofracture >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::lagrangianContactMechanics:
      return buildStrategyYaml< hypre::mgr::LagrangianContactMechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::augmentedLagrangianContactMechanics:
      return buildStrategyYaml< hypre::mgr::AugmentedLagrangianContactMechanics >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::lagrangianContactMechanicsBubbleStab:
      return buildStrategyYaml< hypre::mgr::LagrangianContactMechanicsBubbleStabilization >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::solidMechanicsEmbeddedFractures:
      return buildStrategyYaml< hypre::mgr::SolidMechanicsEmbeddedFractures >( params, labelNames, numComponentsPerField, preconditionerYaml );
    case StrategyType::invalid:
      return false;
  }

  return false;
}

bool buildGeneratedInputArgsParseTarget( LinearSolverParameters const & params,
                                         stdVector< string > const & fieldNames,
                                         arrayView1d< int const > const & numComponentsPerField,
                                         hypre::hypredrive::InputArgsParseTarget & target )
{
  std::string linearSystemYaml;
  std::string preconditionerYaml;

  switch( params.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::none:
      preconditionerYaml = "preconditioner: none\n";
      break;
    case LinearSolverParameters::PreconditionerType::amg:
      if( !buildAMGPreconditionerYaml( params, preconditionerYaml ) )
      {
        return false;
      }
      break;
    case LinearSolverParameters::PreconditionerType::mgr:
      if( !buildMGRPreconditionerYaml( params, fieldNames, numComponentsPerField, linearSystemYaml, preconditionerYaml ) )
      {
        return false;
      }
      break;
    case LinearSolverParameters::PreconditionerType::iluk:
    case LinearSolverParameters::PreconditionerType::ilut:
      if( !buildILUPreconditionerYaml( params, preconditionerYaml ) )
      {
        return false;
      }
      break;
    default:
      return false;
  }

  std::ostringstream stream;
  if( !linearSystemYaml.empty() )
  {
    stream << linearSystemYaml;
  }
  appendSolverYaml( stream, params );
  stream << preconditionerYaml;

  target.source = hypre::hypredrive::InputSource::generatedFallback;
  target.argument = stream.str();
  return true;
}

}

namespace hypre
{

namespace hypredrive
{

bool shouldUse( LinearSolverParameters const & params )
{
  if( envFlagEnabled( "GEOS_HYPREDRV_FORCE_LEGACY" ) )
  {
    return false;
  }

  if( !params.hypredriveInputFile.empty() )
  {
    return true;
  }

  if( params.krylov.useAdaptiveTol )
  {
    return false;
  }

  return getSolverName( params.solverType ) != nullptr &&
         supportsGeneratedPreconditioner( params.preconditionerType );
}

bool buildInputArgsParseTarget( LinearSolverParameters const & params,
                                InputArgsParseTarget & target )
{
  stdVector< string > fieldNames;
  array1d< int > numComponentsPerField;
  return buildInputArgsParseTarget( params, fieldNames, numComponentsPerField, target );
}

bool buildInputArgsParseTarget( LinearSolverParameters const & params,
                                stdVector< string > const & fieldNames,
                                arrayView1d< int const > const & numComponentsPerField,
                                InputArgsParseTarget & target )
{
  if( !params.hypredriveInputFile.empty() )
  {
    target.source = InputSource::authoritativeFile;
    target.argument = params.hypredriveInputFile;
    return true;
  }

  if( !shouldUse( params ) )
  {
    return false;
  }

  return buildGeneratedInputArgsParseTarget( params, fieldNames, numComponentsPerField, target );
}

std::string formatInputArgsParseTargetYaml( InputArgsParseTarget const & target )
{
  if( target.source == InputSource::authoritativeFile )
  {
    std::optional< std::string > const contents = readTextFileContents( target.argument );
    if( contents.has_value() )
    {
      return GEOS_FMT( "# hypredrive authoritative file: {}\n{}",
                       target.argument,
                       *contents );
    }

    return GEOS_FMT( "# hypredrive authoritative file: {}\n"
                     "# GEOS was unable to read this YAML file while formatting the hypredrive input dump.",
                     target.argument );
  }

  return GEOS_FMT( "# hypredrive input generated by GEOS\n{}",
                   target.argument );
}

bool wasInputArgsParseTargetLogged( InputArgsParseTarget const & target )
{
  return loggedInputArgsParseTargets().count( makeInputArgsParseTargetSignature( target ) ) > 0;
}

void markInputArgsParseTargetLogged( InputArgsParseTarget const & target )
{
  loggedInputArgsParseTargets().insert( makeInputArgsParseTargetSignature( target ) );
}

void logInputArgsParseTarget( LinearSolverParameters const & params,
                              InputArgsParseTarget const & target )
{
  if( params.logLevel < 1 )
  {
    return;
  }

  if( wasInputArgsParseTargetLogged( target ) )
  {
    return;
  }

  markInputArgsParseTargetLogged( target );

  if( target.source == InputSource::authoritativeFile )
  {
    std::optional< std::string > const contents = readTextFileContents( target.argument );
    if( contents.has_value() )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "        hypredrive input | authoritative file | {}\n{}",
                                 target.argument,
                                 *contents ) );
    }
    else
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "        hypredrive input | authoritative file | {}\n"
                                 "        hypredrive input | unable to read authoritative YAML from GEOS",
                                 target.argument ) );
    }
  }
  else
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "        hypredrive input | generated fallback\n{}",
                               target.argument ) );
  }
}

void initializeRuntime()
{
  checkHypredriveCall( HYPREDRV_Initialize(), "HYPREDRV_Initialize" );
}

void finalizeRuntime()
{
  checkHypredriveCall( HYPREDRV_Finalize(), "HYPREDRV_Finalize" );
}

}

}

namespace
{

std::string makeConfigurationSignature( hypre::hypredrive::InputArgsParseTarget const & target )
{
  return GEOS_FMT( "{}\n{}",
                   target.source == hypre::hypredrive::InputSource::authoritativeFile ? "authoritativeFile" : "generatedFallback",
                   hypre::hypredrive::formatInputArgsParseTargetYaml( target ) );
}

std::string makeStructureSignature( HypreMatrix const & mat,
                                    stdVector< string > const & fieldNames,
                                    array1d< int > const & numComponentsPerField,
                                    array1d< int > const & pointMarkers,
                                    LinearSolverExecutionContext const * const executionContext )
{
  std::ostringstream signature;
  signature << "systemSetupTimestamp: "
            << ( executionContext != nullptr ? executionContext->systemSetupTimestamp : 0 ) << '\n';
  signature << "localRows: " << mat.numLocalRows() << '\n';
  signature << "localCols: " << mat.numLocalCols() << '\n';
  signature << "globalRows: " << mat.numGlobalRows() << '\n';
  signature << "globalCols: " << mat.numGlobalCols() << '\n';
  signature << "fieldNames:";
  for( string const & fieldName : fieldNames )
  {
    signature << ' ' << fieldName;
  }
  signature << '\n';
  signature << "numComponentsPerField:";
  for( int const numComponents : numComponentsPerField )
  {
    signature << ' ' << numComponents;
  }
  signature << '\n';
  signature << "pointMarkers:";
  for( int const pointMarker : pointMarkers )
  {
    signature << ' ' << pointMarker;
  }
  signature << '\n';

  return signature.str();
}

std::string makeTimestepScopeName( LinearSolverExecutionContext const & context )
{
  return GEOS_FMT( "timestep-{}-{}-{}",
                   context.cycleNumber,
                   context.timeStepAttempt,
                   context.configurationAttempt );
}

std::string makeNewtonScopeName( LinearSolverExecutionContext const & context )
{
  return GEOS_FMT( "newton-{}-{}-{}-{}",
                   context.cycleNumber,
                   context.timeStepAttempt,
                   context.configurationAttempt,
                   context.nonlinearIteration );
}

}

HypredriveSolver::HypredriveSolver( LinearSolverParameters parameters )
  : Base( std::move( parameters ) )
{}

HypredriveSolver::~HypredriveSolver()
{
  HypredriveSolver::clear();
}

void HypredriveSolver::setExecutionContext( LinearSolverExecutionContext const & context )
{
  m_executionContext = context;
  m_hasExecutionContext = true;
}

void HypredriveSolver::setNearNullKernel( arrayView1d< HypreVector const > const & nearNullKernel )
{
  m_nearNullKernel = nearNullKernel;
}

void HypredriveSolver::setup( HypreMatrix const & mat )
{
  Base::setup( mat );

  bool const configured = configureHypredrive( mat );

  // Falling back to the legacy solver must be an all-or-nothing decision: the two paths
  // execute different sequences of collective operations, so if some ranks kept hypredrive
  // while others fell back, the next solve would deadlock. Configuration can legitimately
  // fail on a subset of ranks (e.g. a rank that owns no rows of the target region has no
  // dof components to describe), so agree on the outcome across the communicator.
  int const allRanksConfigured = MpiWrapper::min( configured ? 1 : 0, mat.comm() );

  if( !allRanksConfigured )
  {
    resetHypredriveState();
    setupLegacy( mat );
  }
}

void HypredriveSolver::createHypredrive( HypreMatrix const & mat,
                                         hypre::hypredrive::InputArgsParseTarget const & parseTarget,
                                         std::string const & configurationSignature,
                                         std::string const & structureSignature,
                                         arrayView1d< int > const & pointMarkers )
{
  resetHypredriveState();

  checkHypredriveCall( HYPREDRV_Create( mat.comm(), &m_hypredrive ), "HYPREDRV_Create" );
  checkHypredriveCall( HYPREDRV_SetLibraryMode( m_hypredrive ), "HYPREDRV_SetLibraryMode" );

  char * argv[] = { const_cast< char * >( parseTarget.argument.c_str() ) };
  std::string const parseContext =
    parseTarget.source == hypre::hypredrive::InputSource::authoritativeFile
    ? GEOS_FMT( "authoritative YAML file '{}'", parseTarget.argument )
    : "YAML generated by GEOS";
  checkHypredriveCall( HYPREDRV_InputArgsParse( 1, argv, m_hypredrive ),
                       "HYPREDRV_InputArgsParse",
                       parseContext );
  if( parseTarget.source == hypre::hypredrive::InputSource::generatedFallback &&
      m_hasExecutionContext &&
      !m_executionContext.solverName.empty() )
  {
    checkHypredriveCall( HYPREDRV_ObjectSetName( m_hypredrive,
                                                 m_executionContext.solverName.c_str() ),
                         "HYPREDRV_ObjectSetName" );
  }

  ++m_hypredriveGeneration;
  m_configurationSignature = configurationSignature;
  m_structureSignature = structureSignature;

  refreshBoundObjects( mat, pointMarkers );
}

void HypredriveSolver::refreshBoundObjects( HypreMatrix const & mat,
                                            arrayView1d< int > const & pointMarkers )
{
  checkHypredriveCall( HYPREDRV_LinearSystemSetMatrix( m_hypredrive,
                                                       toHypreMatrix( mat.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetMatrix" );

  // HypreVector::create is collective, so the decision to recreate the bound vectors
  // must be made globally: a repartitioning can leave some ranks with an unchanged
  // local size while others change, and letting each rank decide on its own would
  // have only a subset of them enter the collective.
  int const rankNeedsRecreate =
    ( !m_dummyRhs.ready() || m_dummyRhs.localSize() != mat.numLocalRows() ) ? 1 : 0;
  if( MpiWrapper::max( rankNeedsRecreate, mat.comm() ) )
  {
    m_dummyRhs.reset();
    m_dummySol.reset();
    m_residual.reset();
    m_dummyRhs.create( mat.numLocalRows(), mat.comm() );
    m_dummySol.create( mat.numLocalRows(), mat.comm() );
    m_residual.create( mat.numLocalRows(), mat.comm() );
  }

  checkHypredriveCall( HYPREDRV_LinearSystemSetRHS( m_hypredrive,
                                                    reinterpret_cast< HYPRE_Vector >( m_dummyRhs.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetRHS" );
  checkHypredriveCall( HYPREDRV_LinearSystemSetSolution( m_hypredrive,
                                                         reinterpret_cast< HYPRE_Vector >( m_dummySol.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetSolution" );

  // This must be called on every rank, including ranks that own no local rows:
  // hypredrive builds the global dof-label set collectively inside this call, so
  // skipping it where the local dofmap is empty leaves those ranks out of the
  // collective and the subsequent solver creation fails there.
  //
  // Library-mode hypredrive uses these labels to configure its MGR hierarchy.
  // The legacy HypreSolver path deliberately keeps its setup dummy untagged because
  // its solve receives caller-owned rhs and solution vectors.
  checkHypredriveCall( HYPREDRV_LinearSystemSetDofmap( m_hypredrive,
                                                       LvArray::integerConversion< int >( pointMarkers.size() ),
                                                       pointMarkers.data() ),
                       "HYPREDRV_LinearSystemSetDofmap" );

  if( !m_nearNullKernel.empty() )
  {
    localIndex const numEntries = mat.numLocalRows();
    localIndex const numModes = m_nearNullKernel.size();
    array1d< HYPRE_Complex > values;
    values.resizeWithoutInitializationOrDestruction( hypre::memorySpace, numEntries * numModes );

    for( localIndex mode = 0; mode < numModes; ++mode )
    {
      GEOS_ERROR_IF( m_nearNullKernel[mode].localSize() != numEntries,
                     "HypreDrive near-null-space vector size does not match the matrix local size" );
      GEOS_LAI_CHECK_ERROR(
        HYPRE_IJVectorGetValues( m_nearNullKernel[mode].unwrappedIJ(),
                                 LvArray::integerConversion< HYPRE_Int >( numEntries ),
                                 nullptr,
                                 values.data() + mode * numEntries ) );
    }

    values.registerTouch( hypre::memorySpace );
    values.move( hostMemorySpace );

    checkHypredriveCall(
      HYPREDRV_LinearSystemSetNearNullSpace(
        m_hypredrive,
        LvArray::integerConversion< int >( numEntries ),
        LvArray::integerConversion< int >( numModes ),
        values.data() ),
      "HYPREDRV_LinearSystemSetNearNullSpace" );
  }
}

void HypredriveSolver::setupLegacy( HypreMatrix const & mat )
{
  if( !m_legacySolver )
  {
    m_legacySolver = std::make_unique< HypreSolver >( m_params );
  }

  m_legacySolver->setup( mat );
  syncLegacyResult();
}

bool HypredriveSolver::configureHypredrive( HypreMatrix const & mat )
{
  if( !hypre::hypredrive::shouldUse( m_params ) )
  {
    return false;
  }

  DofManager const * const dofManager = mat.dofManager();
  stdVector< string > fieldNames;
  array1d< int > numComponentsPerField;
  array1d< int > pointMarkers;
  if( dofManager != nullptr )
  {
    fieldNames = dofManager->fieldNames();
    numComponentsPerField = dofManager->numComponentsPerField();
  }
  hypre::fillKrylovDofLabels( mat, pointMarkers );

  hypre::hypredrive::InputArgsParseTarget parseTarget;
  if( !hypre::hypredrive::buildInputArgsParseTarget( m_params,
                                                     fieldNames,
                                                     numComponentsPerField,
                                                     parseTarget ) )
  {
    if( !m_reportedGeneratedYamlFailure )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "Warning: {}: hypredrive input generation failed for the current "
                                 "linear-solver configuration; falling back to the legacy hypre solver",
                                 ( m_hasExecutionContext && !m_executionContext.solverName.empty() )
                                   ? m_executionContext.solverName : "linear solver" ) );
      m_reportedGeneratedYamlFailure = true;
    }
    return false;
  }

  std::string const configurationSignature = makeConfigurationSignature( parseTarget );
  std::string const structureSignature =
    makeStructureSignature( mat,
                            fieldNames,
                            numComponentsPerField,
                            pointMarkers,
                            m_hasExecutionContext ? &m_executionContext : nullptr );

  hypre::hypredrive::logInputArgsParseTarget( m_params, parseTarget );

  if( m_legacySolver )
  {
    m_legacySolver->clear();
    m_legacySolver.reset();
  }

  // Legacy HypreSolver destroys and recreates MGR on every setup. Reusing a
  // HYPREDRV handle across Newton steps can keep MGR/AMG/ILU bookkeeping even
  // when YAML reuse is off, which shows up as later-solve iteration drift.
  // AMG-only generated configs keep the handle so unit tests can assert reuse.
  bool const recreateHandle = ( m_hypredrive == nullptr ) ||
                              ( m_configurationSignature != configurationSignature ) ||
                              ( m_structureSignature != structureSignature ) ||
                              ( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr );

  if( recreateHandle )
  {
    createHypredrive( mat,
                      parseTarget,
                      configurationSignature,
                      structureSignature,
                      pointMarkers.toView() );
  }
  else
  {
    refreshBoundObjects( mat, pointMarkers.toView() );
  }

  syncExecutionAnnotations();

  // As on the legacy hypre path, floating point exceptions must be disabled while
  // hypre builds the preconditioner (benign FPEs occur, e.g. in MGR interpolation setup).
  LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );

  if( m_linearSolverCreated )
  {
    // HYPREDRV_LinearSolverCreate does not replace an existing solver and
    // preconditioner. Destroy them before rebuilding a reused handle so
    // repeated setup cycles do not retain Hypre/MGR/ILU allocations.
    checkHypredriveCall( HYPREDRV_LinearSolverDestroy( m_hypredrive ),
                         "HYPREDRV_LinearSolverDestroy" );
    m_linearSolverCreated = false;
  }

  checkHypredriveCall( HYPREDRV_AnnotateBegin( m_hypredrive, "system", -1 ),
                       "HYPREDRV_AnnotateBegin" );
  checkHypredriveCall( HYPREDRV_LinearSolverCreate( m_hypredrive ), "HYPREDRV_LinearSolverCreate" );
  m_linearSolverCreated = true;
  checkHypredriveCall( HYPREDRV_LinearSolverSetup( m_hypredrive ), "HYPREDRV_LinearSolverSetup" );
  checkHypredriveCall( HYPREDRV_AnnotateEnd( m_hypredrive, "system", -1 ),
                       "HYPREDRV_AnnotateEnd" );
  checkHypredriveCall( HYPREDRV_LinearSolverGetSetupTime( m_hypredrive, &m_result.setupTime ),
                       "HYPREDRV_LinearSolverGetSetupTime" );

  return true;
}

void HypredriveSolver::apply( HypreVector const & src,
                              HypreVector & dst ) const
{
  if( m_legacySolver )
  {
    m_legacySolver->apply( src, dst );
    return;
  }

  applyHypredrive( src, dst );
}

void HypredriveSolver::applyHypredrive( HypreVector const & rhs,
                                        HypreVector & sol ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( rhs.ready() );
  GEOS_LAI_ASSERT( sol.ready() );

  m_dummySol.copy( sol );
  checkHypredriveCall( HYPREDRV_LinearSystemSetInitialGuess( m_hypredrive,
                                                             reinterpret_cast< HYPRE_Vector >( m_dummySol.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetInitialGuess" );
  checkHypredriveCall( HYPREDRV_LinearSystemSetRHS( m_hypredrive,
                                                    reinterpret_cast< HYPRE_Vector >( rhs.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetRHS" );
  checkHypredriveCall( HYPREDRV_LinearSystemSetSolution( m_hypredrive,
                                                         reinterpret_cast< HYPRE_Vector >( sol.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetSolution" );
  checkHypredriveCall( HYPREDRV_LinearSystemResetInitialGuess( m_hypredrive ),
                       "HYPREDRV_LinearSystemResetInitialGuess" );

  {
    // As during setup, hypre's Krylov and MGR kernels can raise benign floating point
    // exceptions while solving (e.g. norms of empty local blocks on ranks that own no
    // rows), so they must be masked here as well.
    LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );
    checkHypredriveCall( HYPREDRV_LinearSolverApply( m_hypredrive ), "HYPREDRV_LinearSolverApply" );
  }

  sol.touch();
}

void HypredriveSolver::solve( HypreVector const & rhs,
                              HypreVector & sol ) const
{
  if( m_legacySolver )
  {
    m_legacySolver->solve( rhs, sol );
    syncLegacyResult();
    return;
  }

  if( isZero( rhs.norm2(), 0.0 ) )
  {
    sol.zero();
    m_result.numIterations = 0;
    m_result.residualReduction = 0.0;
    m_result.solveTime = 0.0;
    m_result.status = LinearSolverResult::Status::Success;
    return;
  }

  applyHypredrive( rhs, sol );

  int numIterations = 0;
  checkHypredriveCall( HYPREDRV_LinearSolverGetNumIter( m_hypredrive, &numIterations ),
                       "HYPREDRV_LinearSolverGetNumIter" );
  m_result.numIterations = numIterations;

  checkHypredriveCall( HYPREDRV_LinearSolverGetSolveTime( m_hypredrive, &m_result.solveTime ),
                       "HYPREDRV_LinearSolverGetSolveTime" );

  // Reuse the residual vector bound during setup: HypreVector::create is collective,
  // so allocating one per solve puts a collective in the solve path.
  matrix().residual( sol, rhs, m_residual );
  real64 const rhsNorm = rhs.norm2();
  real64 const denominator = rhsNorm > 0.0 ? rhsNorm : 1.0;
  m_result.residualReduction = m_residual.norm2() / denominator;
  m_result.status = ( m_result.residualReduction <= m_params.krylov.relTolerance )
                    ? LinearSolverResult::Status::Success
                    : LinearSolverResult::Status::NotConverged;

  if( m_params.logLevel >= 1 )
  {
    HYPRE_BigInt globalNumRows, globalNumNonzeros;
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixGetGlobalInfo( matrix().unwrappedIJ(),
                                                       &globalNumRows,
                                                       &globalNumRows,
                                                       &globalNumNonzeros ) );

    GEOS_LOG_RANK_0( GEOS_FMT( "        Linear Solver | {} | Unknowns: {} | Nonzeros: {} | Iterations: {} | Final Rel Res: {:.4e} | Setup Time: {:.3f} s | Solve Time: {:.3f} s",
                               m_result.status,
                               stringutilities::addCommaSeparators( globalNumRows ),
                               stringutilities::addCommaSeparators( globalNumNonzeros ),
                               m_result.numIterations,
                               m_result.residualReduction,
                               m_result.setupTime,
                               m_result.solveTime ) );
  }
}

void HypredriveSolver::syncLegacyResult() const
{
  if( m_legacySolver )
  {
    m_result = m_legacySolver->result();
  }
}

void HypredriveSolver::syncExecutionAnnotations()
{
  if( m_hypredrive == nullptr || !m_hasExecutionContext )
  {
    return;
  }

  std::string const desiredTimestepScope = makeTimestepScopeName( m_executionContext );
  std::string const desiredNewtonScope = makeNewtonScopeName( m_executionContext );

  if( m_newtonScopeActive && m_activeNewtonScope != desiredNewtonScope )
  {
    checkHypredriveCall( HYPREDRV_AnnotateLevelEnd( m_hypredrive, 1,
                                                    m_activeNewtonScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelEnd" );
    m_newtonScopeActive = false;
    m_activeNewtonScope.clear();
  }

  if( m_timestepScopeActive && m_activeTimestepScope != desiredTimestepScope )
  {
    checkHypredriveCall( HYPREDRV_AnnotateLevelEnd( m_hypredrive, 0,
                                                    m_activeTimestepScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelEnd" );
    m_timestepScopeActive = false;
    m_activeTimestepScope.clear();
  }

  if( !m_timestepScopeActive )
  {
    checkHypredriveCall( HYPREDRV_AnnotateLevelBegin( m_hypredrive, 0,
                                                      desiredTimestepScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelBegin" );
    m_timestepScopeActive = true;
    m_activeTimestepScope = desiredTimestepScope;
  }

  if( !m_newtonScopeActive )
  {
    checkHypredriveCall( HYPREDRV_AnnotateLevelBegin( m_hypredrive, 1,
                                                      desiredNewtonScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelBegin" );
    m_newtonScopeActive = true;
    m_activeNewtonScope = desiredNewtonScope;
  }
}

void HypredriveSolver::closeExecutionAnnotations()
{
  if( m_hypredrive != nullptr )
  {
    if( m_newtonScopeActive )
    {
      checkHypredriveCall( HYPREDRV_AnnotateLevelEnd( m_hypredrive, 1,
                                                      m_activeNewtonScope.c_str(), -1 ),
                           "HYPREDRV_AnnotateLevelEnd" );
    }

    if( m_timestepScopeActive )
    {
      checkHypredriveCall( HYPREDRV_AnnotateLevelEnd( m_hypredrive, 0,
                                                      m_activeTimestepScope.c_str(), -1 ),
                           "HYPREDRV_AnnotateLevelEnd" );
    }
  }

  m_newtonScopeActive = false;
  m_timestepScopeActive = false;
  m_activeNewtonScope.clear();
  m_activeTimestepScope.clear();
}

void HypredriveSolver::destroyHypredrive()
{
  if( m_hypredrive != nullptr )
  {
    if( m_linearSolverCreated )
    {
      checkHypredriveCall( HYPREDRV_LinearSolverDestroy( m_hypredrive ),
                           "HYPREDRV_LinearSolverDestroy" );
      m_linearSolverCreated = false;
    }
    checkHypredriveCall( HYPREDRV_Destroy( &m_hypredrive ), "HYPREDRV_Destroy" );
    m_hypredrive = nullptr;
  }
}

void HypredriveSolver::resetHypredriveState()
{
  closeExecutionAnnotations();
  destroyHypredrive();
  m_dummyRhs.reset();
  m_dummySol.reset();
  m_configurationSignature.clear();
  m_structureSignature.clear();
}

void HypredriveSolver::clear()
{
  Base::clear();

  resetHypredriveState();

  if( m_legacySolver )
  {
    m_legacySolver->clear();
    m_legacySolver.reset();
  }
}

}
