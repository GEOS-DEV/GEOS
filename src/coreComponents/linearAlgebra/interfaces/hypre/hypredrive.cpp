#include "linearAlgebra/interfaces/hypre/hypredrive.hpp"

#include "common/GeosxConfig.hpp"
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

void checkHypreDriveCall( uint32_t const errorCode,
                          char const * const call )
{
  if( errorCode != 0 )
  {
    HYPREDRV_ErrorCodeDescribe( errorCode );
    GEOS_ERROR( GEOS_FMT( "Error in call to {}", call ) );
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

std::string makeInputArgsParseTargetSignature( hypre::hypreDrive::InputArgsParseTarget const & target )
{
  return GEOS_FMT( "{}\n{}",
                   target.source == hypre::hypreDrive::InputSource::authoritativeFile ? "authoritativeFile" : "generatedFallback",
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
      // hypreDrive follows hypre's MGR global smoother table, where code 2 is exposed as mixed-gs.
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

std::string joinLabelNames( std::vector< HYPRE_Int > const & values,
                            stdVector< string > const & labelNames )
{
  std::ostringstream stream;
  for( size_t i = 0; i < values.size(); ++i )
  {
    HYPRE_Int const label = values[i];
    GEOS_ERROR_IF( label < 0 || label >= LvArray::integerConversion< HYPRE_Int >( labelNames.size() ),
                   GEOS_FMT( "Invalid dof label {} while generating hypreDrive MGR YAML", label ) );
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
                 GEOS_FMT( "Unsupported generated hypreDrive solver {}", params.solverType ) );

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
  if( params.dofsPerNode > 1 )
  {
    appendLine( stream, 2, "reordering: 0" );
  }

  yaml = stream.str();
  return true;
}

enum class AMGFlavor
{
  pressure,
  pressureTemperature,
  displacementFiltered,
  displacement
};

struct LevelAMGBlock
{
  HYPRE_Int level;
  AMGFlavor flavor;
};

struct MGRSpecialization
{
  AMGFlavor coarseFlavor = AMGFlavor::pressure;
  std::vector< LevelAMGBlock > fRelaxAMGLevels;
  HYPRE_Int pmax = 0;
  HYPRE_Int coarseMinCoarseSize = -1;
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
      specialization.fRelaxAMGLevels.push_back( { 0, AMGFlavor::displacementFiltered } );
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
      specialization.fRelaxAMGLevels.push_back( { 1, AMGFlavor::displacement } );
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
}

void appendDisplacementAMG( std::ostringstream & stream,
                            integer const indentLevel,
                            integer const separateComponents,
                            bool const filterFunctions )
{
  appendAMGHeader( stream, indentLevel );
  appendLine( stream, indentLevel + 1, "coarsening:" );
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
  appendLine( stream, indentLevel + 2, "order: 1" );
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
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  appendLine( stream, indentLevel + 1, "coarsening:" );
  if( minCoarseSize >= 0 )
  {
    appendLine( stream, indentLevel + 2, GEOS_FMT( "min_coarse_size: {}", minCoarseSize ) );
  }
  appendLine( stream, indentLevel + 2, "type: 8" );
  appendLine( stream, indentLevel + 2, "max_row_sum: 1.0" );
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "down_type: 18" );
  appendLine( stream, indentLevel + 2, "up_type: 18" );
  appendLine( stream, indentLevel + 2, "coarse_type: 18" );
  appendLine( stream, indentLevel + 2, "num_sweeps: 2" );
#else
  if( minCoarseSize >= 0 )
  {
    appendLine( stream, indentLevel + 1, "coarsening:" );
    appendLine( stream, indentLevel + 2, GEOS_FMT( "min_coarse_size: {}", minCoarseSize ) );
  }
  appendLine( stream, indentLevel + 1, "relaxation:" );
  appendLine( stream, indentLevel + 2, "order: 1" );
#endif
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
  appendLine( stream, 0, "preconditioner:" );
  appendLine( stream, 1, "mgr:" );
  appendLine( stream, 2, "tolerance: 0.0" );
  appendLine( stream, 2, "max_iter: 1" );
  appendLine( stream, 2, GEOS_FMT( "print_level: {}", getMGRPrintLevel( params.logLevel ) ) );
  appendLine( stream, 2, GEOS_FMT( "non_c_to_f: {}", 1 ) );
  appendLine( stream, 2, GEOS_FMT( "nonglk_max_elmts: {}", 1 ) );
  appendLine( stream, 2, GEOS_FMT( "pmax: {}", specialization.pmax ) );
  appendLine( stream, 2, GEOS_FMT( "coarse_th: {}", strategy.m_coarseGridThreshold ) );
  appendLine( stream, 2, GEOS_FMT( "num_levels: {}", MGRStrategyProbe< STRATEGY >::numReductionLevels + 1 ) );
  appendLine( stream, 2, "level:" );

  std::vector< HYPRE_Int > activeLabels( static_cast< size_t >( strategy.m_numBlocks ) );
  std::iota( activeLabels.begin(), activeLabels.end(), 0 );

  for( HYPRE_Int level = 0; level < MGRStrategyProbe< STRATEGY >::numReductionLevels; ++level )
  {
    std::vector< HYPRE_Int > cLabels( strategy.m_labels[level].begin(), strategy.m_labels[level].end() );
    std::vector< HYPRE_Int > fLabels;

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
                     GEOS_FMT( "Missing hypreDrive AMG block description for MGR strategy {} level {}",
                               params.mgr.strategy, level ) );
      appendLine( stream, 4, "f_relaxation:" );
      appendAMGByFlavor( stream, 5, *amgFlavor, params.mgr );
    }
    else
    {
      char const * const fRelaxationName = getGeneratedMGRFRelaxationName( strategy.m_levelFRelaxType[level] );
      if( fRelaxationName == nullptr )
      {
        return false;
      }
      appendLine( stream, 4, GEOS_FMT( "f_relaxation: {}", fRelaxationName ) );
    }

    char const * const gRelaxationName = getGeneratedMGRGlobalSmootherName( strategy.m_levelGlobalSmootherType[level] );
    if( gRelaxationName == nullptr )
    {
      return false;
    }
    appendLine( stream, 4, GEOS_FMT( "g_relaxation: {}", gRelaxationName ) );
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
                                         hypre::hypreDrive::InputArgsParseTarget & target )
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

  target.source = hypre::hypreDrive::InputSource::generatedFallback;
  target.argument = stream.str();
  return true;
}

}

namespace hypre
{

namespace hypreDrive
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
      return GEOS_FMT( "# hypreDrive authoritative file: {}\n{}",
                       target.argument,
                       *contents );
    }

    return GEOS_FMT( "# hypreDrive authoritative file: {}\n"
                     "# GEOS was unable to read this YAML file while formatting the hypreDrive input dump.",
                     target.argument );
  }

  return GEOS_FMT( "# hypreDrive input generated by GEOS\n{}",
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
      GEOS_LOG_RANK_0( GEOS_FMT( "        hypreDrive input | authoritative file | {}\n{}",
                                 target.argument,
                                 *contents ) );
    }
    else
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "        hypreDrive input | authoritative file | {}\n"
                                 "        hypreDrive input | unable to read authoritative YAML from GEOS",
                                 target.argument ) );
    }
  }
  else
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "        hypreDrive input | generated fallback\n{}",
                               target.argument ) );
  }
}

void initializeRuntime()
{
  checkHypreDriveCall( HYPREDRV_Initialize(), "HYPREDRV_Initialize" );
}

void finalizeRuntime()
{
  checkHypreDriveCall( HYPREDRV_Finalize(), "HYPREDRV_Finalize" );
}

}

}

namespace
{

std::string makeConfigurationSignature( hypre::hypreDrive::InputArgsParseTarget const & target )
{
  return GEOS_FMT( "{}\n{}",
                   target.source == hypre::hypreDrive::InputSource::authoritativeFile ? "authoritativeFile" : "generatedFallback",
                   hypre::hypreDrive::formatInputArgsParseTargetYaml( target ) );
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

HypreDriveSolver::HypreDriveSolver( LinearSolverParameters parameters )
  : Base( std::move( parameters ) )
{}

HypreDriveSolver::~HypreDriveSolver()
{
  HypreDriveSolver::clear();
}

void HypreDriveSolver::setExecutionContext( LinearSolverExecutionContext const & context )
{
  m_executionContext = context;
  m_hasExecutionContext = true;
}

void HypreDriveSolver::setup( HypreMatrix const & mat )
{
  Base::setup( mat );

  if( !configureHypreDrive( mat ) )
  {
    resetHypreDriveState();
    setupLegacy( mat );
  }
}

void HypreDriveSolver::createHypreDrive( HypreMatrix const & mat,
                                         hypre::hypreDrive::InputArgsParseTarget const & parseTarget,
                                         std::string const & configurationSignature,
                                         std::string const & structureSignature,
                                         arrayView1d< int > const & pointMarkers )
{
  resetHypreDriveState();

  checkHypreDriveCall( HYPREDRV_Create( mat.comm(), &m_hypreDrive ), "HYPREDRV_Create" );
  checkHypreDriveCall( HYPREDRV_SetLibraryMode( m_hypreDrive ), "HYPREDRV_SetLibraryMode" );

  char * argv[] = { const_cast< char * >( parseTarget.argument.c_str() ) };
  checkHypreDriveCall( HYPREDRV_InputArgsParse( 1, argv, m_hypreDrive ),
                       "HYPREDRV_InputArgsParse" );
  if( parseTarget.source == hypre::hypreDrive::InputSource::generatedFallback &&
      m_hasExecutionContext &&
      !m_executionContext.solverName.empty() )
  {
    checkHypreDriveCall( HYPREDRV_ObjectSetName( m_hypreDrive,
                                                 m_executionContext.solverName.c_str() ),
                         "HYPREDRV_ObjectSetName" );
  }

  ++m_hypreDriveGeneration;
  m_configurationSignature = configurationSignature;
  m_structureSignature = structureSignature;

  refreshBoundObjects( mat, pointMarkers );
}

void HypreDriveSolver::refreshBoundObjects( HypreMatrix const & mat,
                                            arrayView1d< int > const & pointMarkers )
{
  checkHypreDriveCall( HYPREDRV_LinearSystemSetMatrix( m_hypreDrive,
                                                       toHypreMatrix( mat.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetMatrix" );

  if( !m_dummyRhs.ready() || m_dummyRhs.localSize() != mat.numLocalRows() )
  {
    m_dummyRhs.reset();
    m_dummySol.reset();
    m_dummyRhs.create( mat.numLocalRows(), mat.comm() );
    m_dummySol.create( mat.numLocalRows(), mat.comm() );
  }

  checkHypreDriveCall( HYPREDRV_LinearSystemSetRHS( m_hypreDrive,
                                                    reinterpret_cast< HYPRE_Vector >( m_dummyRhs.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetRHS" );
  checkHypreDriveCall( HYPREDRV_LinearSystemSetSolution( m_hypreDrive,
                                                         reinterpret_cast< HYPRE_Vector >( m_dummySol.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetSolution" );

  if( pointMarkers.size() > 0 )
  {
    checkHypreDriveCall( HYPREDRV_LinearSystemSetDofmap( m_hypreDrive,
                                                         LvArray::integerConversion< int >( pointMarkers.size() ),
                                                         pointMarkers.data() ),
                         "HYPREDRV_LinearSystemSetDofmap" );
  }
}

void HypreDriveSolver::setupLegacy( HypreMatrix const & mat )
{
  if( !m_legacySolver )
  {
    m_legacySolver = std::make_unique< HypreSolver >( m_params );
  }

  m_legacySolver->setup( mat );
  syncLegacyResult();
}

bool HypreDriveSolver::configureHypreDrive( HypreMatrix const & mat )
{
  if( !hypre::hypreDrive::shouldUse( m_params ) )
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
    dofManager->getLocalDofComponentLabels( pointMarkers );
  }

  hypre::hypreDrive::InputArgsParseTarget parseTarget;
  if( !hypre::hypreDrive::buildInputArgsParseTarget( m_params,
                                                     fieldNames,
                                                     numComponentsPerField,
                                                     parseTarget ) )
  {
    return false;
  }

  std::string const configurationSignature = makeConfigurationSignature( parseTarget );
  std::string const structureSignature =
    makeStructureSignature( mat,
                            fieldNames,
                            numComponentsPerField,
                            pointMarkers,
                            m_hasExecutionContext ? &m_executionContext : nullptr );

  hypre::hypreDrive::logInputArgsParseTarget( m_params, parseTarget );

  if( m_legacySolver )
  {
    m_legacySolver->clear();
    m_legacySolver.reset();
  }

  bool const recreateHandle = ( m_hypreDrive == nullptr ) ||
                              ( m_configurationSignature != configurationSignature ) ||
                              ( m_structureSignature != structureSignature );

  if( recreateHandle )
  {
    createHypreDrive( mat,
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

  checkHypreDriveCall( HYPREDRV_AnnotateBegin( m_hypreDrive, "system", -1 ),
                       "HYPREDRV_AnnotateBegin" );
  checkHypreDriveCall( HYPREDRV_LinearSolverCreate( m_hypreDrive ), "HYPREDRV_LinearSolverCreate" );
  checkHypreDriveCall( HYPREDRV_LinearSolverSetup( m_hypreDrive ), "HYPREDRV_LinearSolverSetup" );
  checkHypreDriveCall( HYPREDRV_AnnotateEnd( m_hypreDrive, "system", -1 ),
                       "HYPREDRV_AnnotateEnd" );
  checkHypreDriveCall( HYPREDRV_LinearSolverGetSetupTime( m_hypreDrive, &m_result.setupTime ),
                       "HYPREDRV_LinearSolverGetSetupTime" );

  return true;
}

void HypreDriveSolver::apply( HypreVector const & src,
                              HypreVector & dst ) const
{
  if( m_legacySolver )
  {
    m_legacySolver->apply( src, dst );
    return;
  }

  applyHypreDrive( src, dst );
}

void HypreDriveSolver::applyHypreDrive( HypreVector const & rhs,
                                        HypreVector & sol ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( rhs.ready() );
  GEOS_LAI_ASSERT( sol.ready() );

  m_dummySol.copy( sol );
  checkHypreDriveCall( HYPREDRV_LinearSystemSetInitialGuess( m_hypreDrive,
                                                             reinterpret_cast< HYPRE_Vector >( m_dummySol.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetInitialGuess" );
  checkHypreDriveCall( HYPREDRV_LinearSystemSetRHS( m_hypreDrive,
                                                    reinterpret_cast< HYPRE_Vector >( rhs.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetRHS" );
  checkHypreDriveCall( HYPREDRV_LinearSystemSetSolution( m_hypreDrive,
                                                         reinterpret_cast< HYPRE_Vector >( sol.unwrappedIJ() ) ),
                       "HYPREDRV_LinearSystemSetSolution" );
  checkHypreDriveCall( HYPREDRV_LinearSystemResetInitialGuess( m_hypreDrive ),
                       "HYPREDRV_LinearSystemResetInitialGuess" );
  checkHypreDriveCall( HYPREDRV_LinearSolverApply( m_hypreDrive ), "HYPREDRV_LinearSolverApply" );

  sol.touch();
}

void HypreDriveSolver::solve( HypreVector const & rhs,
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

  applyHypreDrive( rhs, sol );

  int numIterations = 0;
  checkHypreDriveCall( HYPREDRV_LinearSolverGetNumIter( m_hypreDrive, &numIterations ),
                       "HYPREDRV_LinearSolverGetNumIter" );
  m_result.numIterations = numIterations;

  checkHypreDriveCall( HYPREDRV_LinearSolverGetSolveTime( m_hypreDrive, &m_result.solveTime ),
                       "HYPREDRV_LinearSolverGetSolveTime" );

  HypreVector residual( rhs );
  matrix().residual( sol, rhs, residual );
  real64 const rhsNorm = rhs.norm2();
  real64 const denominator = rhsNorm > 0.0 ? rhsNorm : 1.0;
  m_result.residualReduction = residual.norm2() / denominator;
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

void HypreDriveSolver::syncLegacyResult() const
{
  if( m_legacySolver )
  {
    m_result = m_legacySolver->result();
  }
}

void HypreDriveSolver::syncExecutionAnnotations()
{
  if( m_hypreDrive == nullptr || !m_hasExecutionContext )
  {
    return;
  }

  std::string const desiredTimestepScope = makeTimestepScopeName( m_executionContext );
  std::string const desiredNewtonScope = makeNewtonScopeName( m_executionContext );

  if( m_newtonScopeActive && m_activeNewtonScope != desiredNewtonScope )
  {
    checkHypreDriveCall( HYPREDRV_AnnotateLevelEnd( m_hypreDrive, 1,
                                                    m_activeNewtonScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelEnd" );
    m_newtonScopeActive = false;
    m_activeNewtonScope.clear();
  }

  if( m_timestepScopeActive && m_activeTimestepScope != desiredTimestepScope )
  {
    checkHypreDriveCall( HYPREDRV_AnnotateLevelEnd( m_hypreDrive, 0,
                                                    m_activeTimestepScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelEnd" );
    m_timestepScopeActive = false;
    m_activeTimestepScope.clear();
  }

  if( !m_timestepScopeActive )
  {
    checkHypreDriveCall( HYPREDRV_AnnotateLevelBegin( m_hypreDrive, 0,
                                                      desiredTimestepScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelBegin" );
    m_timestepScopeActive = true;
    m_activeTimestepScope = desiredTimestepScope;
  }

  if( !m_newtonScopeActive )
  {
    checkHypreDriveCall( HYPREDRV_AnnotateLevelBegin( m_hypreDrive, 1,
                                                      desiredNewtonScope.c_str(), -1 ),
                         "HYPREDRV_AnnotateLevelBegin" );
    m_newtonScopeActive = true;
    m_activeNewtonScope = desiredNewtonScope;
  }
}

void HypreDriveSolver::closeExecutionAnnotations()
{
  if( m_hypreDrive != nullptr )
  {
    if( m_newtonScopeActive )
    {
      checkHypreDriveCall( HYPREDRV_AnnotateLevelEnd( m_hypreDrive, 1,
                                                      m_activeNewtonScope.c_str(), -1 ),
                           "HYPREDRV_AnnotateLevelEnd" );
    }

    if( m_timestepScopeActive )
    {
      checkHypreDriveCall( HYPREDRV_AnnotateLevelEnd( m_hypreDrive, 0,
                                                      m_activeTimestepScope.c_str(), -1 ),
                           "HYPREDRV_AnnotateLevelEnd" );
    }
  }

  m_newtonScopeActive = false;
  m_timestepScopeActive = false;
  m_activeNewtonScope.clear();
  m_activeTimestepScope.clear();
}

void HypreDriveSolver::destroyHypreDrive()
{
  if( m_hypreDrive != nullptr )
  {
    checkHypreDriveCall( HYPREDRV_Destroy( &m_hypreDrive ), "HYPREDRV_Destroy" );
    m_hypreDrive = nullptr;
  }
}

void HypreDriveSolver::resetHypreDriveState()
{
  closeExecutionAnnotations();
  destroyHypreDrive();
  m_dummyRhs.reset();
  m_dummySol.reset();
  m_configurationSignature.clear();
  m_structureSignature.clear();
}

void HypreDriveSolver::clear()
{
  Base::clear();

  resetHypreDriveState();

  if( m_legacySolver )
  {
    m_legacySolver->clear();
    m_legacySolver.reset();
  }
}

}
