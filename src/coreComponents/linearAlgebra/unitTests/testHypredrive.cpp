#include "common/GeosxConfig.hpp"
#include <sstream>
#include "common/format/Format.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#ifdef GEOS_USE_HYPRE
#include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"
#include "linearAlgebra/interfaces/hypre/mgrStrategies/SinglePhasePoromechanicsConformingFracturesALM.hpp"
#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"

#ifdef GEOS_USE_HYPREDRV
#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/interfaces/hypre/hypredrive.hpp"
#endif
#endif

namespace geos
{

#if defined(GEOS_USE_HYPRE) && defined(GEOS_USE_HYPREDRV)
class HypredriveSolverTestPeer
{
public:
  static size_t generation( HypredriveSolver const & solver )
  {
    return solver.m_hypredriveGeneration;
  }

  static HYPREDRV_t handle( HypredriveSolver const & solver )
  {
    return solver.m_hypredrive;
  }
};

namespace
{

stdVector< string > makeFieldNames()
{
  stdVector< string > fieldNames;
  for( integer i = 0; i < 10; ++i )
  {
    fieldNames.emplace_back( GEOS_FMT( "field{}", i ) );
  }
  return fieldNames;
}

array1d< int > makeMgrNumComponentsPerField()
{
  array1d< int > result( 10 );
  for( localIndex i = 0; i < result.size(); ++i )
  {
    // Use a generous synthetic field size so every MGR strategy constructor
    // has enough labels to build its reduction hierarchy during the probe.
    result[i] = 3;
  }
  return result;
}

LinearSolverParameters makeMgrParameters( LinearSolverParameters::MGR::StrategyType const strategy )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::mgr;
  params.mgr.strategy = strategy;
  params.mgr.separateComponents = 1;
  params.mgr.areWellsShut = 0;
  return params;
}

class ScopedEnvVar final
{
public:
  ScopedEnvVar( char const * const name,
                char const * const value )
    : m_name( name )
  {
    char const * const previous = std::getenv( name );
    if( previous != nullptr )
    {
      m_hadPrevious = true;
      m_previous = previous;
    }

    setenv( m_name.c_str(), value, 1 );
  }

  ~ScopedEnvVar()
  {
    if( m_hadPrevious )
    {
      setenv( m_name.c_str(), m_previous.c_str(), 1 );
    }
    else
    {
      unsetenv( m_name.c_str() );
    }
  }

private:
  std::string m_name;
  std::string m_previous;
  bool m_hadPrevious = false;
};

class ScopedCoutCapture final
{
public:
  ScopedCoutCapture()
    : m_previous( std::cout.rdbuf( m_stream.rdbuf() ) )
  {}

  ~ScopedCoutCapture()
  {
    std::cout.rdbuf( m_previous );
  }

  std::string str() const
  {
    return m_stream.str();
  }

private:
  std::ostringstream m_stream;
  std::streambuf * m_previous;
};

LinearSolverParameters makeReusableAMGParameters()
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.logLevel = 0;
  return params;
}

LinearSolverExecutionContext makeExecutionContext( Timestamp const systemSetupTimestamp,
                                                   integer const nonlinearIteration,
                                                   std::string const & solverName = {} )
{
  LinearSolverExecutionContext context;
  context.solverName = solverName;
  context.cycleNumber = 7;
  context.timeStepAttempt = 0;
  context.configurationAttempt = 0;
  context.nonlinearIteration = nonlinearIteration;
  context.systemSetupTimestamp = systemSetupTimestamp;
  return context;
}

TEST( HypredriveYaml, BuildsGeneratedFallbackForAMG )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_EQ( target.source, hypre::hypredrive::InputSource::generatedFallback );
  EXPECT_NE( target.argument.find( "solver:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "gmres:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "preconditioner:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "amg:" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "linear_system:" ), std::string::npos );
}

TEST( HypredriveYaml, AmgIluSmootherDisablesRcm )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::iluk;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_NE( target.argument.find( "ilu:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "reordering: 0" ), std::string::npos );
}

TEST( HypredriveYaml, MapsFGSAndBGSRelaxationNames )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.coarseType = LinearSolverParameters::AMG::CoarseType::default_;

  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::fgs;
  {
    hypre::hypredrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
    EXPECT_NE( target.argument.find( "down_type: forward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "up_type: forward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "coarse_type: forward-hl1gs" ), std::string::npos );
  }

  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::bgs;
  {
    hypre::hypredrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
    EXPECT_NE( target.argument.find( "down_type: backward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "up_type: backward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "coarse_type: backward-hl1gs" ), std::string::npos );
  }

  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::l1sgs;
  params.amg.coarseType = LinearSolverParameters::AMG::CoarseType::direct;
  {
    hypre::hypredrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
    EXPECT_NE( target.argument.find( "down_type: l1-hsgs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "up_type: l1-hsgs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "coarse_type: ge" ), std::string::npos );
    EXPECT_EQ( target.argument.find( "coarse_type: l1-hsgs" ), std::string::npos );
  }
}

TEST( HypredriveYaml, UsesCanonicalL1JacobiRelaxationName )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::l1jacobi;
  params.amg.coarseType = LinearSolverParameters::AMG::CoarseType::l1jacobi;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_NE( target.argument.find( "down_type: l1-jacobi" ), std::string::npos );
  EXPECT_NE( target.argument.find( "up_type: l1-jacobi" ), std::string::npos );
  EXPECT_NE( target.argument.find( "coarse_type: l1-jacobi" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "l1jacobi" ), std::string::npos );
}

TEST( HypredriveYaml, BuildsGeneratedYamlForEveryMGRStrategy )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  stdVector< string > const fieldNames = makeFieldNames();
  array1d< int > const numComponentsPerField = makeMgrNumComponentsPerField();

  for( integer value = static_cast< integer >( StrategyType::singlePhaseReservoirFVM );
       value <= static_cast< integer >( StrategyType::solidMechanicsEmbeddedFractures );
       ++value )
  {
    StrategyType const strategy = static_cast< StrategyType >( value );
    stdVector< string > strategyFieldNames = fieldNames;
    array1d< int > strategyNumComponentsPerField = numComponentsPerField;
    if( strategy == StrategyType::singlePhasePoromechanicsConformingFracturesALM )
    {
      strategyFieldNames = { "field0", "field1", "field2" };
      strategyNumComponentsPerField.resize( 3 );
    }
    else if( strategy == StrategyType::singlePhasePoromechanicsConformingFracturesALMReservoirFVM )
    {
      strategyFieldNames = { "field0", "field1", "field2", "field3" };
      strategyNumComponentsPerField.resize( 4 );
    }
    hypre::hypredrive::InputArgsParseTarget target;

    try
    {
      EXPECT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( makeMgrParameters( strategy ),
                                                                 strategyFieldNames,
                                                                 strategyNumComponentsPerField,
                                                                 target ) )
        << static_cast< int >( value );
      EXPECT_EQ( target.source, hypre::hypredrive::InputSource::generatedFallback ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "linear_system:" ), std::string::npos ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "dof_labels:" ), std::string::npos ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "field0_0: 0" ), std::string::npos ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "preconditioner:" ), std::string::npos ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "mgr:" ), std::string::npos ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "coarsest_level:" ), std::string::npos ) << static_cast< int >( value );
      EXPECT_NE( target.argument.find( "f_dofs: [field" ), std::string::npos ) << static_cast< int >( value );
    }
    catch( std::exception const & error )
    {
      FAIL() << static_cast< int >( value ) << ": " << error.what();
    }
  }
}

TEST( HypredriveYaml, BuildsSelectedALMPoromechanicsMGRStrategy )
{
  stdVector< string > const fieldNames = { "totalDisplacement", "totalBubbleDisplacement", "pressure" };
  array1d< int > numComponentsPerField( 3 );
  numComponentsPerField[0] = 3;
  numComponentsPerField[1] = 3;
  numComponentsPerField[2] = 1;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget(
                 makeMgrParameters( LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFracturesALM ),
                 fieldNames,
                 numComponentsPerField,
                 target ) );
  // The outer MGR has one reduction level and its F-relaxation is a nested
  // two-level MGR for the displacement and bubble-displacement blocks.
  EXPECT_EQ( target.argument.find( "num_levels: 3" ), std::string::npos );
  std::string::size_type const outerNumLevels = target.argument.find( "num_levels: 2" );
  ASSERT_NE( outerNumLevels, std::string::npos );
  EXPECT_NE( target.argument.find( "num_levels: 2", outerNumLevels + 1 ), std::string::npos );
  EXPECT_NE( target.argument.find( "cycle: v(1,0)" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [totalDisplacement_0, totalDisplacement_1, totalDisplacement_2, "
                                   "totalBubbleDisplacement_0, totalBubbleDisplacement_1, totalBubbleDisplacement_2]" ),
             std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [totalDisplacement_0, totalDisplacement_1, totalDisplacement_2]" ),
             std::string::npos );
  EXPECT_NE( target.argument.find( "filter_functions: 1" ), std::string::npos );
  EXPECT_NE( target.argument.find( "aggressive:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "max_nnz_row: 10" ), std::string::npos );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  EXPECT_NE( target.argument.find( "down_type: 16" ), std::string::npos );
  EXPECT_NE( target.argument.find( "up_type: 16" ), std::string::npos );
#else
  EXPECT_NE( target.argument.find( "down_type: l1sym-hgs" ), std::string::npos );
  EXPECT_NE( target.argument.find( "up_type: l1sym-hgs" ), std::string::npos );
  EXPECT_NE( target.argument.find( "coarse_type: ge" ), std::string::npos );
  EXPECT_NE( target.argument.find( "order: 0" ), std::string::npos );
#endif
}

TEST( HypreMGR, SetsUpFullyCoupledSinglePhaseALM )
{
  array1d< int > numComponentsPerField( 3 );
  numComponentsPerField[0] = 3;
  numComponentsPerField[1] = 3;
  numComponentsPerField[2] = 1;

  HypreMatrix matrix;
  testing::computeIdentity( MPI_COMM_GEOS, 7, matrix );

  LinearSolverParameters params = makeMgrParameters( LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFracturesALM );
  HyprePrecWrapper precond;
  ASSERT_EQ( HYPRE_MGRCreate( &precond.ptr ), 0 );
  HypreMGRData mgrData;
  mgrData.pointMarkers.resize( 7 );
  for( HYPRE_Int label = 0; label < 7; ++label )
  {
    mgrData.pointMarkers[label] = label;
  }

  hypre::mgr::SinglePhasePoromechanicsConformingFracturesALM strategy( numComponentsPerField.toView() );
  strategy.setup( params.mgr, precond, mgrData );

  EXPECT_EQ( HYPRE_MGRSetup( precond.ptr, matrix.unwrapped(), nullptr, nullptr ), 0 );

  HypreVector rhs;
  HypreVector solution;
  rhs.create( matrix.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );
  solution.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  solution.zero();
  EXPECT_EQ( HYPRE_MGRSolve( precond.ptr, matrix.unwrapped(), rhs.unwrapped(), solution.unwrapped() ), 0 );

  EXPECT_EQ( HYPRE_MGRDestroy( precond.ptr ), 0 );
  EXPECT_EQ( mgrData.coarseSolver.destroy( mgrData.coarseSolver.ptr ), 0 );
  EXPECT_EQ( mgrData.nestedSolver.destroy( mgrData.nestedSolver.ptr ), 0 );
}

size_t countSubstrings( std::string const & text, std::string const & token )
{
  size_t count = 0;
  for( std::string::size_type pos = 0;
       ( pos = text.find( token, pos ) ) != std::string::npos;
       pos += token.size() )
  {
    ++count;
  }
  return count;
}

TEST( HypredriveYaml, NestedMgrAmgUsesHypreBoomerAMGCreateMaxCoarseSize )
{
  stdVector< string > const fieldNames = makeFieldNames();
  array1d< int > const numComponentsPerField = makeMgrNumComponentsPerField();

  {
    hypre::hypredrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget(
                   makeMgrParameters( LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics ),
                   fieldNames,
                   numComponentsPerField,
                   target ) );
    std::string::size_type const coarsest = target.argument.find( "coarsest_level:" );
    ASSERT_NE( coarsest, std::string::npos );
    EXPECT_NE( target.argument.substr( 0, coarsest ).find( "max_coarse_size: 9" ), std::string::npos );
    EXPECT_NE( target.argument.substr( coarsest ).find( "max_coarse_size: 9" ), std::string::npos );
    EXPECT_EQ( countSubstrings( target.argument, "max_coarse_size: 9" ), 2 );
    EXPECT_EQ( target.argument.find( "relax_type:" ), std::string::npos );
    EXPECT_EQ( countSubstrings( target.argument, "type: schwarz" ), 2 );
    EXPECT_EQ( countSubstrings( target.argument, "num_levels: 0" ), 2 );
    EXPECT_EQ( countSubstrings( target.argument, "reordering: 0" ), 2 );
  }

  {
    hypre::hypredrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget(
                   makeMgrParameters( LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanics ),
                   fieldNames,
                   numComponentsPerField,
                   target ) );
    std::string::size_type const coarsest = target.argument.find( "coarsest_level:" );
    ASSERT_NE( coarsest, std::string::npos );
    EXPECT_NE( target.argument.substr( 0, coarsest ).find( "max_coarse_size: 9" ), std::string::npos );
    EXPECT_NE( target.argument.substr( coarsest ).find( "max_coarse_size: 9" ), std::string::npos );
    EXPECT_EQ( countSubstrings( target.argument, "max_coarse_size: 9" ), 2 );
  }

  {
    hypre::hypredrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget(
                   makeMgrParameters( LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM ),
                   fieldNames,
                   numComponentsPerField,
                   target ) );
    std::string::size_type const coarsest = target.argument.find( "coarsest_level:" );
    ASSERT_NE( coarsest, std::string::npos );
    EXPECT_EQ( target.argument.substr( 0, coarsest ).find( "max_coarse_size: 9" ), std::string::npos );
    EXPECT_NE( target.argument.substr( coarsest ).find( "max_coarse_size: 9" ), std::string::npos );
    EXPECT_EQ( countSubstrings( target.argument, "max_coarse_size: 9" ), 1 );
    EXPECT_EQ( target.argument.find( "relax_type:" ), std::string::npos );
    EXPECT_EQ( countSubstrings( target.argument, "type: schwarz" ), 1 );
    EXPECT_EQ( countSubstrings( target.argument, "num_levels: 0" ), 1 );
    EXPECT_NE( target.argument.find( "reordering: 0" ), std::string::npos );
  }
}

TEST( HypredriveYaml, NestedMgrIluMatchesHypreMgrType16Defaults )
{
  stdVector< string > const fieldNames = makeFieldNames();
  array1d< int > const numComponentsPerField = makeMgrNumComponentsPerField();

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget(
                 makeMgrParameters( LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM ),
                 fieldNames,
                 numComponentsPerField,
                 target ) );

  std::string::size_type const gRelax = target.argument.find( "g_relaxation:" );
  ASSERT_NE( gRelax, std::string::npos );
  std::string const fromGRelax = target.argument.substr( gRelax );
  EXPECT_NE( fromGRelax.find( "reordering: 0" ), std::string::npos );
  EXPECT_NE( fromGRelax.find( "max_row_nnz: 1000" ), std::string::npos );
  EXPECT_NE( fromGRelax.find( "fill_level: 0" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "g_relaxation: ilu" ), std::string::npos );

  EXPECT_NE( target.argument.find( "type: \"none\"" ), std::string::npos );
  EXPECT_NE( target.argument.find( "num_sweeps: 0" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "relax_type:" ), std::string::npos );
}

TEST( HypredriveYaml, UsesSemanticCompositionalLabelsWhenFieldNamesAreAvailable )
{
  stdVector< string > const fieldNames = { "totalDisplacement", "compositionalVariables" };
  array1d< int > numComponentsPerField( 2 );
  numComponentsPerField[0] = 3;
  numComponentsPerField[1] = 3;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( makeMgrParameters( LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics ),
                                                             fieldNames,
                                                             numComponentsPerField,
                                                             target ) );
  EXPECT_NE( target.argument.find( "totalDisplacement_0: 0" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "totaldisplacement_0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "pressure: 3" ), std::string::npos );
  EXPECT_NE( target.argument.find( "density_0: 4" ), std::string::npos );
  EXPECT_NE( target.argument.find( "density_1: 5" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "compositionalvariables_0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [density_1]" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [density_0]" ), std::string::npos );
}

TEST( HypredriveYaml, UsesSemanticCompositionalWellLabelsWhenFieldNamesAreAvailable )
{
  stdVector< string > const fieldNames = { "totalDisplacement", "compositionalVariables", "compositionalWellVars" };
  array1d< int > numComponentsPerField( 3 );
  numComponentsPerField[0] = 3;
  numComponentsPerField[1] = 3;
  numComponentsPerField[2] = 4;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( makeMgrParameters( LinearSolverParameters::MGR::StrategyType::multiphasePoromechanicsReservoirFVM ),
                                                             fieldNames,
                                                             numComponentsPerField,
                                                             target ) );
  EXPECT_NE( target.argument.find( "totalDisplacement_0: 0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellPressure: 6" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellDensity_0: 7" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellDensity_1: 8" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellRate: 9" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "compositionalwellvars_0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [wellPressure, wellDensity_0, wellDensity_1, wellRate]" ), std::string::npos );
}

TEST( HypredriveYaml, UsesTemperatureLabelForThermalCompositionalBlocks )
{
  stdVector< string > const fieldNames = { "compositionalVariables" };
  array1d< int > numComponentsPerField( 1 );
  numComponentsPerField[0] = 4;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( makeMgrParameters( LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseFVM ),
                                                             fieldNames,
                                                             numComponentsPerField,
                                                             target ) );
  EXPECT_NE( target.argument.find( "pressure: 0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "density_0: 1" ), std::string::npos );
  EXPECT_NE( target.argument.find( "density_1: 2" ), std::string::npos );
  EXPECT_NE( target.argument.find( "temperature: 3" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [density_1]" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [density_0]" ), std::string::npos );
}

TEST( HypredriveYaml, UsesWellTemperatureLabelForThermalReservoirWellBlocks )
{
  stdVector< string > const fieldNames = { "compositionalVariables", "compositionalWellVars" };
  array1d< int > numComponentsPerField( 2 );
  numComponentsPerField[0] = 4;
  numComponentsPerField[1] = 5;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( makeMgrParameters( LinearSolverParameters::MGR::StrategyType::thermalCompositionalMultiphaseReservoirFVM ),
                                                             fieldNames,
                                                             numComponentsPerField,
                                                             target ) );
  EXPECT_NE( target.argument.find( "wellPressure: 4" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellDensity_0: 5" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellDensity_1: 6" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellRate: 7" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellTemperature: 8" ), std::string::npos );
}

TEST( HypredriveYaml, BuildsThermalSinglePhasePoromechanicsReservoirYaml )
{
  stdVector< string > const fieldNames = { "totalDisplacement", "singlePhaseVariables", "singlePhaseWellVars" };
  array1d< int > numComponentsPerField( 3 );
  numComponentsPerField[0] = 3;
  numComponentsPerField[1] = 2;
  numComponentsPerField[2] = 3;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( makeMgrParameters( LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanicsReservoirFVM ),
                                                             fieldNames,
                                                             numComponentsPerField,
                                                             target ) );
  EXPECT_NE( target.argument.find( "totalDisplacement_0: 0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "pressure: 3" ), std::string::npos );
  EXPECT_NE( target.argument.find( "temperature: 4" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellPressure: 5" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellRate: 6" ), std::string::npos );
  EXPECT_NE( target.argument.find( "wellTemperature: 7" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "singlePhaseVariables_0" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "singlePhaseWellVars_0" ), std::string::npos );
  EXPECT_NE( target.argument.find( "f_dofs: [wellPressure, wellRate, wellTemperature]" ), std::string::npos );
  EXPECT_NE( target.argument.find( "num_functions: 2" ), std::string::npos );
}

TEST( HypredriveYaml, UsesAuthoritativeFileWhenProvided )
{
  std::string const authoritativeFile = "/tmp/geos-hypredrive-authoritative.yml";

  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::direct;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::none;
  params.hypredriveInputFile = Path( authoritativeFile.c_str() );

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_EQ( target.source, hypre::hypredrive::InputSource::authoritativeFile );
  EXPECT_EQ( target.argument, authoritativeFile );
  EXPECT_TRUE( hypre::hypredrive::shouldUse( params ) );

  auto solver = HypreInterface::createSolver( params );
  EXPECT_NE( dynamic_cast< HypredriveSolver * >( solver.get() ), nullptr );
}

TEST( HypredriveLogging, LogLevelGatesGeneratedYamlDump )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );

  {
    ScopedCoutCapture capture;
    params.logLevel = 0;
    hypre::hypredrive::logInputArgsParseTarget( params, target );
    EXPECT_TRUE( capture.str().empty() );
  }

  {
    ScopedCoutCapture capture;
    params.logLevel = 1;
    hypre::hypredrive::logInputArgsParseTarget( params, target );
    EXPECT_NE( capture.str().find( "generated fallback" ), std::string::npos );
    EXPECT_NE( capture.str().find( "preconditioner:" ), std::string::npos );
  }
}

TEST( HypredriveLogging, LogsAuthoritativeFileContents )
{
  LinearSolverParameters params;
  params.hypredriveInputFile = "/tmp/geos-hypredrive-authoritative.yml";
  params.logLevel = 1;

  {
    std::ofstream output( params.hypredriveInputFile );
    ASSERT_TRUE( output.good() );
    output << "solver:\n"
              "  gmres:\n"
              "    max_iter: 3\n";
  }

  hypre::hypredrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypredrive::buildInputArgsParseTarget( params, target ) );

  ScopedCoutCapture capture;
  hypre::hypredrive::logInputArgsParseTarget( params, target );
  EXPECT_NE( capture.str().find( "authoritative file" ), std::string::npos );
  EXPECT_NE( capture.str().find( params.hypredriveInputFile ), std::string::npos );
  EXPECT_NE( capture.str().find( "max_iter: 3" ), std::string::npos );

  EXPECT_EQ( std::remove( params.hypredriveInputFile.c_str() ), 0 );
}

TEST( HypredriveLogging, PrintsStatisticsSummaryWhenHandleIsDestroyed )
{
  std::string const authoritativeFile = "/tmp/geos-hypredrive-statistics-authoritative.yml";
  {
    std::ofstream output( authoritativeFile );
    ASSERT_TRUE( output.good() );
    output << "general:\n"
              "  statistics: 1\n"
              "solver:\n"
              "  gmres:\n"
              "    max_iter: 5\n"
              "preconditioner:\n"
              "  amg:\n"
              "    print_level: 0\n";
  }

  LinearSolverParameters params;
  params.hypredriveInputFile = Path( authoritativeFile.c_str() );

  HypreMatrix matrix;
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix );

  HypreVector rhs;
  HypreVector sol;
  rhs.create( matrix.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );
  sol.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  sol.zero();

  HypredriveSolver solver( params );
  solver.setExecutionContext( makeExecutionContext( 51, 0 ) );
  solver.setup( matrix );
  solver.solve( rhs, sol );

  ::testing::internal::CaptureStdout();
  solver.clear();
  std::string const output = ::testing::internal::GetCapturedStdout();

  EXPECT_NE( output.find( "STATISTICS SUMMARY" ), std::string::npos );
  EXPECT_EQ( std::remove( authoritativeFile.c_str() ), 0 );
}

TEST( HypredriveLogging, GeneratedFallbackNamesStatisticsSummary )
{
  LinearSolverParameters params = makeReusableAMGParameters();

  HypreMatrix matrix;
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix );

  HypreVector rhs;
  HypreVector sol;
  rhs.create( matrix.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );
  sol.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  sol.zero();

  HypredriveSolver solver( params );
  solver.setExecutionContext( makeExecutionContext( 61, 0, "namedGeneratedSolver" ) );
  solver.setup( matrix );
  solver.solve( rhs, sol );

  ::testing::internal::CaptureStdout();
  solver.clear();
  std::string const output = ::testing::internal::GetCapturedStdout();

  EXPECT_NE( output.find( "STATISTICS SUMMARY for namedGeneratedSolver:" ), std::string::npos );
}

TEST( HypredriveSolverSelection, UsesAdapterWhenAvailable )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  auto solver = HypreInterface::createSolver( params );
  EXPECT_NE( dynamic_cast< HypredriveSolver * >( solver.get() ), nullptr );
}

TEST( HypredriveSolverSelection, HonorsLegacyOverride )
{
  ScopedEnvVar legacyOverride( "GEOS_HYPREDRV_FORCE_LEGACY", "1" );

  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  auto solver = HypreInterface::createSolver( params );
  EXPECT_NE( dynamic_cast< HypreSolver * >( solver.get() ), nullptr );
}

void compareHypredriveAndLegacySolutions( LinearSolverParameters const & params,
                                          int const numDofTags = 1 )
{
  struct ScopedKrylovDofLabels
  {
    ~ScopedKrylovDofLabels()
    {
      hypre::testing::clearKrylovDofLabels();
    }
  };

  ScopedKrylovDofLabels const clearLabelsOnExit;

  HypreMatrix matrix;
  testing::compute2DLaplaceOperator( MPI_COMM_GEOS, 100, matrix );

  array1d< int > labels;
  if( numDofTags > 1 )
  {
    labels.resize( matrix.numLocalRows() );
    for( localIndex i = 0; i < labels.size(); ++i )
    {
      labels[i] = LvArray::integerConversion< int >( i % numDofTags );
    }
    hypre::testing::setKrylovDofLabels( labels.toViewConst() );
  }

  HypreVector rhs;
  rhs.create( matrix.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );

  HypreVector solHypredrive;
  HypreVector solLegacy;
  solHypredrive.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  solLegacy.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  solHypredrive.zero();
  solLegacy.zero();

  HypredriveSolver hypredriveSolver( params );
  hypredriveSolver.setup( matrix );
  hypredriveSolver.solve( rhs, solHypredrive );
  ASSERT_TRUE( hypredriveSolver.result().success() );
  hypredriveSolver.clear();

  HypreSolver legacySolver( params );
  legacySolver.setup( matrix );
  legacySolver.solve( rhs, solLegacy );
  ASSERT_TRUE( legacySolver.result().success() );
  legacySolver.clear();

  EXPECT_EQ( hypredriveSolver.result().numIterations, legacySolver.result().numIterations );

  // Both solutions satisfy the same tolerance; their difference is bounded by
  // the solve tolerance amplified by the operator conditioning.
  HypreVector difference( solHypredrive );
  difference.axpy( -1.0, solLegacy );
  EXPECT_LE( difference.norm2() / solLegacy.norm2(), 1e-9 );

  // Both must satisfy the requested tolerance against the true residual.
  HypreVector residual( rhs );
  matrix.residual( solHypredrive, rhs, residual );
  EXPECT_LE( residual.norm2() / rhs.norm2(), 2.0 * params.krylov.relTolerance );
}

TEST( HypredriveNumerics, MatchesLegacyHypreSolverOnLaplaceGmresAmg )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.krylov.relTolerance = 1e-12;
  params.krylov.maxIterations = 300;
  params.logLevel = 0;

  compareHypredriveAndLegacySolutions( params );
}

TEST( HypredriveNumerics, MatchesLegacyHypreSolverOnLaplaceGmresAmgWithMultipleDofTags )
{
  // hypredrive library mode exercises a tagged dofmap. The legacy solver uses
  // untagged caller-owned vectors and should still produce the same solution.
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.krylov.relTolerance = 1e-12;
  params.krylov.maxIterations = 300;
  params.logLevel = 0;

  compareHypredriveAndLegacySolutions( params, 2 );
}

TEST( HypredriveNumerics, MatchesLegacyHypreSolverOnLaplaceBicgstabIlu )
{
  // Second solver/preconditioner combination supported by the generated-YAML path.
  // (The 2D elasticity test operator is assembled without Dirichlet conditions and is
  // singular, so equivalence is exercised on the Laplace operator only.)
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::bicgstab;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::iluk;
  params.krylov.relTolerance = 1e-12;
  params.krylov.maxIterations = 300;
  params.logLevel = 0;

  compareHypredriveAndLegacySolutions( params );
}

TEST( HypredriveLogging, WarnsOnceWhenGeneratedYamlFallsBackToLegacy )
{
  // AMG with rigid-body-mode near null space is not representable in generated
  // hypredrive YAML, so the adapter must warn once and fall back to the legacy
  // solver, which treats an absent null space benignly.
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.nullSpaceType = LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes;
  params.krylov.relTolerance = 1e-10;
  params.logLevel = 0;

  HypreMatrix matrix;
  testing::compute2DLaplaceOperator( MPI_COMM_GEOS, 20, matrix );

  HypreVector rhs;
  HypreVector sol;
  rhs.create( matrix.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );
  sol.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  sol.zero();

  HypredriveSolver solver( params );

  ScopedCoutCapture capture;
  solver.setup( matrix );
  solver.setup( matrix );
  std::string const output = capture.str();

  std::string const warning = "falling back to the legacy hypre solver";
  std::string::size_type const first = output.find( warning );
  EXPECT_NE( first, std::string::npos );
  EXPECT_EQ( output.find( warning, first + warning.size() ), std::string::npos );

  solver.solve( rhs, sol );
  EXPECT_TRUE( solver.result().success() );
  solver.clear();
}

TEST( HypredriveSolverReuse, ReusesHandleAcrossCompatibleSetupCycles )
{
  HypreMatrix matrix1;
  HypreMatrix matrix2;
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix1 );
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix2 );

  HypreVector rhs;
  HypreVector sol;
  rhs.create( matrix1.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );
  sol.create( matrix1.numLocalCols(), MPI_COMM_GEOS );
  sol.zero();

  HypredriveSolver solver( makeReusableAMGParameters() );
  solver.setExecutionContext( makeExecutionContext( 11, 0 ) );
  solver.setup( matrix1 );

  HYPREDRV_t const handle1 = HypredriveSolverTestPeer::handle( solver );
  size_t const generation1 = HypredriveSolverTestPeer::generation( solver );
  ASSERT_NE( handle1, nullptr );
  solver.solve( rhs, sol );

  sol.zero();
  solver.setExecutionContext( makeExecutionContext( 11, 1 ) );
  solver.setup( matrix2 );

  EXPECT_EQ( HypredriveSolverTestPeer::handle( solver ), handle1 );
  EXPECT_EQ( HypredriveSolverTestPeer::generation( solver ), generation1 );
  solver.solve( rhs, sol );

  solver.clear();
}

TEST( HypredriveSolverReuse, RecreatesHandleWhenStructureChanges )
{
  HypreMatrix matrix1;
  HypreMatrix matrix2;
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix1 );
  testing::computeIdentity( MPI_COMM_GEOS, 5, matrix2 );

  HypredriveSolver solver( makeReusableAMGParameters() );
  solver.setExecutionContext( makeExecutionContext( 11, 0 ) );
  solver.setup( matrix1 );
  size_t const generation1 = HypredriveSolverTestPeer::generation( solver );

  solver.setExecutionContext( makeExecutionContext( 12, 0 ) );
  solver.setup( matrix2 );

  EXPECT_GT( HypredriveSolverTestPeer::generation( solver ), generation1 );
  solver.clear();
}

TEST( HypredriveSolverReuse, RecreatesHandleWhenAuthoritativeYamlChanges )
{
  std::string const authoritativeFile = "/tmp/geos-hypredrive-reuse-authoritative.yml";
  {
    std::ofstream output( authoritativeFile );
    ASSERT_TRUE( output.good() );
    output << "solver:\n"
              "  gmres:\n"
              "    max_iter: 5\n"
              "preconditioner:\n"
              "  reuse:\n"
              "    enabled: yes\n"
              "    frequency: 1\n"
              "  amg:\n"
              "    print_level: 0\n";
  }

  LinearSolverParameters params;
  params.hypredriveInputFile = Path( authoritativeFile.c_str() );

  HypreMatrix matrix;
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix );

  HypredriveSolver solver( params );
  solver.setExecutionContext( makeExecutionContext( 21, 0 ) );
  solver.setup( matrix );
  size_t const generation1 = HypredriveSolverTestPeer::generation( solver );

  {
    std::ofstream output( authoritativeFile );
    ASSERT_TRUE( output.good() );
    output << "solver:\n"
              "  gmres:\n"
              "    max_iter: 7\n"
              "preconditioner:\n"
              "  reuse:\n"
              "    enabled: yes\n"
              "    frequency: 1\n"
              "  amg:\n"
              "    print_level: 0\n";
  }

  solver.setExecutionContext( makeExecutionContext( 21, 1 ) );
  solver.setup( matrix );
  EXPECT_GT( HypredriveSolverTestPeer::generation( solver ), generation1 );

  solver.clear();
  EXPECT_EQ( std::remove( authoritativeFile.c_str() ), 0 );
}

TEST( HypredriveSolverReuse, KeepsMultipleSolverHandlesIndependent )
{
  HypreMatrix matrix1;
  HypreMatrix matrix2;
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix1 );
  testing::computeIdentity( MPI_COMM_GEOS, 4, matrix2 );

  HypreVector rhs1;
  HypreVector rhs2;
  HypreVector sol1;
  HypreVector sol2;
  rhs1.create( matrix1.numLocalRows(), MPI_COMM_GEOS );
  rhs2.create( matrix2.numLocalRows(), MPI_COMM_GEOS );
  rhs1.set( 1.0 );
  rhs2.set( 2.0 );
  sol1.create( matrix1.numLocalCols(), MPI_COMM_GEOS );
  sol2.create( matrix2.numLocalCols(), MPI_COMM_GEOS );
  sol1.zero();
  sol2.zero();

  HypredriveSolver solver1( makeReusableAMGParameters() );
  HypredriveSolver solver2( makeReusableAMGParameters() );

  solver1.setExecutionContext( makeExecutionContext( 31, 0 ) );
  solver1.setup( matrix1 );
  HYPREDRV_t const handle1 = HypredriveSolverTestPeer::handle( solver1 );
  size_t const generation1 = HypredriveSolverTestPeer::generation( solver1 );

  solver2.setExecutionContext( makeExecutionContext( 41, 0 ) );
  solver2.setup( matrix2 );
  HYPREDRV_t const handle2 = HypredriveSolverTestPeer::handle( solver2 );
  size_t const generation2 = HypredriveSolverTestPeer::generation( solver2 );

  ASSERT_NE( handle1, nullptr );
  ASSERT_NE( handle2, nullptr );
  EXPECT_NE( handle1, handle2 );

  solver1.solve( rhs1, sol1 );
  solver2.solve( rhs2, sol2 );

  sol1.zero();
  sol2.zero();

  solver1.setExecutionContext( makeExecutionContext( 31, 1 ) );
  solver1.setup( matrix1 );
  EXPECT_EQ( HypredriveSolverTestPeer::handle( solver1 ), handle1 );
  EXPECT_EQ( HypredriveSolverTestPeer::generation( solver1 ), generation1 );

  solver2.setExecutionContext( makeExecutionContext( 41, 1 ) );
  solver2.setup( matrix2 );
  EXPECT_EQ( HypredriveSolverTestPeer::handle( solver2 ), handle2 );
  EXPECT_EQ( HypredriveSolverTestPeer::generation( solver2 ), generation2 );

  solver1.setExecutionContext( makeExecutionContext( 32, 0 ) );
  solver1.setup( matrix1 );
  EXPECT_GT( HypredriveSolverTestPeer::generation( solver1 ), generation1 );

  size_t const refreshedGeneration2 = HypredriveSolverTestPeer::generation( solver2 );
  solver2.setExecutionContext( makeExecutionContext( 41, 2 ) );
  solver2.setup( matrix2 );
  EXPECT_EQ( HypredriveSolverTestPeer::handle( solver2 ), handle2 );
  EXPECT_EQ( HypredriveSolverTestPeer::generation( solver2 ), refreshedGeneration2 );

  solver1.clear();
  solver2.clear();
}

}

#endif

}

#if defined(GEOS_USE_HYPRE) && defined(GEOS_USE_HYPREDRV)
int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
#endif
