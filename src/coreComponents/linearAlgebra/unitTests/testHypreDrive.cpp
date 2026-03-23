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

#ifdef GEOS_USE_HYPREDRV
#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/hypredrive.hpp"
#endif
#endif

namespace geos
{

#if defined(GEOS_USE_HYPRE) && defined(GEOS_USE_HYPREDRV)
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

TEST( HypreDriveYaml, BuildsGeneratedFallbackForAMG )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  hypre::hypreDrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_EQ( target.source, hypre::hypreDrive::InputSource::generatedFallback );
  EXPECT_NE( target.argument.find( "solver:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "gmres:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "preconditioner:" ), std::string::npos );
  EXPECT_NE( target.argument.find( "amg:" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "linear_system:" ), std::string::npos );
}

TEST( HypreDriveYaml, MapsFGSAndBGSRelaxationNames )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.coarseType = LinearSolverParameters::AMG::CoarseType::default_;

  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::fgs;
  {
    hypre::hypreDrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );
    EXPECT_NE( target.argument.find( "down_type: forward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "up_type: forward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "coarse_type: forward-hl1gs" ), std::string::npos );
  }

  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::bgs;
  {
    hypre::hypreDrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );
    EXPECT_NE( target.argument.find( "down_type: backward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "up_type: backward-hl1gs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "coarse_type: backward-hl1gs" ), std::string::npos );
  }

  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::l1sgs;
  params.amg.coarseType = LinearSolverParameters::AMG::CoarseType::direct;
  {
    hypre::hypreDrive::InputArgsParseTarget target;
    ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );
    EXPECT_NE( target.argument.find( "down_type: l1-hsgs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "up_type: l1-hsgs" ), std::string::npos );
    EXPECT_NE( target.argument.find( "coarse_type: ge" ), std::string::npos );
    EXPECT_EQ( target.argument.find( "coarse_type: l1-hsgs" ), std::string::npos );
  }
}

TEST( HypreDriveYaml, UsesCanonicalL1JacobiRelaxationName )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  params.amg.smootherType = LinearSolverParameters::AMG::SmootherType::l1jacobi;
  params.amg.coarseType = LinearSolverParameters::AMG::CoarseType::l1jacobi;

  hypre::hypreDrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_NE( target.argument.find( "down_type: l1-jacobi" ), std::string::npos );
  EXPECT_NE( target.argument.find( "up_type: l1-jacobi" ), std::string::npos );
  EXPECT_NE( target.argument.find( "coarse_type: l1-jacobi" ), std::string::npos );
  EXPECT_EQ( target.argument.find( "l1jacobi" ), std::string::npos );
}

TEST( HypreDriveYaml, BuildsGeneratedYamlForEveryMGRStrategy )
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;

  stdVector< string > const fieldNames = makeFieldNames();
  array1d< int > const numComponentsPerField = makeMgrNumComponentsPerField();

  HypreInterface::initialize();

  for( integer value = static_cast< integer >( StrategyType::singlePhaseReservoirFVM );
       value <= static_cast< integer >( StrategyType::solidMechanicsEmbeddedFractures );
       ++value )
  {
    StrategyType const strategy = static_cast< StrategyType >( value );
    hypre::hypreDrive::InputArgsParseTarget target;

    try
    {
      EXPECT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( makeMgrParameters( strategy ),
                                                                 fieldNames,
                                                                 numComponentsPerField,
                                                                 target ) )
        << static_cast< int >( value );
      EXPECT_EQ( target.source, hypre::hypreDrive::InputSource::generatedFallback ) << static_cast< int >( value );
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

  HypreInterface::finalize();
}

TEST( HypreDriveYaml, UsesAuthoritativeFileWhenProvided )
{
  std::string const authoritativeFile = "/tmp/geos-hypredrive-authoritative.yml";

  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::direct;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::none;
  params.hypredriveInputFile = authoritativeFile;

  hypre::hypreDrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );
  EXPECT_EQ( target.source, hypre::hypreDrive::InputSource::authoritativeFile );
  EXPECT_EQ( target.argument, authoritativeFile );
  EXPECT_TRUE( hypre::hypreDrive::shouldUse( params ) );

  auto solver = HypreInterface::createSolver( params );
  EXPECT_NE( dynamic_cast< HypreDriveSolver * >( solver.get() ), nullptr );
}

TEST( HypreDriveLogging, LogLevelGatesGeneratedYamlDump )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  hypre::hypreDrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );

  {
    ScopedCoutCapture capture;
    params.logLevel = 0;
    hypre::hypreDrive::logInputArgsParseTarget( params, target );
    EXPECT_TRUE( capture.str().empty() );
  }

  {
    ScopedCoutCapture capture;
    params.logLevel = 1;
    hypre::hypreDrive::logInputArgsParseTarget( params, target );
    EXPECT_NE( capture.str().find( "generated fallback" ), std::string::npos );
    EXPECT_NE( capture.str().find( "preconditioner:" ), std::string::npos );
  }
}

TEST( HypreDriveLogging, LogsAuthoritativeFileContents )
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

  hypre::hypreDrive::InputArgsParseTarget target;
  ASSERT_TRUE( hypre::hypreDrive::buildInputArgsParseTarget( params, target ) );

  ScopedCoutCapture capture;
  hypre::hypreDrive::logInputArgsParseTarget( params, target );
  EXPECT_NE( capture.str().find( "authoritative file" ), std::string::npos );
  EXPECT_NE( capture.str().find( params.hypredriveInputFile ), std::string::npos );
  EXPECT_NE( capture.str().find( "max_iter: 3" ), std::string::npos );

  EXPECT_EQ( std::remove( params.hypredriveInputFile.c_str() ), 0 );
}

TEST( HypreDriveSolverSelection, UsesAdapterWhenAvailable )
{
  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  auto solver = HypreInterface::createSolver( params );
  EXPECT_NE( dynamic_cast< HypreDriveSolver * >( solver.get() ), nullptr );
}

TEST( HypreDriveSolverSelection, HonorsLegacyOverride )
{
  ScopedEnvVar legacyOverride( "GEOS_HYPREDRV_FORCE_LEGACY", "1" );

  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::gmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::amg;

  auto solver = HypreInterface::createSolver( params );
  EXPECT_NE( dynamic_cast< HypreSolver * >( solver.get() ), nullptr );
}

}
#endif

}
