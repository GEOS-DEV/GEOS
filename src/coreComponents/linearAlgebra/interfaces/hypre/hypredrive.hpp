#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPRE_HYPREDRIVE_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPRE_HYPREDRIVE_HPP_

#include "common/LinearSolverBase.hpp"
#include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"

#include <HYPREDRV.h>

#include <memory>
#include <string>

namespace geos
{

class HypreSolver;
class HypreDriveSolverTestPeer;

namespace hypre
{

namespace hypreDrive
{

enum class InputSource
{
  authoritativeFile,
  generatedFallback
};

struct InputArgsParseTarget
{
  InputSource source = InputSource::generatedFallback;
  std::string argument;
};

bool shouldUse( LinearSolverParameters const & params );

bool buildInputArgsParseTarget( LinearSolverParameters const & params,
                                InputArgsParseTarget & target );

bool buildInputArgsParseTarget( LinearSolverParameters const & params,
                                stdVector< string > const & fieldNames,
                                arrayView1d< int const > const & numComponentsPerField,
                                InputArgsParseTarget & target );

std::string formatInputArgsParseTargetYaml( InputArgsParseTarget const & target );

bool wasInputArgsParseTargetLogged( InputArgsParseTarget const & target );

void markInputArgsParseTargetLogged( InputArgsParseTarget const & target );

void logInputArgsParseTarget( LinearSolverParameters const & params,
                              InputArgsParseTarget const & target );

void initializeRuntime();

void finalizeRuntime();

}

}

class HypreDriveSolver final : public LinearSolverBase< HypreInterface >
{
public:

  using Base = LinearSolverBase< HypreInterface >;

  explicit HypreDriveSolver( LinearSolverParameters parameters );

  ~HypreDriveSolver() override;

  void setExecutionContext( LinearSolverExecutionContext const & context ) override;

  void setup( HypreMatrix const & mat ) override;

  void apply( HypreVector const & src,
              HypreVector & dst ) const override;

  void solve( HypreVector const & rhs,
              HypreVector & sol ) const override;

  void clear() override;

private:

  bool configureHypreDrive( HypreMatrix const & mat );

  void createHypreDrive( HypreMatrix const & mat,
                         hypre::hypreDrive::InputArgsParseTarget const & parseTarget,
                         std::string const & configurationSignature,
                         std::string const & structureSignature,
                         arrayView1d< int > const & pointMarkers );

  void refreshBoundObjects( HypreMatrix const & mat,
                            arrayView1d< int > const & pointMarkers );

  void setupLegacy( HypreMatrix const & mat );

  void applyHypreDrive( HypreVector const & rhs,
                        HypreVector & sol ) const;

  void syncExecutionAnnotations();

  void closeExecutionAnnotations();

  void syncLegacyResult() const;

  void destroyHypreDrive();

  void resetHypreDriveState();

  using Base::m_params;
  using Base::m_result;

  friend class HypreDriveSolverTestPeer;

  LinearSolverExecutionContext m_executionContext{};
  std::string m_configurationSignature;
  std::string m_structureSignature;
  std::string m_activeTimestepScope;
  std::string m_activeNewtonScope;
  bool m_hasExecutionContext = false;
  bool m_timestepScopeActive = false;
  bool m_newtonScopeActive = false;
  size_t m_hypreDriveGeneration = 0;
  HYPREDRV_t m_hypreDrive{};
  mutable HypreVector m_dummyRhs;
  mutable HypreVector m_dummySol;
  std::unique_ptr< HypreSolver > m_legacySolver;
};

}

#endif
