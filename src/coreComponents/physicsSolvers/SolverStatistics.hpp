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

/*
 * @file SolverStatistics.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
#define GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP

#include "dataRepository/Group.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"

namespace geos
{

/**
 * @brief Class containing iteration data for a time-step
 */
class IterationsStatistics : public dataRepository::Group
{
public:

/**
 * @brief Constructor for SolverStatistics Objects.
 * @param[in] name the name of this instantiation of SolverStatistics in the repository
 * @param[in] parent the parent group of this instantiation of SolverStatistics
 * @note Currently, we register only the iteration statistics as the convergence value will be lost during the solving
 */
  IterationsStatistics( string const & name,
                        dataRepository::Group * const parent );

  /// Number of time steps
  integer m_numTimeSteps = 0;

  /// Number of time step cuts
  integer m_numTimeStepCuts = 0;


  /// Number of outer loop iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumOuterLoopIterations = 0;

  /// Number of nonlinear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumNonlinearIterations = 0;

  /// Number of linear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumLinearIterations = 0;


  /// Cumulative number of successful outer loop iterations
  integer m_numSuccessfulOuterLoopIterations = 0;

  /// Cumulative number of successful nonlinear iterations
  integer m_numSuccessfulNonlinearIterations = 0;

  /// Cumulative number of successful linear iterations
  integer m_numSuccessfulLinearIterations = 0;


  /// Cumulative number of discarded outer loop iterations
  integer m_numDiscardedOuterLoopIterations = 0;

  /// Cumulative number of discarded nonlinear iterations
  integer m_numDiscardedNonlinearIterations = 0;

  /// Cumulative number of discarded linear iterations
  integer m_numDiscardedLinearIterations = 0;

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the number of time steps
    static constexpr char const * numTimeStepsString() { return "numTimeSteps"; }
    /// String key for the number of time step cuts
    static constexpr char const * numTimeStepCutsString() { return "numTimeStepCuts"; }

    /// String key for the successful number of outer loop iterations
    static constexpr char const * numSuccessfulOuterLoopIterationsString() { return "numSuccessfulOuterLoopIterations"; }
    /// String key for the successful number of nonlinear iterations
    static constexpr char const * numSuccessfulNonlinearIterationsString() { return "numSuccessfulNonlinearIterations"; }
    /// String key for the successful number of linear iterations
    static constexpr char const * numSuccessfulLinearIterationsString() { return "numSuccessfulLinearIterations"; }

    /// String key for the discarded number of outer loop iterations
    static constexpr char const * numDiscardedOuterLoopIterationsString() { return "numDiscardedOuterLoopIterations"; }
    /// String key for the discarded number of nonlinear iterations
    static constexpr char const * numDiscardedNonlinearIterationsString() { return "numDiscardedNonlinearIterations"; }
    /// String key for the discarded number of linear iterations
    static constexpr char const * numDiscardedLinearIterationsString() { return "numDiscardedLinearIterations"; }
  };
/**
 * @
 *
 */
  /**
   * @brief Initialize the counters used for an individual time step
   */
  virtual void resetCurrentTimeStepStatistics();

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @param[in] numLinearIterations the number of linear iterations done by the linear solver
   * @detail This function is well suited for Newton's method, or for single-physics solvers in sequential schemes
   */
  virtual void updateNonlinearIteration( integer const numLinearIterations );

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @detail This function is well suited for the outer loop in sequential schemes
   */
  virtual void updateNonlinearIteration();

  /**
   * @brief Tell the solverStatistics that we are doing an outer loop iteration
   */
  virtual void incrementNonlinearIteration();

  /**
   * @brief Tell the solverStatistics that there is a time step cut
   */
  virtual void updateTimeStepCut();

  /**
   * @brief Save the statistics for the individual time step and increment the cumulative stats
   */
  virtual void iterateTimeStepStatistics();

  /**
   * @brief Register the corresponding solver statistics to the TableData
   */
  virtual void writeStatsToTable();

  /**
   * @brief  Set the filename output file.
   * @param filename The filename as a string_view.
   */
  virtual void setFilename( string_view filename )
  { m_iterationsFilename = filename; }

  /**
   * @brief Output the statistics to the console and csv file if needed
   * @param writeCSV Indicate if we output to CSV FILE
   */
  virtual void outputStatistics( bool writeCSV );

  /**
   * @brief Set the Residual Norms filename
   * @param solverName The solverName as a string_view.
   */
  virtual void setSolverName( string_view solverName )
  {
    m_solverName = solverName;
  }

private:
  ///
  string m_solverName;
  /// Table containing statistics relative to non linear parameter
  TableData m_iterationData;
  /// Table Layout contenaning header for both CSV and log
  std::unique_ptr< TableLayout > m_iterationCSVLayout;
  /// Format the CSV iterations Data
  std::unique_ptr< TableCSVFormatter > m_iterationCSVFormatter;
  /// Filename for the convergence CSV.
  string m_iterationsFilename;
  // Ouput stream for each timestep
  std::ofstream logStream;
};

/**
 * @brief An empty class used for all sub-solvers instances
 */
class NullIterationsStatistics : public IterationsStatistics
{
public:

  NullIterationsStatistics( string const & name,
                            dataRepository::Group * const parent ): IterationsStatistics( name, parent ){}

  /**
   * @brief Group key associated with NullIterationsStatistics.
   */
  struct groupKeyStruct
  {
    /// @return string for the NullIterationsStatistics wrapper
    static constexpr char const * NullIterationsStatisticsString() { return "NullIterationsStatistics"; }
  };
///@{
/// @cond DO_NOT_DOCUMENT
  void resetCurrentTimeStepStatistics() override {}

  void updateNonlinearIteration( integer const numLinearIterations )override {}

  void updateNonlinearIteration() override {}

  void incrementNonlinearIteration() override {}

  void updateTimeStepCut() override {}

  void iterateTimeStepStatistics() override {}

  void writeStatsToTable() override {}

  void outputStatistics( bool writeCSV ) {}

  void setSolverName( string_view solverName ){}
  /// @endcond
///@}
};

/**
 * @brief Class containing convergence information given a time-step
 */
class ConvergenceStatistics
{
public:

  /**
   * @brief Construct a new Convergence Statistics object
   */
  ConvergenceStatistics();

  /// Current cycle number
  integer m_cycleNumber = 0;

  /// Number of time steps
  integer m_numTimeSteps = 0;

  /// Maximum number of current Newton iterations.
  integer m_currentNewtonIter = 0;

  /// Maximum value for residual mass.
  real64 m_residualMass = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual volume.
  real64 m_residualVol = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual energy.
  real64 m_residualEnergy = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual flow.
  real64 m_residualFlow = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual bubble displacement.
  real64 m_residualBubbleDisp = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual fracture.
  real64 m_residualFracture = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual stick.
  real64 m_residualStick = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual slip.
  real64 m_residualSlip = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual open.
  real64 m_residualOpen = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual solid.
  real64 m_residualSolid = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual contact.
  real64 m_residualContact = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual proppant.
  real64 m_residualProppant = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual well.
  real64 m_residualWell = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual damage.
  real64 m_residualDamage = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum total residual value.
  real64 m_totalResidual = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual norm at the end of line search
  real64 m_residualNormT = std::numeric_limits< real64 >::quiet_NaN();

  /**
   * @brief Output the cumulative statistics to the terminal
   * @param writeCSV Indicates if the output should be written to a CSV file.
   */
  void outputResidualNorm( bool writeCSV );

  /**
   * @brief Prepare the layout and register the corresponding residuals norms to the TableData
   */
  void writeResidualNormToTable();

  /**
   * @brief Remove the last residual norms when a configuration did not converge.
   * @note This is done based on the number of Newton iterations.
   */
  void removeInvalidResidualNorms();

  void updateCycleNumber( integer const cycleNumber )
  { m_cycleNumber = cycleNumber; }

  /**
   * @brief Save the current newton iteration
   * @param currentNewtonIter The current newton iteration performed by the the linear solver
   */
  void updateNewtonIter( integer currentNewtonIter );

  /**
   * @brief  Set the filename output file.
   * @param filename The filename as a string_view.
   */
  void setFilename( string_view filename )
  { m_convergenceFilename = filename; }

private:
  /// Table containing statistics related  to non linear norms
  std::unique_ptr< TableLayout > m_convergenceLayout;
  /// Table containing statistics data for nonlinear norms.
  TableData m_convergenceData;
  /// Format the CSV convergence Data
  std::unique_ptr< TableCSVFormatter > m_convergenceFormatter;
  /// Filename for the convergence CSV.
  string m_convergenceFilename;
  // Ouput stream for each timestep
  std::ofstream logStream;
};

/**
 * @class SolverStatistics
 * @brief This class records solver statistics for each time step.
 */
class SolverStatistics : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for SolverStatistics Objects.
   * @param[in] name the name of this instantiation of SolverStatistics in the repository
   * @param[in] parent the parent group of this instantiation of SolverStatistics
   * @noteCurrently we register only the iteration statistics as the convergence value will be lost during the solving
   */
  SolverStatistics( string const & name,
                    dataRepository::Group * const parent );

  /**
   * @brief Group key associated with IterationsStatistics.
   */
  struct groupKeyStruct
  {
    /// @return string for the IterationsStatistics wrapper
    static constexpr char const * IterationsStatisticsString() { return "IterationsStatistics"; }
  };

  /**
   * @return The output directory where all statistics related to the solver are atored
   */
  string_view getOutputDir() const
  { return m_outputDir; }

  /**
   * @brief Set the Residual Norms filename
   * @param solverName The solverName as a string_view.
   */
  void setOutputFilesName( string_view solverName );

  /// Contain iteration data given a time step
  IterationsStatistics m_iterationsStats;
  /// Contain convergence data given a time step
  ConvergenceStatistics m_convergenceStats;

private:
  /// Name of the directory containing solvers statistics csv
  string m_directoryName = "convergence";

  /// Output directory for solver statistics (CSV), passed in the constructor.
  string m_outputDir;

  std::ofstream logStream;
};



} //namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
