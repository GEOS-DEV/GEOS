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
#include "common/format/table/TableLayout.hpp"

namespace geos
{

/**
 * @class SolverStatistics
 * @brief This class is used to log the solver statistics
 */
class SolverStatistics : public dataRepository::Group
{
public:

  /// Number of time steps
  integer m_numTimeSteps;

  /// Number of time step cuts
  integer m_numTimeStepCuts;


  /// Number of outer loop iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumOuterLoopIterations;

  /// Number of nonlinear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumNonlinearIterations;

  /// Number of linear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumLinearIterations;


  /// Cumulative number of successful outer loop iterations
  integer m_numSuccessfulOuterLoopIterations;

  /// Cumulative number of successful nonlinear iterations
  integer m_numSuccessfulNonlinearIterations;

  /// Cumulative number of successful linear iterations
  integer m_numSuccessfulLinearIterations;


  /// Cumulative number of discarded outer loop iterations
  integer m_numDiscardedOuterLoopIterations;

  /// Cumulative number of discarded nonlinear iterations
  integer m_numDiscardedNonlinearIterations;

  /// Cumulative number of discarded linear iterations
  integer m_numDiscardedLinearIterations;


  integer m_currentNewtonIter = std::numeric_limits< integer >::max();

  real64 m_residualMass = std::numeric_limits< real64 >::max();

  real64 m_residualVol = std::numeric_limits< real64 >::max();

  real64 m_residualEnergy =  std::numeric_limits< real64 >::max();

  real64 m_residualFlow =  std::numeric_limits< real64 >::max();

  real64 m_residualBubbleDisp =  std::numeric_limits< real64 >::max();

  real64 m_residualFracture =  std::numeric_limits< real64 >::max();

  real64 m_residualStick =  std::numeric_limits< real64 >::max();

  real64 m_residualSlip =  std::numeric_limits< real64 >::max();

  real64 m_residualOpen =  std::numeric_limits< real64 >::max();

  real64 m_residualSolid =  std::numeric_limits< real64 >::max();

  real64 m_residualContact =  std::numeric_limits< real64 >::max();

  real64 m_residualProppant =  std::numeric_limits< real64 >::max();

  real64 m_residualWell =  std::numeric_limits< real64 >::max();

  real64 m_residualDamage =  std::numeric_limits< real64 >::max();

  real64 m_totalResidual =  std::numeric_limits< real64 >::max();

  /**
   * @brief Constructor for SolverStatistics Objects.
   * @param[in] name the name of this instantiation of SolverStatistics in the repository
   * @param[in] parent the parent group of this instantiation of SolverStatistics
   */
  SolverStatistics( string const & name,
                    dataRepository::Group * const parent );

  /**
   * @brief Prepare the residual table layout whether it's thermal or not
   * @param isThermal Value indicating the thermal state
   */
  void prepareResidualTableLayout();

  /**
   * @brief Initialize the counters used for an individual time step
   */
  void initializeTimeStepStatistics();

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @param[in] numLinearIterations the number of linear iterations done by the linear solver
   * @detail This function is well suited for Newton's method, or for single-physics solvers in sequential schemes
   */
  void logNonlinearIteration( integer const numLinearIterations );

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @detail This function is well suited for the outer loop in sequential schemes
   */
  void logNonlinearIteration();

  /**
   * @brief Tell the solverStatistics that we are doing an outer loop iteration
   */
  void logOuterLoopIteration();

  /**
   * @brief Tell the solverStatistics that there is a time step cut
   */
  void logTimeStepCut();

  void removeInvalidResidualNorms();

  /**
   * @brief Save the current newton iteration
   * @param currentNewtonIter The current newton iteration done by the the linear solver
   */
  void logNewtonIter( integer currentNewtonIter );

  /**
   * @brief Save the current residuals norms to the corresponding TableData
   * @param residualNorm
   * @param residualMass
   * @param residualVol
   */

  /**
   * @brief Save the statistics for the individual time step and increment the cumulative stats
   */
  void saveTimeStepStatistics();

  /**
   * @brief Register the corresponding solver statistics to the TableData
   */
  void registerStatsToTable();

  /**
   * @brief Register the corresponding thermal residuals norms to the TableData
   */
  void registerResidualNormToTable();

  /**
   * @brief Output the cumulative statistics to the terminal
   */
  void outputStatistics( bool writeCSV );

  /**
   * @brief Output the cumulative statistics to the terminal
   */
  void outputResidualNorm( bool writeCSV );

  void setResidualNormsFileName( string_view filename )
  { m_residualNormsFileName = filename; }

  string getOutputDir()
  { return m_outputDir; }


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

private:
  /// Table containing statistics relative to non linear parameter
  TableData m_nonLinearData;

  /// Table containing statistics relative to non linear norms
  std::unique_ptr< TableLayout > m_nonLinearNormsLayout;

  /// Table containing statistics relative to non linear norms
  TableData m_nonLinearNormsData;

  /// Residuals norms csv filename
  string m_residualNormsFileName;

  /// Output directory for solver statiscis csv
  string m_outputDir;

};

} //namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
