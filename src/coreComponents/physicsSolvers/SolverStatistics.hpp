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
  /**
   * @brief Constructor for SolverStatistics Objects.
   * @param[in] name the name of this instantiation of SolverStatistics in the repository
   * @param[in] parent the parent group of this instantiation of SolverStatistics
   * @note For now, we register only the iteration statistics as the convergence value will be lost during the solving
   */
  SolverStatistics( string const & name,
                    dataRepository::Group * const parent );

  /**
   * @return The output directory where all statistics related to the solver are atored
   */
  string getOutputDir()
  { return m_outputDir; }

  class IterationsStatistics {} m_iterationsStats;

  class ConvergenceStatistics {} m_convergenceStats;

private:
  /// Output directory for solver statiscis csv
  string m_outputDir;

};

class SolverStatistics::IterationsStatistics : public SolverStatistics
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

  /**
   * @brief Save the statistics for the individual time step and increment the cumulative stats
   */
  void saveTimeStepStatistics();

  /**
   * @brief Register the corresponding solver statistics to the TableData
   */
  void registerStatsToTable();

  /**
   * @brief Output the statistics to the console and csv file if needed
   * @param writeCSV Indicate if we output to CSV FILE
   */
  void outputStatistics( bool writeCSV );


private:
  /// Table containing statistics relative to non linear parameter
  TableData m_nonLinearData;
};

class SolverStatistics::ConvergenceStatistics : public SolverStatistics
{
public:
  /// Maximum number of current Newton iterations.
  integer m_currentNewtonIter;

  /// Maximum value for residual mass.
  real64 m_residualMass;

  /// Maximum value for residual volume.
  real64 m_residualVol;

  /// Maximum value for residual energy.
  real64 m_residualEnergy;

  /// Maximum value for residual flow.
  real64 m_residualFlow;

  /// Maximum value for residual bubble displacement.
  real64 m_residualBubbleDisp;

  /// Maximum value for residual fracture.
  real64 m_residualFracture;

  /// Maximum value for residual stick.
  real64 m_residualStick;

  /// Maximum value for residual slip.
  real64 m_residualSlip;

  /// Maximum value for residual open.
  real64 m_residualOpen;

  /// Maximum value for residual solid.
  real64 m_residualSolid;

  /// Maximum value for residual contact.
  real64 m_residualContact;

  /// Maximum value for residual proppant.
  real64 m_residualProppant;

  /// Maximum value for residual well.
  real64 m_residualWell;

  /// Maximum value for residual damage.
  real64 m_residualDamage;

  /// Maximum total residual value.
  real64 m_totalResidual;

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the current number of Newton iterations.
    static constexpr char const * currentNewtonIterationsString() { return "currentNewtonIterations"; }
    /// String key for the residual mass.
    static constexpr char const * residualMassString() { return "residualMass"; }
    /// String key for the residual volume.
    static constexpr char const * residualVolumeString() { return "residualVolume"; }
    /// String key for the residual energy.
    static constexpr char const * residualEnergyString() { return "residualEnergy"; }
    /// String key for the residual flow.
    static constexpr char const * residualFlowString() { return "residualFlow"; }
    /// String key for the residual bubble displacement.
    static constexpr char const * residualBubbleDispString() { return "residualBubbleDisp"; }
    /// String key for the residual fracture.
    static constexpr char const * residualFractureString() { return "residualFracture"; }
    /// String key for the residual stick.
    static constexpr char const * residualStickString() { return "residualStick"; }
    /// String key for the residual slip.
    static constexpr char const * residualSlipString() { return "residualSlip"; }
    /// String key for the residual open.
    static constexpr char const * residualOpenString() { return "residualOpen"; }
    /// String key for the residual solid.
    static constexpr char const * residualSolidString() { return "residualSolid"; }
    /// String key for the residual contact.
    static constexpr char const * residualContactString() { return "residualContact"; }
    /// String key for the residual proppant.
    static constexpr char const * residualProppantString() { return "residualProppant"; }
    /// String key for the residual well.
    static constexpr char const * residualWellString() { return "residualWell"; }
    /// String key for the residual damage.
    static constexpr char const * residualDamageString() { return "residualDamage"; }
    /// String key for the total residual value.
    static constexpr char const * totalResidualString() { return "totalResidual"; }
  };

  /**
   * @brief Output the cumulative statistics to the terminal
   * @param Indicate if we output to CSV FILE
   */
  void outputResidualNorm( bool writeCSV );

  /**
   * @brief Set the Residual Norms filename
   * @param filename The filename string
   */
  void setResidualNormsFileName( string_view filename )
  { m_residualNormsFileName = filename; }

  /**
   * @brief Register the corresponding residuals norms to the TableData
   */
  void registerResidualNormToTable();

  /**
   * @brief When a configuration did not converge, removes the last residual norms by the number of
   * the newton iterations
   */
  void removeInvalidResidualNorms();

  /**
   * @brief Save the current newton iteration
   * @param currentNewtonIter The current newton iteration done by the the linear solver
   */
  void logNewtonIter( integer currentNewtonIter );

private:
  /// Table containing statistics relative to non linear norms
  std::unique_ptr< TableLayout > m_nonLinearNormsLayout;

  /// Table containing statistics relative to non linear norms
  TableData m_nonLinearNormsData;

  /// Residuals norms csv filename
  string m_residualNormsFileName;

};

} //namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
