/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file SolverStatistics.hpp
 */


#ifndef GEOSX_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
#define GEOSX_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP

#include "dataRepository/Group.hpp"

namespace geosx
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
   */
  SolverStatistics( string const & name,
                    dataRepository::Group * const parent );

  /**
   * @brief Default destructor
   */
  virtual ~SolverStatistics() = default;

  /**
   * @brief Deleted default constructor.
   */
  SolverStatistics() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  SolverStatistics( SolverStatistics const & ) = delete;

  /**
   * @brief Default move Constructor
   * @param The source object of the move.
   */
  SolverStatistics( SolverStatistics && ) = default;

  /**
   * @brief Deleted assignment operator.
   */
  SolverStatistics & operator=( SolverStatistics const & ) = delete;

  /**
   * @brief Deleted move operator.
   */
  SolverStatistics & operator=( SolverStatistics && ) = delete;

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
   * @brief Tell the solverStatistics that there is a time step cut
   */
  void logTimeStepCut();

  /**
   * @brief Save the statistics for the individual time step and increment the cumulative stats
   */
  void saveTimeStepStatistics();

  /**
   * @brief Output the cumulative statistics to the terminal
   */
  void outputStatistics() const;

private:

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the number of time steps
    static constexpr char const * numTimeStepsString() { return "numTimeSteps"; }
    /// String key for the number of time steps
    static constexpr char const * numTimeStepCutsString() { return "numTimeStepCuts"; }

    /// String key for the successful number of nonlinear iterations
    static constexpr char const * numSuccessfulNonlinearIterationsString() { return "numSuccessfulNonlinearIterations"; }
    /// String key for the successful number of linear iterations
    static constexpr char const * numSuccessfulLinearIterationsString() { return "numSuccessfulLinearIterations"; }
    /// String key for the failed number of nonlinear iterations
    static constexpr char const * numFailedNonlinearIterationsString() { return "numFailedNonlinearIterations"; }
    /// String key for the failed number of linear iterations
    static constexpr char const * numFailedLinearIterationsString() { return "numFailedLinearIterations"; }
  };

  /// Number of time steps
  integer m_numTimeSteps;

  /// Number of time step cuts
  integer m_numTimeStepCuts;

  /// Number of nonlinear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumNonlinearIterations;

  /// Number of linear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumLinearIterations;

  /// Cumulative number of successful nonlinear iterations
  integer m_numSuccessfulNonlinearIterations;

  /// Cumulative number of successful linear iterations
  integer m_numSuccessfulLinearIterations;

  /// Cumulative number of failed nonlinear iterations
  integer m_numFailedNonlinearIterations;

  /// Cumulative number of failed linear iterations
  integer m_numFailedLinearIterations;

};

} //namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
