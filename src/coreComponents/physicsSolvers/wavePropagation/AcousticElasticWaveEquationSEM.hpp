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


/**
 * @file AcousticElasticWaveEquation.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_

#include "physicsSolvers/wavePropagation/CoupledWaveSolver.hpp"
#include "physicsSolvers/wavePropagation/ElasticWaveEquationSEM.hpp"
#include "physicsSolvers/wavePropagation/AcousticWaveEquationSEM.hpp"

namespace geos
{

class AcousticElasticWaveEquationSEM : public CoupledWaveSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >
{
public:
  using Base = CoupledWaveSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >;
  using Base::m_solvers;

  enum class SolverType : integer
  {
    AcousticWaveEquationSEM = 0,
    ElasticWaveEquationSEM = 1
  };

  /// String used to form the solverName used to register solvers in CoupledWaveSolver
  static string coupledSolverAttributePrefix() { return "acousticelastic"; }

  /**
   * @brief main constructor for AcousticElasticWaveEquationSEM objects
   * @param name the name of this instantiation of AcousticElasticWaveEquationSEM in the repository
   * @param parent the parent group of this instantiation of AcousticElasticWaveEquationSEM
   */
  AcousticElasticWaveEquationSEM( const string & name,
                                  Group * const parent )
    : Base( name, parent )
  { }

  /// Destructor for the class
  ~AcousticElasticWaveEquationSEM() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new AcousticElasticWaveEquationSEM object through the object catalog.
   */
  static string catalogName() { return "AcousticElasticSEM"; }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  AcousticWaveEquationSEM * acousticSolver() const
  {
    return std::get< toUnderlying( SolverType::AcousticWaveEquationSEM ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  ElasticWaveEquationSEM * elasticSolver() const
  {
    return std::get< toUnderlying( SolverType::ElasticWaveEquationSEM ) >( m_solvers );
  }

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  SortedArray< localIndex > m_interfaceNodesSet;
};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_ */
