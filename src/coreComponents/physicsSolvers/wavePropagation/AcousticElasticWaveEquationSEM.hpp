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

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/wavePropagation/AcousticWaveEquationSEM.hpp"
#include "physicsSolvers/wavePropagation/ElasticWaveEquationSEM.hpp"

namespace geos
{

class AcousticElasticWaveEquationSEM : public CoupledSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >
{
public:
  using Base = CoupledSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >;
  using Base::m_solvers;
  // using Base::m_dofManager;
  // using Base::m_localMatrix;
  // using Base::m_rhs;
  // using Base::m_solution;

  enum class SolverType : integer
  {
    AcousticWaveEquationSEM = 0,
    ElasticWaveEquationSEM = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "acousticelastic"; }

  /**
   * @brief main constructor for AcousticElasticWaveEquationSEM objects
   * @param name the name of this instantiation of AcousticElasticWaveEquationSEM in the repository
   * @param parent the parent group of this instantiation of AcousticElasticWaveEquationSEM
   */
  AcousticElasticWaveEquationSEM( const string & name,
                                  Group * const parent );

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
};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_ */
