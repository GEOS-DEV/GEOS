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

#include "physicsSolvers/wavePropagation/AcousticWaveEquationSEM.hpp"
#include "physicsSolvers/wavePropagation/ElasticWaveEquationSEM.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"

namespace geos
{

class AcousticElasticWaveEquationSEM public CoupledSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >
{
public:
  using Base = CoupledSolver< AcousticWaveEquationSEM, ElasticWaveEquationSEM >;
  using Base::m_solvers;
};

}

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEM_HPP_ */
