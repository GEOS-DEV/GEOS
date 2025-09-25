/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A / TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICVTISEM_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICVTISEM_HPP_

#include "AcousticSecondOrderAnisotropicWaveEquationSEM.hpp"


namespace geos
{

/**
 * @brief Intermediate base class for second-order VTI (Vertical Transverse Isotropy) solvers.
 *
 * This class serves as a common parent for all VTI-based solvers.
 * It does not add new functionality on its own but provides a clear inheritance.
 */
class AcousticSecondOrderAnisotropicVTISEM : public AcousticSecondOrderAnisotropicWaveEquationSEM
{
public:
  /**
   * @brief Inherit the constructor from the parent class.
   */
  using AcousticSecondOrderAnisotropicWaveEquationSEM::AcousticSecondOrderAnisotropicWaveEquationSEM;
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICVTISEM_HPP_
