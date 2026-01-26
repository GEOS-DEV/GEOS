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

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICTTISEM_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICTTISEM_HPP_

#include "AcousticSecondOrderAnisotropicWaveEquationSEM.hpp"


namespace geos
{

/**
 * @brief Intermediate base class for second-order TTI (Tilted Transverse Isotropy) solvers.
 *
 * Its primary responsibility is to register and compute the data fields required for
 * TTI anisotropy (Tilt, Azimuth, Dip, etc.), building upon the VTI foundation
 * provided by its parent class.
 */
class AcousticSecondOrderAnisotropicTTISEM : public AcousticSecondOrderAnisotropicWaveEquationSEM
{
public:
  // Inherit constructor from the parent class
  using AcousticSecondOrderAnisotropicWaveEquationSEM::AcousticSecondOrderAnisotropicWaveEquationSEM;

  // Override to add TTI-specific fields
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  // Override to zero out TTI-specific nodal fields
  virtual void postInputInitialization() override;

protected:
  virtual void initializePostInitialConditionsPreSubGroups() override final;
/**
 * @brief Template method for TTI matrix initialization (common VTI stuff + TTI DOF arrays)
 * @tparam DampingComputer The specific damping computer to use (Zhang or Fletcher)
 * @param mesh The mesh level
 * @param regionNames The region names to process
 */
  template< typename DampingComputer >
  void initializeMatricesTemplate( MeshLevel & mesh, string_array const & regionNames );

};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERANISOTROPICTTISEM_HPP_
