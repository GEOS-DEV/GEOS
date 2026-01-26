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

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERVTIFLETCHER_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERVTIFLETCHER_HPP_

#include "AcousticSecondOrderAnisotropicVTISEM.hpp"


namespace geos
{

/**
 * @brief Concrete solver for the second-order VTI wave equation using Fletcher's formulation.
 */
class AcousticSecondOrderVTIFletcher : public AcousticSecondOrderAnisotropicVTISEM
{
public:
  AcousticSecondOrderVTIFletcher( const std::string & name, Group * const parent );
  virtual ~AcousticSecondOrderVTIFletcher() override = default;

  static std::string catalogName() { return "AcousticSecondOrderVTIFletcherSEM"; }
  std::string getCatalogName() const override { return catalogName(); }

  // Override to add fields unique to the Fletcher method
  void registerDataOnMesh( Group & meshBodies ) override;

protected:
  // Implement the pure virtual function with the specific kernel for this solver
  void applyStiffnessKernels( const real64 & dt, MeshLevel & mesh, const string_array & regionNames ) override;
  virtual void initializeMatrices( MeshLevel & mesh, string_array const & regionNames ) override;
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICSECONDORDERVTIFLETCHER_HPP_
