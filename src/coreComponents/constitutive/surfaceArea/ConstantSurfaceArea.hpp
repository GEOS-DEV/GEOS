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

/**
 * @file ConstantSurfaceArea.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SURFACEAREA_CONSTANTSURFACEAREA_HPP_
#define GEOS_CONSTITUTIVE_SURFACEAREA_CONSTANTSURFACEAREA_HPP_

#include "constitutive/surfaceArea/SurfaceAreaBase.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Update class keeping the reactive surface area at its initial value.
 */
class ConstantSurfaceAreaUpdates : public SurfaceAreaBaseUpdates
{
public:

  ConstantSurfaceAreaUpdates( arrayView3d< real64, reactivefluid::USD_SPECIES > const & surfaceArea,
                              arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialSurfaceArea ):
    SurfaceAreaBaseUpdates( surfaceArea, initialSurfaceArea )
  {}

  GEOS_HOST_DEVICE
  virtual void updateFromPorosityAndVolumeFractions( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & porosity,
                                                     real64 const & initialPorosity,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & volumeFractions,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & initialVolumeFractions ) const override final
  {
    GEOS_UNUSED_VAR( porosity, initialPorosity, volumeFractions, initialVolumeFractions );

    for( integer r = 0; r < numKineticReactions(); ++r )
    {
      m_surfaceArea[k][q][r] = m_initialSurfaceArea[k][q][r];
    }
  }
};


/**
 * @brief Model in which the reactive surface area of each mineral does not evolve.
 */
class ConstantSurfaceArea : public SurfaceAreaBase
{
public:

  ConstantSurfaceArea( string const & name, Group * const parent );

  static string catalogName() { return "ConstantSurfaceArea"; }

  virtual string getCatalogName() const override { return catalogName(); }

  using KernelWrapper = ConstantSurfaceAreaUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_surfaceArea, m_initialSurfaceArea );
  }
};

}/* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_SURFACEAREA_CONSTANTSURFACEAREA_HPP_
