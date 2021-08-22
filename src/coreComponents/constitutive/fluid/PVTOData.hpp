/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PVTOData.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTODATA_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTODATA_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

namespace constitutive
{

class PVTOData
{
public:

  class KernelWrapper
  {
public:

    /**
     * @brief The constructor of the kernel wrapper
     */
    KernelWrapper( arrayView1d< real64 const > const & Rs,
                   arrayView1d< real64 const > const & bubblePressure,
                   arrayView1d< real64 const > const & saturatedBo,
                   arrayView1d< real64 const > const & saturatedViscosity,
                   arrayView2d< real64 const > const & undersaturatedPressure,
                   arrayView2d< real64 const > const & undersaturatedBo,
                   arrayView2d< real64 const > const & undersaturatedViscosity,
                   arrayView1d< real64 const > const & surfaceMassDensity,
                   arrayView1d< real64 const > const & surfaceMoleDensity );

    /// Default constructor for the kernel wrapper
    KernelWrapper() = default;

    /// Default copy constructor
    KernelWrapper( KernelWrapper const & ) = default;

    /// Default move constructor
    KernelWrapper( KernelWrapper && ) = default;

    /// Copy assignment operator
    KernelWrapper & operator=( KernelWrapper const & ) = delete;

    /// Deleted move assignment operator
    KernelWrapper & operator=( KernelWrapper && ) = delete;

    /**
     * @brief The function that populates the wrapper
     */
    void create( arrayView1d< real64 const > const & Rs,
                 arrayView1d< real64 const > const & bubblePressure,
                 arrayView1d< real64 const > const & saturatedBo,
                 arrayView1d< real64 const > const & saturatedViscosity,
                 arrayView2d< real64 const > const & undersaturatedPressure,
                 arrayView2d< real64 const > const & undersaturatedBo,
                 arrayView2d< real64 const > const & undersaturatedViscosity,
                 arrayView1d< real64 const > const & surfaceMassDensity,
                 arrayView1d< real64 const > const & surfaceMoleDensity );

    void move( LvArray::MemorySpace const space, bool const touch )
    {
      m_Rs.move( space, touch );
      m_bubblePressure.move( space, touch );

      m_saturatedBo.move( space, touch );
      m_saturatedViscosity.move( space, touch );

      m_undersaturatedPressure2d.move( space, touch );
      m_undersaturatedBo2d.move( space, touch );
      m_undersaturatedViscosity2d.move( space, touch );

      m_surfaceMassDensity.move( space, touch );
      m_surfaceMoleDensity.move( space, touch );
    }

public:

    /// Rs array
    arrayView1d< real64 const > m_Rs;
    /// Bubble-point pressure array
    arrayView1d< real64 const > m_bubblePressure;

    // Saturated data (free gas phase present)

    /// Saturated oil phase formation volume factor
    arrayView1d< real64 const > m_saturatedBo;
    /// Saturated oil phase viscosity
    arrayView1d< real64 const > m_saturatedViscosity;

    // Undersaturated data (no free gas)

    /// Undersaturated pressure
    arrayView2d< real64 const > m_undersaturatedPressure2d;
    /// Undersaturated oil phase formation volume factor
    arrayView2d< real64 const > m_undersaturatedBo2d;
    /// Undersaturated oil phase viscosity
    arrayView2d< real64 const > m_undersaturatedViscosity2d;

    /// Surface mass density
    arrayView1d< real64 const > m_surfaceMassDensity;
    /// Surface mole density
    arrayView1d< real64 const > m_surfaceMoleDensity;

  };

public:

  PVTOData()
    :
    m_maxRelativePressure( 0.0 ),
    m_nSaturatedPoints( 0 )
  {}

  KernelWrapper createKernelWrapper() const;

  /// Max pressure used for the construction of the undersaturated tables
  real64 m_maxRelativePressure;

  /// Rs array
  array1d< real64 > m_Rs;
  /// Bubble-point pressure array
  array1d< real64 > m_bubblePressure;

  // Saturated data (free gas phase present)

  /// Number of saturated points
  localIndex m_nSaturatedPoints;
  /// Saturated oil phase formation volume factor
  array1d< real64 > m_saturatedBo;
  /// Saturated oil phase viscosity
  array1d< real64 > m_saturatedViscosity;

  // Undersaturated data (no free gas)

  /// Undersaturated pressure
  array2d< real64 > m_undersaturatedPressure2d;
  /// Undersaturated oil phase formation volume factor
  array2d< real64 > m_undersaturatedBo2d;
  /// Undersaturated oil phase viscosity
  array2d< real64 > m_undersaturatedViscosity2d;

  /// Surface mass density
  array1d< real64 > m_surfaceMassDensity;
  /// Surface mole density
  array1d< real64 > m_surfaceMoleDensity;

  // temporary arrays used in the construction of the undersaturated arrays
  // can be discarded after construction
  array1d< array1d< real64 > > m_undersaturatedPressure;
  array1d< array1d< real64 > > m_undersaturatedBo;
  array1d< array1d< real64 > > m_undersaturatedViscosity;

  KernelWrapper m_kernelWrapper;

};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTODATA_HPP_
