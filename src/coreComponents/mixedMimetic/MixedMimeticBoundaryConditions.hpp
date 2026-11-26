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
 * @file MixedMimeticBoundaryConditions.hpp
 *
 * Boundary conditions of an operator of the mixed formulation. Each operator relates the normal flux f through a
 * face to the trace x_f of its primal variable x: M ( sigma f ) = x_K - x_f in the face equations, or
 * A f = x_L - x_R on a face condensed into a two-point flux relation. A boundary condition is a relation between
 * the trace and the normal flux, with sigma f the outward normal flux and |f| the face area:
 *
 *   Dirichlet   x_f = g                               g the prescribed trace
 *   Neumann     sigma f = g |f|                       g the prescribed outward normal flux per unit area
 *   Robin       sigma f = alpha |f| ( x_f - g )       alpha > 0 the transfer coefficient, g the exterior value
 *
 * The Robin condition contains the other two as limits: Dirichlet for alpha -> infinity and the homogeneous
 * Neumann condition for alpha -> 0. In the mixed formulation the Dirichlet and Robin conditions are natural and
 * the Neumann condition is essential: the flux degree of freedom is prescribed, sigma f - g |f| = 0. On a
 * non-condensed face the Robin trace is eliminated, x_f = g + sigma f / ( alpha |f| ), which adds one diagonal
 * term to the face equation; on a condensed face the resistance 1 / ( alpha |f| ) is added to A, with x_R = g.
 * The homogeneous Neumann condition applies to an impervious or adiabatic face. Data ( g, alpha ) that depend on the
 * time or on the solution history, as for a Carter-Tracy aquifer, are updated at every time step.
 */

#ifndef GEOS_MIXEDMIMETIC_MIXEDMIMETICBOUNDARYCONDITIONS_HPP_
#define GEOS_MIXEDMIMETIC_MIXEDMIMETICBOUNDARYCONDITIONS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace mixedMimeticBoundary
{

/// type of the condition on a face
struct BoundaryType
{
  static constexpr integer interior = 0;   ///< interior face, no condition
  static constexpr integer dirichlet = 1;  ///< prescribed trace
  static constexpr integer neumann = 2;    ///< prescribed outward normal flux per unit area
  static constexpr integer robin = 3;      ///< linear relation between the outward normal flux and the trace
};

/**
 * @struct FaceBoundaryCondition
 * @brief Boundary condition of one operator on one face.
 */
struct FaceBoundaryCondition
{
  integer type = BoundaryType::interior;   ///< BoundaryType
  real64 value = 0.0;                      ///< g: prescribed trace (Dirichlet), outward normal flux per unit area (Neumann), exterior value
                                           ///< (Robin)
  real64 coefficient = 0.0;                ///< transfer coefficient alpha of the Robin condition

  /// essential condition: the flux degree of freedom is prescribed
  GEOS_HOST_DEVICE
  bool isEssential() const { return type == BoundaryType::neumann; }

  /// natural condition: the face may be condensed into a two-point flux relation
  GEOS_HOST_DEVICE
  bool isCondensable() const { return type == BoundaryType::dirichlet || type == BoundaryType::robin; }

  /**
   * @brief Prescribed value of the flux degree of freedom of an essential face, from sigma f = g |f|.
   * @param sigma orientation of the face relative to the cell
   * @param area face area
   */
  GEOS_HOST_DEVICE
  real64 essentialFlux( real64 const sigma, real64 const area ) const
  {
    return sigma * value * area;
  }

  /**
   * @brief Trace in the equation of a non-condensed face: g (Dirichlet), g + sigma f / ( alpha |f| ) (Robin).
   * @param sigma orientation of the face relative to the cell
   * @param area face area
   * @param outwardFlux the outward normal flux sigma f
   * @param dTrace_dOutwardFlux the derivative of the trace with respect to sigma f
   */
  GEOS_HOST_DEVICE
  real64 trace( real64 const sigma, real64 const area, real64 const outwardFlux, real64 & dTrace_dOutwardFlux ) const
  {
    GEOS_UNUSED_VAR( sigma );
    dTrace_dOutwardFlux = ( type == BoundaryType::robin ) ? 1.0 / ( coefficient * area ) : 0.0;
    return value + dTrace_dOutwardFlux * outwardFlux;
  }

  /**
   * @brief Resistance added to the two-point flux relation of a condensed face: 1 / ( alpha |f| ) (Robin), 0 otherwise.
   * @param area face area
   */
  GEOS_HOST_DEVICE
  real64 resistance( real64 const area ) const
  {
    return ( type == BoundaryType::robin ) ? 1.0 / ( coefficient * area ) : 0.0;
  }
};

/**
 * @struct FaceBoundaryView
 * @brief Device view on the boundary conditions of one operator, stored per face.
 */
struct FaceBoundaryView
{
  arrayView1d< integer const > type;
  arrayView1d< real64 const > value;
  arrayView1d< real64 const > coefficient;

  GEOS_HOST_DEVICE
  FaceBoundaryCondition operator[]( localIndex const kf ) const
  {
    return { type[kf], value[kf], coefficient[kf] };
  }
};

} // namespace mixedMimeticBoundary

} // namespace geos

#endif // GEOS_MIXEDMIMETIC_MIXEDMIMETICBOUNDARYCONDITIONS_HPP_
