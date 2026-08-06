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
 * @file SurfaceAreaBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SURFACEAREA_SURFACEAREABASE_HPP_
#define GEOS_CONSTITUTIVE_SURFACEAREA_SURFACEAREABASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable update routines for the reactive surface area.
 *
 * The base implementation keeps the surface area at its initial value. Derived models
 * override updateFromPorosityAndVolumeFractions to describe how the reactive surface
 * area of each mineral evolves as the rock dissolves or precipitates.
 */
class SurfaceAreaBaseUpdates
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_surfaceArea.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_surfaceArea.size( 1 ); }

  /**
   * @brief Get the number of kinetic reactions.
   * @return number of kinetic reactions
   */
  GEOS_HOST_DEVICE
  integer numKineticReactions() const { return LvArray::integerConversion< integer >( m_surfaceArea.size( 2 ) ); }

  /**
   * @brief Update the reactive surface area of all the minerals in one quadrature point.
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the current porosity
   * @param[in] initialPorosity the initial porosity
   * @param[in] volumeFractions the current mineral volume fractions
   * @param[in] initialVolumeFractions the initial mineral volume fractions
   */
  GEOS_HOST_DEVICE
  virtual void updateFromPorosityAndVolumeFractions( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & porosity,
                                                     real64 const & initialPorosity,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & volumeFractions,
                                                     arraySlice1d< real64 const, reactivefluid::USD_SPECIES - 2 > const & initialVolumeFractions ) const
  {
    GEOS_UNUSED_VAR( k, q, porosity, initialPorosity, volumeFractions, initialVolumeFractions );
  }

protected:

  SurfaceAreaBaseUpdates( arrayView3d< real64, reactivefluid::USD_SPECIES > const & surfaceArea,
                          arrayView3d< real64 const, reactivefluid::USD_SPECIES > const & initialSurfaceArea ):
    m_surfaceArea( surfaceArea ),
    m_initialSurfaceArea( initialSurfaceArea )
  {}

  /// Reactive surface area of each mineral
  arrayView3d< real64, reactivefluid::USD_SPECIES > m_surfaceArea;

  /// Initial reactive surface area of each mineral
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > m_initialSurfaceArea;
};


/**
 * @brief Base class for all the models describing the reactive surface area of the minerals.
 *
 * The model owns the reactive surface area used by the kinetic reactions, one value per
 * mineral and per quadrature point. The number of minerals is deduced from the size of the
 * defaultInitialSurfaceArea input.
 */
class SurfaceAreaBase : public ConstitutiveBase
{
public:

  SurfaceAreaBase( string const & name, dataRepository::Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                dataRepository::Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Number of minerals (kinetic reactions) handled by this model.
   * @return the number of kinetic reactions
   */
  integer numKineticReactions() const { return m_numKineticReactions; }

  /**
   * @brief Const accessor for the current reactive surface area.
   * @return the reactive surface area of each mineral
   */
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > getSurfaceArea() const { return m_surfaceArea; }

  /**
   * @brief Const accessor for the initial reactive surface area.
   * @return the initial reactive surface area of each mineral
   */
  arrayView3d< real64 const, reactivefluid::USD_SPECIES > getInitialSurfaceArea() const { return m_initialSurfaceArea; }

  /**
   * @brief Set the surface area to its initial value.
   */
  virtual void initializeState() const;

  struct viewKeyStruct
  {
    static constexpr char const * defaultInitialSurfaceAreaString() { return "defaultInitialSurfaceArea"; }
  } viewKeys;

protected:

  virtual void postInputInitialization() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts );

  /// Default initial surface area, one value per mineral
  array1d< real64 > m_defaultInitialSurfaceArea;

  /// Initial reactive surface area of each mineral
  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_initialSurfaceArea;

  /// Current reactive surface area of each mineral
  array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES > m_surfaceArea;

  /// Number of minerals (kinetic reactions)
  integer m_numKineticReactions;
};

}/* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_SURFACEAREA_SURFACEAREABASE_HPP_
