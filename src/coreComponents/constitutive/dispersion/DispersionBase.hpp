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
 * @file DispersionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONBASE_HPP
#define GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONBASE_HPP

#include "common/DataLayouts.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/dispersion/Layout.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief The abstract base class to perform the dispersion update
 */
class DispersionBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_dispersivity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_dispersivity.size( 1 ); }

protected:

  /**
   * @brief Constructor for the class performing the dispersion updates
   * @param dispersivity the array of cell-wise dispersion in the subregion
   */
  DispersionBaseUpdate( arrayView4d< real64 > const & dispersivity )
    : m_dispersivity( dispersivity )
  {}

  /// View on the cell-wise dispersivity
  arrayView4d< real64 > const m_dispersivity;

private:

  /**
   * @brief Pointwise update function called from the solver
   * @param[in] k index of the cell in the subRegion
   * @param[in] q constitutive index (equal to one in this class)
   * @param[in] laggedTotalVelocityComponents the components of the total velocity
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice2d< real64 const > const & laggedTotalVelocityComponents,
                       arraySlice1d< real64 const > const & phaseDensity ) const = 0;
};

/**
 * @brief The abstract base class for dispersion
 */
class DispersionBase : public ConstitutiveBase
{
public:

  /// Max number of phases allowed in the class
  static constexpr integer MAX_NUM_PHASES = 3;

  /**
   * @brief Constructor for the abstract base class
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  DispersionBase( string const & name, dataRepository::Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Getter for the number of fluid phases
   * @return the number of fluid phases
   */
  integer numFluidPhases() const { return LvArray::integerConversion< integer >( m_phaseNames.size() ); }

  /**
   * @brief Getter for the dispersivities in the subRegion
   * @return an arrayView of dispersivities
   */
  arrayView4d< real64 const > dispersivity() const { return m_dispersivity; }

  /**
   * @brief Getter for phase Velolcity in the subRegion
   * @return an arrayView of the phase velocities
   */
  //TODO reduce gauss point dim
  arrayView4d< real64 const > phaseVelocity() const { return m_phaseVelocity; }

  /**
   * @brief Initialize the velocity state (needed because dispersion depends on total velocity)
   * @param[in] initialVelocity the initial velocity field after reservoir initialization
   *
   * Note: this is needed because for now, the velocity field is treated **explicitly** in the dispersion tensor
   */
  virtual void initializeVelocityState( arrayView4d< real64 const > const & initialVelocity, arrayView3d< real64 const > const & phaseDensity ) const
  { GEOS_UNUSED_VAR( initialVelocity ); }

  /**
   * @brief Save the velocity state (needed because dispersion depends on total velocity)
   * @param[in] convergedVelocity the converged velocity field
   *
   * Note: this is needed because for now, the velocity is treated **explicitly** in the dispersion tensor
   */
  virtual void saveConvergedVelocityState( arrayView4d< real64 const > const & convergedVelocity, arrayView3d< real64 const > const & phaseDensity ) const
  { GEOS_UNUSED_VAR( convergedVelocity ); }


  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * phaseNamesString() { return "phaseNames"; }
  };

protected:

  virtual void postProcessInput() override;

  /// cell-wise dispersivity in the subregion
  /// TODO: support full tensor if linear isotropic diffusion is no longer enough
  array4d< real64 > m_dispersivity;

  /// phase names read from input
  string_array m_phaseNames;
  // misc
  array4d< real64, dispersion::LAYOUT_PHASE_VELOCITY > m_phaseVelocity;

};

} // namespace constitutive

} // namespace geos


#endif // GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONBASE_HPP
