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
 * @file LinearIsotropicDispersion.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DISPERSION_LINEARISOTROPICDISPERSION_HPP_
#define GEOS_CONSTITUTIVE_DISPERSION_LINEARISOTROPICDISPERSION_HPP_

#include "constitutive/dispersion/Layout.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/DataLayouts.hpp"
#include "constitutive/dispersion/DispersionBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief The update class for linear isotropic diffusion
 */
class LinearIsotropicDispersionUpdate : public DispersionBaseUpdate
{
public:

  /**
   * @brief Constructor for the class performing the dispersion updates
   * @param dispersivity the array of cell-wise dispersivities in the subregion
   * @param longitunidalDispersivity longitudinal dispersivity in the subregion
   */
  LinearIsotropicDispersionUpdate(arrayView3d<real64> const &dispersivity,
                                  real64 const &longitudinalDispersivity,
                                  arrayView3d<real64> const &phaseVelocityNorm)
    : DispersionBaseUpdate( dispersivity ),
    m_longitudinalDispersivity( longitudinalDispersivity ),
    m_phaseVelocityNorm(phaseVelocityNorm)
  {}

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       arraySlice2d< real64 const > const & laggedTotalVelocityComponents ) const override
  {
    for( integer i = 0; i < 3; ++i )
    {
        m_dispersivity[k][q][i] = m_longitudinalDispersivity;
      for( int ip = 0; ip < laggedTotalVelocityComponents.size(0); ++ip )
      {
        real64 const velocityNorm = LvArray::tensorOps::l2Norm< 3 >( laggedTotalVelocityComponents[ip] );
        m_phaseVelocityNorm[k][q][ip] = velocityNorm;
      }
    }
  }

protected:

  /// Longitudinal dispersivity
  real64 const m_longitudinalDispersivity;
  arrayView3d< real64 > const m_phaseVelocityNorm;

};

/**
 * @brief The class for constant dispersion
 */
class LinearIsotropicDispersion : public DispersionBase
{
public:

  /**
   * @brief Constructor for the class storing constant dispersion
   * @param[in] name the name of the class
   * @param[in] parent pointer to the parent Group
   */
  LinearIsotropicDispersion( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  static string catalogName() { return "LinearIsotropicDispersion"; }

  virtual void allocateConstitutiveData( dataRepository::Group & parent, const geos::localIndex numConstitutivePointsPerParentIndex ) override;

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void initializeVelocityState( arrayView4d< real64 const > const & initialVelocity ) const override;

  virtual void saveConvergedVelocityState( arrayView4d< real64 const > const & convergedVelocity ) const override;

    /**
      * @brief Getter for phase Velolcity norm in the subRegion
      * @return an arrayView of the phase velocities
      * @note it is a *lagged* version of the velocities and not only the norm of it
      */
    arrayView3d< real64 const > phaseVelocityNorm() const { return m_phaseVelocityNorm; }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = LinearIsotropicDispersionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_dispersivity,
                          m_longitudinalDispersivity,
                          m_phaseVelocityNorm );
  }

  struct viewKeyStruct : public DispersionBase::viewKeyStruct
  {
    static constexpr char const * longitudinalDispersivityString() { return "longitudinalDispersivity"; }
  } viewKeys;

protected:

  virtual void postProcessInput() override;

private:

  /// Longitudinal dispersivity
  real64 m_longitudinalDispersivity;
  array3d< real64, dispersion::LAYOUT_PHASE_VELOCITY_NORM > m_phaseVelocityNorm;

};

} // namespace constitutive

} // namespace geos


#endif // GEOS_CONSTITUTIVE_DISPERSION_LINEARISOTROPICDISPERSION_HPP_
