/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParticleFluidBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

class ParticleFluidBase : public ConstitutiveBase
{
public:

  ParticleFluidBase( std::string const & name, Group * const parent );

  virtual ~ParticleFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const override = 0;

  virtual void AllocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** ParticleFluidBase-specific interface

  virtual void PointUpdate( localIndex const NC, real64 const & proppantConcentration, arraySlice1d< real64 const > const & ComponentConcentration,
                            arraySlice1d< real64 const > const & nIndex, arraySlice1d< real64 const > const & KIndex, real64 const & fluidDensity,
                            real64 const & dFluidDensity_dPressure, arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                            localIndex const k ) = 0;

  virtual void PointUpdate( localIndex const NC, real64 const & proppantConcentration, real64 const & fluidDensity, real64 const & dFluidDensity_dPressure,
                            arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration, real64 const & fluidViscosity,
                            real64 const & dFluidViscosity_dPressure, arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration,
                            localIndex const k ) = 0;

  virtual void BatchUpdate( arrayView1d< real64 const > const & concentration ) = 0;

  static constexpr localIndex MAX_NUM_COMPONENTS = 4;

  // *** Data repository keys

  struct viewKeyStruct
  {

    static constexpr auto settlingFactorString    = "settlingFactor";
    static constexpr auto dSettlingFactor_dPressureString  = "dSettlingFactor_dPressure";
    static constexpr auto dSettlingFactor_dProppantConcentrationString  = "dSettlingFactor_dProppantConcentration";
    static constexpr auto dSettlingFactor_dComponentConcentrationString  = "dSettlingFactor_dComponentConcentration";

    static constexpr auto collisionFactorString    = "collisionFactor";
    static constexpr auto dCollisionFactor_dProppantConcentrationString  = "dCollisionFactor_dProppantConcentration";

    static constexpr auto maxProppantConcentrationString    = "maxProppantConcentration";

    static constexpr auto isCollisionalSlipString    = "isCollisionalSlip";

    static constexpr auto proppantPackPermeabilityString    = "proppantPackPermeability";
  } viewKeysParticleFluidBase;

protected:

  virtual void PostProcessInput() override;

  array1d< real64 > m_settlingFactor;
  array1d< real64 > m_dSettlingFactor_dPressure;
  array1d< real64 > m_dSettlingFactor_dProppantConcentration;
  array2d< real64 > m_dSettlingFactor_dComponentConcentration;

  array1d< real64 > m_collisionFactor;
  array1d< real64 > m_dCollisionFactor_dProppantConcentration;

  array1d< real64 > m_proppantPackPermeability;

  integer m_isCollisionalSlip;

  real64 m_maxProppantConcentration;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_PARTICLEFLUIDBASE_HPP
