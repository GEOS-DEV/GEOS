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
 * @file SlipDependentPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_SLIPDEPENDENTPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_SLIPDEPENDENTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class SlipDependentPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  SlipDependentPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                           arrayView3d< real64 > const & dPerm_dPressure,
                                           arrayView4d< real64 > const & dPerm_dDispJump,
                                           real64 const shearDispThreshold,
                                           real64 const maxPermMultiplier,
                                           arrayView3d< real64 > const & initialPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dDispJump( dPerm_dDispJump ),
    m_shearDispThreshold( shearDispThreshold ),
    m_maxPermMultiplier( maxPermMultiplier ),
    m_initialPermeability( initialPermeability )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const ( &dispJump )[3],
                arraySlice1d<real64> const & initialPermeability,
                arraySlice1d< real64 > const & permeability,
                arraySlice2d< real64 > const & dPerm_dDispJump ) const;

  GEOSX_HOST_DEVICE
  virtual void updateFromApertureAndShearDisplacement( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & oldHydraulicAperture,
                                                       real64 const & newHydraulicAperture,
                                                       real64 const ( &dispJump )[3] ) const override
  {
    GEOSX_UNUSED_VAR( q, oldHydraulicAperture, newHydraulicAperture );

    compute( dispJump,
             m_initialPermeability[k][0],
             m_permeability[k][0],
             m_dPerm_dDispJump[k][0] );
  }

private:

  /// Derivative of fracture permeability to shear displacement jump between fracture surfaces
  arrayView4d< real64 > m_dPerm_dDispJump;

  /// Threshold of shear displacement
  real64 m_shearDispThreshold;

  /// Maximum permeability multiplier
  real64 m_maxPermMultiplier;

  /// Initial permeability tensor
  arrayView3d< real64 > m_initialPermeability;

};


class SlipDependentPermeability : public PermeabilityBase
{
public:

  SlipDependentPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "SlipDependentPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void initializeState() const override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = SlipDependentPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dDispJump,
                          m_shearDispThreshold,
                          m_maxPermMultiplier,
                          m_initialPermeability );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dDispJumpString() { return "dPerm_dDispJump"; }
    static constexpr char const * shearDispThresholdString() { return "shearDispThreshold"; }
    static constexpr char const * maxPermMultiplierString() { return "maxPermMultiplier"; }
    static constexpr char const * initialPermeabilityString() { return "iniPermeability"; }
  } ;

private:

  /// Derivative of fracture permeability to shear displacement jump between fracture surfaces
  array4d< real64 > m_dPerm_dDispJump;

  /// Threshold of shear displacement
  real64 m_shearDispThreshold;

  /// Maximum permeability multiplier
  real64 m_maxPermMultiplier;

  /// Initial permeability tensor
  array3d< real64 > m_initialPermeability;

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SlipDependentPermeabilityUpdate::compute( real64 const ( &dispJump )[3],
                                               arraySlice1d< real64 > const & initialPermeability,
                                               arraySlice1d< real64 > const & permeability,
                                               arraySlice2d< real64 > const & dPerm_dDispJump ) const
{ 
  real64 const shearMag = std::sqrt( dispJump[1]*dispJump[1] + dispJump[2]*dispJump[2] );
  
  real64 const tmpTanh = std::tanh ( 3.0 * shearMag/m_shearDispThreshold );
  
  real64 const permMultiplier = ( m_maxPermMultiplier - 1.0 ) * tmpTanh + 1.0;
  
  real64 const dpermMultiplier_dshearMag = ( m_maxPermMultiplier - 1.0 ) * ( 1.0 - tmpTanh * tmpTanh ) * 3.0/m_shearDispThreshold;  

  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = permMultiplier * initialPermeability[i];
    dPerm_dDispJump[i][0] = 0.0;
    real64 const tmpValue = shearMag > 0.0 ? initialPermeability[i] * dpermMultiplier_dshearMag /shearMag : 0.0;
    dPerm_dDispJump[i][1] = tmpValue * dispJump[1];
    dPerm_dDispJump[i][2] = tmpValue * dispJump[2];
  }
}



} /* namespace constitutive */

} /* namespace geosx */

#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
