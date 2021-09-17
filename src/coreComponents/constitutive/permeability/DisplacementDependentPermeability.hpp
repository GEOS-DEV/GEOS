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
 * @file DisplacementDependentPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_DISPLACEMENTDEPENDENTPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_DISPLACEMENTDEPENDENTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class DisplacementDependentPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  DisplacementDependentPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                           arrayView3d< real64 > const & dPerm_dPressure,
                                           arrayView3d< real64 > const & dPerm_dDisplacement,
                                           real64 const shearDispThreshold,
                                           real64 const maxPermMultiplier,
                                           arrayView3d< real64 > const & iniPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dStrain( dPerm_dDisplacement ),
    m_shearDispThreshold( shearDispThreshold ),
    m_maxPermMultiplier( maxPermMultiplier ),
    m_iniPermeability( iniPermeability )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const ( &displacementJump )[3],
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dDisplacement ) const;

  GEOSX_HOST_DEVICE
  virtual void updateFromShearDisplacement( localIndex const k,
                                            localIndex const q,
                                            real64 const ( &displacementJump )[3] ) const override
  {
    compute( displacementJump,
             m_permeability[k][q],
             m_dPerm_dDisplacement[k][q] );
  }

private:

  /// dPermeability_dDisplacement
  arrayView3d< real64 > m_dPerm_dDisplacement;

  /// Threshold of shear displacement
  real64 m_shearDispThreshold;

  /// Maximum permeability multiplier
  real64 m_maxPermMultiplier;

  /// Initial permeability tensor
  arrayView3d< real64 > m_iniPermeability;

};


class DisplacementDependentPermeability : public PermeabilityBase
{
public:

  DisplacementDependentPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "DisplacementDependentPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = DisplacementDependentPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dDisplacement,
                          m_shearDispThreshold,
                          m_maxPermMultiplier,
                          m_iniPermeability );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dDisplacementString() { return "dPerm_dDisplacement"; }
    static constexpr char const * shearDispThresholdString() { return "shearDispThreshold"; }
    static constexpr char const * maxPermMultiplierString() { return "maxPermMultiplier"; }
    static constexpr char const * iniPermeabilityString() { return "iniPermeability"; }
  } viewKeys;

private:

  /// dPermeability_dStrain
  arrayView3d< real64 > m_dPerm_dDisplacement;

  /// Threshold of shear strain
  real64 m_shearDispThreshold;

  /// Maximum permeability multiplier
  real64 m_maxPermMultiplier;

  /// Initial permeability tensor
  arrayView3d< real64 > m_iniPermeability;

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void DisplacementDependentPermeabilityUpdate::compute( real64 const ( &displacementJump )[3],
                                                       arraySlice1d< real64 > const & permeability,
                                                       arraySlice1d< real64 > const & dPerm_dDisplacement  ) const
{ 
  real64 const shearDisp = std::max( abs(displacementJump[1]), abs(displacementJump[2]) )
   
  if (  shearDisp < m_shearDispThreshold )
    {
        real64 const permMultiplier = (m_maxPermMultiplier - 1) * shearDisp/m_shearDispThreshold + 1;

        real64 const dPerm_dDisplacementValue = (m_maxPermMultiplier - 1)/m_shearDispThreshold;
    }
    else
    {
        real64 const permMultiplier = m_maxPermMultiplier;

        real64 const dPerm_dDisplacement = 0;
    }

  
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = permMultiplier * m_iniPermeability[i];
    dPerm_dDisplacement[i] = dPerm_dDisplacementValue * m_iniPermeability[i];
  }
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
