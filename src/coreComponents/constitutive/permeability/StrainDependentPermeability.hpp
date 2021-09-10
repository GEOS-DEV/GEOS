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
 * @file StrainDependentPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class StrainDependentPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  StrainDependentPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                     arrayView3d< real64 > const & dPerm_dPressure,
                                     arrayView3d< real64 > const & dPerm_dStrain,
                                     real64 const strainThreshold,
                                     real64 const maxPermMultiplier,
                                     arrayView3d< real64 > const & iniPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dStrain( dPerm_dStrain ),
    m_strainThreshold( strainThreshold ),
    m_maxPermMultiplier( maxPermMultiplier ),
    m_iniPermeability( iniPermeability )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const ( &totalStrain )[6],
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dStrain ) const;

  GEOSX_HOST_DEVICE
  virtual void updateFromStrain( localIndex const k,
                                 localIndex const q,
                                 real64 const ( &totalStrain )[6] ) const override
  {
    compute( strain,
             m_permeability[k][q],
             m_dPerm_dStrain[k][q] );
  }

private:

  /// dPermeability_dStrain
  arrayView3d< real64 > m_dPerm_dStrain;

  /// Threshold of shear strain
  real64 m_strainThreshold;

  /// Maximum permeability multiplier
  real64 m_maxPermMultiplier;

  /// Initial permeability tensor
  arrayView3d< real64 > m_iniPermeability;

};


class StrainDependentPermeability : public PermeabilityBase
{
public:

  StrainDependentPermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "StrainDependentPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = StrainDependentPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dStrain,
                          m_strainThreshold,
                          m_maxPermMultiplier,
                          m_iniPermeability );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dStrainString() { return "dPerm_dStrain"; }
    static constexpr char const * strainThresholdString() { return "strainThreshold"; }
    static constexpr char const * maxPermMultiplierString() { return "maxPermMultiplier"; }
    static constexpr char const * iniPermeabilityString() { return "iniPermeability"; }
  } viewKeys;

private:

  /// dPermeability_dStrain
  arrayView3d< real64 > m_dPerm_dStrain;

  /// Threshold of shear strain
  real64 m_strainThreshold;

  /// Maximum permeability multiplier
  real64 m_maxPermMultiplier;

  /// Initial permeability tensor
  arrayView3d< real64 > m_iniPermeability;

};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void StrainDependentPermeabilityUpdate::compute( real64 const ( &totalStrain )[6],
                                                 arraySlice1d< real64 > const & permeability,
                                                 arraySlice1d< real64 > const & dPerm_dStrain  ) const
{ 
  real64 const shearStrain = std::max( abs(totalStrain[3]), abs(totalStrain[4]), abs(totalStrain[5]) )
   
  if ( shearStrain < m_strainThreshold )
    {
        real64 const permMultiplier = m_maxPermMultiplier * shearStrain/m_strainThreshold;

        real64 const dPerm_dStrainValue = m_maxPermMultiplier/m_strainThreshold;
    }
    else
    {
        real64 const permMultiplier = m_maxPermMultiplier;

        real64 const dPerm_dStrainValue = 0;
    }

  
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = permMultiplier * m_iniPermeability[i];
    dPerm_dStrain[i] = dPerm_dStrainValue * m_iniPermeability[i];
  }
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
