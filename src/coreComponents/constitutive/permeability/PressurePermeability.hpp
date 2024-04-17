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
 * @file PressurePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_PRESSUREPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_PRESSUREPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

class PressurePermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  PressurePermeabilityUpdate( R1Tensor const pressureDependenceConstants,
                              real64 const & referencePressure,
                              arrayView3d< real64 > const & referencePermeability,
                              arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_pressureDependenceConstants( pressureDependenceConstants ),
    m_referencePressure( referencePressure ),
    m_referencePermeability( referencePermeability )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & deltaPressure,
                R1Tensor const pressureDependenceConstants,
                real64 const (&referencePermeability)[3],
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndPorosity( localIndex const k,
                                              localIndex const q,
                                              real64 const & pressure,
                                              real64 const & porosity ) const override
  {
    GEOS_UNUSED_VAR( q, porosity );

    real64 const deltaPressure = pressure - m_referencePressure;

    real64 referencePermeability[3];

    referencePermeability[0] = m_referencePermeability[k][0][0];
    referencePermeability[1] = m_referencePermeability[k][0][1];
    referencePermeability[2] = m_referencePermeability[k][0][2];

    compute( deltaPressure,
             m_pressureDependenceConstants,
             referencePermeability,
             m_permeability[k][0],
             m_dPerm_dPressure[k][0] );
  }

private:

  /// Pressure dependent coefficients for each permeability component
  R1Tensor m_pressureDependenceConstants;

  /// Reference pressure in the model
  real64 const m_referencePressure;

  arrayView3d< real64 > m_referencePermeability;

};


class PressurePermeability : public PermeabilityBase
{
public:

  PressurePermeability( string const & name, Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PressurePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PressurePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_pressureDependenceConstants,
                          m_referencePressure,
                          m_referencePermeability,
                          m_permeability,
                          m_dPerm_dPressure );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * referencePermeabilityComponentsString() { return "referencePermeabilityComponents"; }
    static constexpr char const * pressureDependenceConstantsString() { return "pressureDependenceConstants"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referencePermeabilityString() { return "referencePermeability"; }
  } viewKeys;

  virtual void initializeState() const override final;

private:

  /// Permeability components at the reference pressure
  R1Tensor m_referencePermeabilityComponents;

  /// Pressure dependent coefficients for each permeability component
  R1Tensor m_pressureDependenceConstants;

  /// Reference pressure in the model
  real64 m_referencePressure;

  array3d< real64 > m_referencePermeability;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PressurePermeabilityUpdate::compute( real64 const & deltaPressure,
                                          R1Tensor const pressureDependenceConstants,
                                          real64 const (&referencePermeability)[3],
                                          arraySlice1d< real64 > const & permeability,
                                          arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const perm = referencePermeability[i] * exp( pressureDependenceConstants[i] * deltaPressure ); // To use ExponentialRelation for
                                                                                                          // this

    permeability[i] = perm;
    dPerm_dPressure[i] = perm * pressureDependenceConstants[i];
  }
}

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_PRESSUREPERMEABILITY_HPP_
