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
 * @file ReactivePressurePermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_REACTIVEPRESSUREPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_REACTIVEPRESSUREPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PressurePermeability.hpp"


namespace geos
{
namespace constitutive
{

class ReactivePressurePermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ReactivePressurePermeabilityUpdate( PressureModelType const & presModelType,
                                      R1Tensor const pressureDependenceConstants,
                                      arrayView1d< real64 const > const & referencePressure,
                                      real64 const & maxPermeability,
                                      real64 const & minPermeability,
                                      arrayView3d< real64 const > const & referencePermeability,
                                      arrayView1d< real64 const > const & referencePorosity,
                                      arrayView3d< real64 > const & permeability,
                                      arrayView3d< real64 const > const & permeability_n,
                                      arrayView3d< real64 > const & dPerm_dPressure,
                                      integer const & explicitFlag )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_presModelType( presModelType ),
    m_pressureDependenceConstants( pressureDependenceConstants ),
    m_referencePressure( referencePressure ),
    m_maxPermeability( maxPermeability ),
    m_minPermeability( minPermeability ),
    m_referencePermeability( referencePermeability ),
    m_referencePorosity( referencePorosity ),
    m_permeability_n( permeability_n ),
    m_explicitFlag( explicitFlag )
  {}

  GEOS_HOST_DEVICE
  void computePressure( real64 const & deltaPressure,
                        R1Tensor const pressureDependenceConstants,
                        real64 const (&referencePermeability)[3],
                        arraySlice1d< real64 > const & permeability,
                        arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  void computePressure( real64 const & deltaPressure,
                        R1Tensor const pressureDependenceConstants,
                        real64 const (&referencePermeability)[3],
                        real64 const maxPermeability,
                        arraySlice1d< real64 > const & permeability,
                        arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  void computeWithPorosity( real64 const & porosity,
                            real64 const & minPermeability,
                            real64 const & referencePorosity,
                            real64 const (&referencePermeability)[3],
                            arraySlice1d< real64 > const & permeability,
                            arraySlice1d< real64 > const & dPerm_dPressure ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndPorosity( localIndex const k,
                                              localIndex const q,
                                              real64 const & pressure,
                                              real64 const & pressure_n,
                                              real64 const & porosity ) const override
  {
    GEOS_UNUSED_VAR( q );

    real64 deltaPressure = (m_explicitFlag)? (pressure_n - m_referencePressure[k]):(pressure - m_referencePressure[k]);

    real64 referencePermeability[3];

    referencePermeability[0] = m_referencePermeability[k][0][0];
    referencePermeability[1] = m_referencePermeability[k][0][1];
    referencePermeability[2] = m_referencePermeability[k][0][2];

    real64 permeability_n[3];

    permeability_n[0] = m_permeability_n[k][0][0];
    permeability_n[1] = m_permeability_n[k][0][1];
    permeability_n[2] = m_permeability_n[k][0][2];

    switch( m_presModelType )
    {
      case PressureModelType::Exponential:
      {
        // Step 1: Update with pressure
        computePressure( deltaPressure,
                         m_pressureDependenceConstants,
                         referencePermeability,
                         m_permeability[k][0],
                         m_dPerm_dPressure[k][0] );

        // Step 2: Update with reactive porosity
        computeWithPorosity( porosity, m_minPermeability, m_referencePorosity[k], referencePermeability, m_permeability[k][0], m_dPerm_dPressure[k][0] );

        if( m_explicitFlag )
        {
          for( localIndex i=0; i < m_permeability[k][0].size(); i++ )
          {
            m_dPerm_dPressure[k][0][i] = 0.0;

            if( m_maxPermeability < 1.0 )
            {
              real64 perm = m_permeability[k][0][i];

              if( perm < referencePermeability[i] )
              {
                perm = referencePermeability[i];
              }

              // Ensure permeability non-decreasing
              if( perm < permeability_n[i] )
              {
                perm = permeability_n[i];
              }

              m_permeability[k][0][i] = (perm > m_maxPermeability)? m_maxPermeability:perm;
            }
          }
        }

        break;
      }
      case PressureModelType::Hyperbolic:
      {
        // Step 1: Update with pressure
        computePressure( deltaPressure,
                         m_pressureDependenceConstants,
                         referencePermeability,
                         m_maxPermeability,
                         m_permeability[k][0],
                         m_dPerm_dPressure[k][0] );

        // Step 2: Update with reactive porosity
        computeWithPorosity( porosity, m_minPermeability, m_referencePorosity[k], referencePermeability, m_permeability[k][0], m_dPerm_dPressure[k][0] );

        if( m_explicitFlag )
        {
          for( localIndex i=0; i < m_permeability[k][0].size(); i++ )
          {
            m_dPerm_dPressure[k][0][i] = 0.0;

            real64 const perm = m_permeability[k][0][i];

            if( perm < referencePermeability[i] )
            {
              m_permeability[k][0][i] = referencePermeability[i];
            }
          }
        }
        break;
      }
      default:
      {
        GEOS_ERROR( "PressureModelType is invalid! It should be either Exponential or Hyperbolic" );
      }
    }
  }

private:

  /// Pressure dependence model type
  PressureModelType m_presModelType;

  /// Pressure dependent coefficients for each permeability component
  R1Tensor m_pressureDependenceConstants;

  /// Reference pressure in the model
  arrayView1d< real64 const >  const m_referencePressure;

  /// Maximum and minimum permeability
  real64 const m_maxPermeability;
  real64 const m_minPermeability;

  arrayView3d< real64 const > const m_referencePermeability;
  arrayView1d< real64 const > const m_referencePorosity;

  arrayView3d< real64 const > const m_permeability_n;

  integer m_explicitFlag;

};


class ReactivePressurePermeability : public PermeabilityBase
{
public:

  ReactivePressurePermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "ReactivePressurePermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ReactivePressurePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_presModelType,
                          m_pressureDependenceConstants,
                          m_referencePressure,
                          m_maxPermeability,
                          m_minPermeability,
                          m_referencePermeability,
                          m_referencePorosity,
                          m_permeability,
                          m_permeability_n,
                          m_dPerm_dPressure,
                          m_explicitFlag );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * referencePermeabilityComponentsString() { return "referencePermeabilityComponents"; }
    static constexpr char const * pressureDependenceConstantsString() { return "pressureDependenceConstants"; }
    static constexpr char const * defaultReferencePressureString() { return "defaultReferencePressure"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * referencePermeabilityString() { return "referencePermeability"; }
    static constexpr char const * defaultReferencePorosityString() { return "defaultReferencePorosity"; }
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * maxPermeabilityString() { return "maxPermeability"; }
    static constexpr char const * minPermeabilityString() { return "minPermeability"; }
    static constexpr char const * pressureModelTypeString() { return "pressureModelType"; }
    static constexpr char const * explicitUpdateFlagString() { return "explicitUpdateFlag"; }
  };

  virtual void initializeState() const override final;

protected:

  virtual void postInputInitialization() override;

private:

  /// Permeability components at the reference pressure
  R1Tensor m_referencePermeabilityComponents;

  real64 m_defaultReferencePorosity;

  /// Pressure dependent coefficients for each permeability component
  R1Tensor m_pressureDependenceConstants;

  /// Reference pressure in the model
  real64 m_defaultReferencePressure;
  array1d< real64 > m_referencePressure;

  /// Maximum and minimum permeability
  real64 m_maxPermeability;
  real64 m_minPermeability;

  array3d< real64 > m_referencePermeability;
  array1d< real64 > m_referencePorosity;

  /// Pressure dependence model type
  PressureModelType m_presModelType;

  /// Explicit update flag
  integer m_explicitFlag;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ReactivePressurePermeabilityUpdate::computePressure( real64 const & deltaPressure,
                                                          R1Tensor const pressureDependenceConstants,
                                                          real64 const (&referencePermeability)[3],
                                                          arraySlice1d< real64 > const & permeability,
                                                          arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const perm = referencePermeability[i] * exp( pressureDependenceConstants[i] * deltaPressure );

    permeability[i] = perm;
    dPerm_dPressure[i] = perm * pressureDependenceConstants[i];
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ReactivePressurePermeabilityUpdate::computePressure( real64 const & deltaPressure,
                                                          R1Tensor const pressureDependenceConstants,
                                                          real64 const (&referencePermeability)[3],
                                                          real64 const maxPermeability,
                                                          arraySlice1d< real64 > const & permeability,
                                                          arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const pressureOffSet = log( maxPermeability/referencePermeability[i] - 1 )/pressureDependenceConstants[i];

    real64 const perm = maxPermeability/( 1 + exp( -pressureDependenceConstants[i]*( deltaPressure - pressureOffSet ) ) );
    permeability[i] = perm;
    dPerm_dPressure[i] = perm*perm/maxPermeability*pressureDependenceConstants[i]*exp( -pressureDependenceConstants[i]*deltaPressure );
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ReactivePressurePermeabilityUpdate::computeWithPorosity( real64 const & porosity,
                                                              real64 const & minPermeability,
                                                              real64 const & referencePorosity,
                                                              real64 const (&referencePermeability)[3],
                                                              arraySlice1d< real64 > const & permeability,
                                                              arraySlice1d< real64 > const & dPerm_dPressure ) const
{
  for( localIndex i=0; i < permeability.size(); i++ )
  {
    real64 const permScale = (minPermeability + ( referencePermeability[i] - minPermeability ) * pow( porosity/referencePorosity, 6 )) / referencePermeability[i];

    real64 const perm = permeability[i];
    permeability[i] = perm * permScale;

    real64 const dPerm_dPres = dPerm_dPressure[i];
    dPerm_dPressure[i] = dPerm_dPres * permScale;
  }
}

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_REACTIVEPRESSUREPERMEABILITY_HPP_
