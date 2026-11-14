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
 * @file StrainDependentPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @brief Kernel wrapper for the volumetric strain dependent permeability model (Chin et al., 2000):
 *        phi(e_v) = 1 - (1 - phi_0) exp( -e_v ),   k_i / k_0,i = ( phi(e_v) / phi_0 )^n_i
 *        where e_v is the total volumetric strain since the initial state (positive in extension).
 */
class StrainDependentPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  StrainDependentPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                     arrayView3d< real64 > const & dPerm_dPressure,
                                     arrayView3d< real64 > const & dPerm_dVolStrain,
                                     arrayView3d< real64 const > const & referencePermeability,
                                     arrayView2d< real64 > const & volStrain,
                                     arrayView2d< real64 const > const & volStrain_n,
                                     R1Tensor const permeabilityExponents )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dVolStrain( dPerm_dVolStrain ),
    m_referencePermeability( referencePermeability ),
    m_volStrain( volStrain ),
    m_volStrain_n( volStrain_n ),
    m_permeabilityExponents( permeabilityExponents )
  {}

  /**
   * @brief Compute the permeability and its derivative w.r.t. the volumetric strain
   * @param[in] initialPorosity the initial porosity phi_0
   * @param[in] referencePermeability the permeability k_0 at the initial state
   * @param[in] volStrain the total volumetric strain e_v since the initial state
   * @param[out] permeability the permeability
   * @param[out] dPerm_dVolStrain the derivative of the permeability w.r.t. e_v
   */
  GEOS_HOST_DEVICE
  void compute( real64 const & initialPorosity,
                arraySlice1d< real64 const > const & referencePermeability,
                real64 const & volStrain,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dVolStrain ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromVolumetricStrainIncrement( localIndex const k,
                                                    localIndex const q,
                                                    real64 const & initialPorosity,
                                                    real64 const & volStrainIncrement ) const override
  {
    m_volStrain[k][q] = m_volStrain_n[k][q] + volStrainIncrement;

    // The permeability is cell-centered (one value per element): it is computed from the
    // element-averaged volumetric strain once all the quadrature points of element k have been updated.
    // This assumes the quadrature points of an element are visited in order, as in the FE kernels.
    localIndex const numQuad = m_volStrain.size( 1 );
    if( q == numQuad - 1 )
    {
      real64 avgVolStrain = 0.0;
      for( localIndex iq = 0; iq < numQuad; ++iq )
      {
        avgVolStrain += m_volStrain[k][iq];
      }
      avgVolStrain /= numQuad;

      compute( initialPorosity,
               m_referencePermeability[k][0],
               avgVolStrain,
               m_permeability[k][0],
               m_dPerm_dVolStrain[k][0] );
    }
  }

private:

  /// Minimum porosity used to keep the porosity ratio defined
  static constexpr real64 minPorosity = 1.0e-6;

  /// Derivative of the permeability w.r.t. the volumetric strain
  arrayView3d< real64 > m_dPerm_dVolStrain;

  /// Permeability at the initial state
  arrayView3d< real64 const > m_referencePermeability;

  /// Total volumetric strain since the initial state
  arrayView2d< real64 > m_volStrain;

  /// Total volumetric strain at the last converged state
  arrayView2d< real64 const > m_volStrain_n;

  /// Exponents of the porosity ratio for each permeability component
  R1Tensor m_permeabilityExponents;
};


class StrainDependentPermeability : public PermeabilityBase
{
public:

  StrainDependentPermeability( string const & name, dataRepository::Group * const parent );

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    dataRepository::Group * const parent ) const override;

  static string catalogName() { return "StrainDependentPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  virtual void saveConvergedState() const override;

  virtual void ignoreConvergedState() const override;

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
                          m_dPerm_dVolStrain,
                          m_referencePermeability,
                          m_volStrain,
                          m_volStrain_n,
                          m_permeabilityExponents );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * referencePermeabilityComponentsString() { return "referencePermeabilityComponents"; }
    static constexpr char const * permeabilityExponentsString() { return "permeabilityExponents"; }
    static constexpr char const * referencePermeabilityString() { return "referencePermeability"; }
    static constexpr char const * dPerm_dVolStrainString() { return "dPerm_dVolStrain"; }
    static constexpr char const * volumetricStrainString() { return "volumetricStrain"; }
    static constexpr char const * volumetricStrain_nString() { return "volumetricStrain_n"; }
    static constexpr char const * permeability_nString() { return "permeability_n"; }
  };

  virtual void initializeState() const override final;

protected:

  virtual void postInputInitialization() override;

private:

  /// Permeability components at the initial state (used where the permeability is not set by a field)
  R1Tensor m_referencePermeabilityComponents;

  /// Exponents of the porosity ratio for each permeability component
  R1Tensor m_permeabilityExponents;

  /// Permeability at the initial state
  array3d< real64 > m_referencePermeability;

  /// Derivative of the permeability w.r.t. the volumetric strain
  array3d< real64 > m_dPerm_dVolStrain;

  /// Total volumetric strain since the initial state
  array2d< real64 > m_volStrain;

  /// Total volumetric strain at the last converged state
  array2d< real64 > m_volStrain_n;

  /// Permeability at the last converged state
  array3d< real64 > m_permeability_n;
};


GEOS_HOST_DEVICE
inline
void StrainDependentPermeabilityUpdate::compute( real64 const & initialPorosity,
                                                 arraySlice1d< real64 const > const & referencePermeability,
                                                 real64 const & volStrain,
                                                 arraySlice1d< real64 > const & permeability,
                                                 arraySlice1d< real64 > const & dPerm_dVolStrain ) const
{
  // without (initial) porosity, the porosity ratio is undefined: keep the reference permeability
  if( initialPorosity <= minPorosity )
  {
    for( localIndex i = 0; i < permeability.size(); ++i )
    {
      permeability[i] = referencePermeability[i];
      dPerm_dVolStrain[i] = 0.0;
    }
    return;
  }

  // porosity from the conservation of the solid volume (incompressible grains),
  // clamped to a small positive value to remain defined under strong compaction
  real64 const solidFraction = ( 1.0 - initialPorosity ) * LvArray::math::exp( -volStrain );
  bool const isClamped = ( 1.0 - solidFraction ) < minPorosity;
  real64 const porosity = isClamped ? minPorosity : 1.0 - solidFraction;
  real64 const logPorosityRatio = LvArray::math::log( porosity / initialPorosity );

  for( localIndex i = 0; i < permeability.size(); ++i )
  {
    permeability[i] = referencePermeability[i] * LvArray::math::exp( m_permeabilityExponents[i] * logPorosityRatio );
    // dk/de_v = n k (1 - phi) / phi (zero where the porosity is clamped)
    dPerm_dVolStrain[i] = isClamped ? 0.0 : m_permeabilityExponents[i] * permeability[i] * solidFraction / porosity;
  }
}


}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_STRAINDEPENDENTPERMEABILITY_HPP_
