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
 * @file PowerLawPermeability.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITY_POWERLAWPERMEABILITY_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITY_POWERLAWPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geos
{
namespace constitutive
{

/**
 * @brief Update class for the power law permeability model.
 *
 * The permeability is scaled from its reference value with the pore space available to the fluid
 *
 *   k = k_min + ( k_0 - k_min ) * f( phi, phi_0, phi_c, n ),   f = ( ( phi - phi_c ) / ( phi_0 - phi_c ) )^n
 *
 * where phi is the porosity, phi_0 the reference porosity at which the permeability is k_0,
 * phi_c the critical (percolation) porosity below which the pore space no longer percolates,
 * and n the exponent of the power law. The residual permeability k_min (zero by default) is the
 * floor the permeability decays to once the pore space stops percolating.
 *
 * @see PFLOTRAN theory guide, reactive transport mode, "Changes in material properties" (Eq. 30-31)
 * https://documentation.pflotran.org/theory_guide/mode_reactive_transport.html#changes-in-material-properties
 */
class PowerLawPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  PowerLawPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & dPerm_dPressure,
                              arrayView3d< real64 > const & dPerm_dPorosity,
                              arrayView3d< real64 const > const & referencePermeability,
                              real64 const referencePorosity,
                              real64 const criticalPorosity,
                              real64 const exponent,
                              real64 const minPermeability )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dPorosity( dPerm_dPorosity ),
    m_referencePermeability( referencePermeability ),
    m_referencePorosity( referencePorosity ),
    m_criticalPorosity( criticalPorosity ),
    m_exponent( exponent ),
    m_minPermeability( minPermeability )
  {}

  GEOS_HOST_DEVICE
  void compute( real64 const & porosity,
                arraySlice1d< real64 const > const & referencePermeability,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPorosity ) const;

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndPorosity( localIndex const k,
                                              localIndex const q,
                                              real64 const & pressure,
                                              real64 const & porosity ) const override
  {
    GEOS_UNUSED_VAR( pressure );

    compute( porosity,
             m_referencePermeability[k][0],
             m_permeability[k][q],
             m_dPerm_dPorosity[k][q] );
  }

private:

  /// dPermeability_dPorosity
  arrayView3d< real64 > m_dPerm_dPorosity;

  /// Permeability at the reference porosity
  arrayView3d< real64 const > m_referencePermeability;

  /// Porosity at which the permeability is the reference permeability
  real64 m_referencePorosity;

  /// Porosity below which the pore space no longer percolates
  real64 m_criticalPorosity;

  /// Exponent of the power law
  real64 m_exponent;

  /// Residual permeability reached when the pore space stops percolating
  real64 m_minPermeability;
};


/**
 * @brief Model in which the permeability follows a power law of the porosity.
 */
class PowerLawPermeability : public PermeabilityBase
{
public:

  PowerLawPermeability( string const & name, dataRepository::Group * const parent );

  static string catalogName() { return "PowerLawPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numPts ) override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PowerLawPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dPorosity,
                          m_referencePermeability,
                          m_referencePorosity,
                          m_criticalPorosity,
                          m_exponent,
                          m_minPermeability );
  }

  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dPorosityString() { return "dPerm_dPorosity"; }
    static constexpr char const * referencePermeabilityComponentsString() { return "referencePermeabilityComponents"; }
    static constexpr char const * referencePermeabilityString() { return "referencePermeability"; }
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * criticalPorosityString() { return "criticalPorosity"; }
    static constexpr char const * exponentString() { return "exponent"; }
    static constexpr char const * minPermeabilityString() { return "minPermeability"; }
  };

  virtual void initializeState() const override final;

protected:

  virtual void postInputInitialization() override;

private:

  /// dPermeability_dPorosity
  array3d< real64 > m_dPerm_dPorosity;

  /// Permeability components at the reference porosity
  R1Tensor m_referencePermeabilityComponents;

  /// Permeability field at the reference porosity
  array3d< real64 > m_referencePermeability;

  /// Porosity at which the permeability is the reference permeability
  real64 m_referencePorosity;

  /// Porosity below which the pore space no longer percolates
  real64 m_criticalPorosity;

  /// Exponent of the power law
  real64 m_exponent;

  /// Residual permeability reached when the pore space stops percolating
  real64 m_minPermeability;
};


GEOS_HOST_DEVICE
inline
void PowerLawPermeabilityUpdate::compute( real64 const & porosity,
                                          arraySlice1d< real64 const > const & referencePermeability,
                                          arraySlice1d< real64 > const & permeability,
                                          arraySlice1d< real64 > const & dPerm_dPorosity ) const
{
  real64 const connectedPorosity = porosity - m_criticalPorosity;

  // below the percolation threshold the pore space is disconnected and only the residual permeability is left
  real64 const permMultiplier = ( connectedPorosity > 0.0 ) ?
                                pow( connectedPorosity / ( m_referencePorosity - m_criticalPorosity ), m_exponent ) : 0.0;

  real64 const dPermMultiplier_dPorosity = ( connectedPorosity > 0.0 ) ?
                                           m_exponent * permMultiplier / connectedPorosity : 0.0;

  for( localIndex i = 0; i < permeability.size(); ++i )
  {
    real64 const permRange = referencePermeability[i] - m_minPermeability;

    permeability[i] = m_minPermeability + permRange * permMultiplier;
    dPerm_dPorosity[i] = permRange * dPermMultiplier_dPorosity;
  }
}

}/* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_PERMEABILITY_POWERLAWPERMEABILITY_HPP_
