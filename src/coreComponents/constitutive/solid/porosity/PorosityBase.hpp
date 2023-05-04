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
 * @file PorosityBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
#define GEOS_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{
namespace constitutive
{

class PorosityBaseUpdates
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_newPorosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_newPorosity.size( 1 ); }

  PorosityBaseUpdates( arrayView2d< real64 > const & newPorosity,
                       arrayView2d< real64 > const & porosity_n,
                       arrayView2d< real64 > const & dPorosity_dPressure,
                       arrayView2d< real64 > const & dPorosity_dTemperature,
                       arrayView2d< real64 > const & initialPorosity,
                       arrayView1d< real64 > const & referencePorosity ):
    m_newPorosity( newPorosity ),
    m_porosity_n( porosity_n ),
    m_dPorosity_dPressure( dPorosity_dPressure ),
    m_dPorosity_dTemperature( dPorosity_dTemperature ),
    m_initialPorosity( initialPorosity ),
    m_referencePorosity ( referencePorosity )
  {}

  /**
   * @brief Helper to save point stress back to m_newPorosity array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] porosity porosity to be saved to m_newPorosity[k][q]
   * @param[in] dPorosity_dPressure porosity derivative w.r.t pressure to be saved to m_dPorosity_dPressure[k][q]
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void savePorosity( localIndex const k,
                     localIndex const q,
                     real64 const & porosity,
                     real64 const & dPorosity_dPressure ) const
  {
    m_newPorosity[k][q] = porosity;
    m_dPorosity_dPressure[k][q] = dPorosity_dPressure;
  }

  /**
   * @brief Helper to save point stress back to m_newPorosity array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] porosity porosity to be saved to m_newPorosity[k][q]
   * @param[in] dPorosity_dPressure porosity derivative w.r.t pressure to be saved to m_dPorosity_dPressure[k][q]
   * @param[in] dPorosity_dTemperature porosity derivative w.r.t temperature to be saved to m_dPorosity_dTemperature[k][q]
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void savePorosity( localIndex const k,
                     localIndex const q,
                     real64 const & porosity,
                     real64 const & dPorosity_dPressure,
                     real64 const & dPorosity_dTemperature ) const
  {
    m_newPorosity[k][q] = porosity;
    m_dPorosity_dPressure[k][q] = dPorosity_dPressure;
    m_dPorosity_dTemperature[k][q] = dPorosity_dTemperature;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 getPorosity( localIndex const k,
                      localIndex const q ) const
  {
    return m_newPorosity[k][q];
  }


  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 getPorosity_n( localIndex const k,
                        localIndex const q ) const
  {
    return m_porosity_n[k][q];
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 getInitialPorosity( localIndex const k,
                             localIndex const q ) const
  {
    return m_initialPorosity[k][q];
  }

  GEOS_HOST_DEVICE
  virtual void updateFromPressureAndTemperature( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & pressure,
                                                 real64 const & pressure_n,
                                                 real64 const & temperature,
                                                 real64 const & temperature_n ) const
  {
    GEOS_UNUSED_VAR( k, q, pressure, pressure_n, temperature, temperature_n );
    GEOS_ERROR( "updateFromPressureAndTemperature is not implemented for porosityBase." );
  }

protected:
  arrayView2d< real64 > m_newPorosity;

  arrayView2d< real64 > m_porosity_n;

  arrayView2d< real64 > m_dPorosity_dPressure;

  arrayView2d< real64 > m_dPorosity_dTemperature;

  arrayView2d< real64 > m_initialPorosity;

  arrayView1d< real64 > m_referencePorosity;
};


class PorosityBase : public ConstitutiveBase
{
public:
  PorosityBase( string const & name, Group * const parent );

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PorosityBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * defaultReferencePorosityString() { return "defaultReferencePorosity"; }
  } viewKeys;

  /**
   * @brief Number of elements storing solid data
   * @return Number of elements
   */
  localIndex numElem() const
  {
    return m_porosity_n.size( 0 );
  }

  /**
   * @brief Number of quadrature points storing solid data
   * @return Number of quadrature points
   */
  localIndex numQuad() const
  {
    return m_porosity_n.size( 1 );
  }


  /**
   * @brief Const accessor for newPorosity.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getPorosity() const { return m_newPorosity; }

  /**
   * @brief Const/non-mutable accessor for porosity_n.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getPorosity_n() const { return m_porosity_n; }


  /**
   * @brief Non-Const/mutable accessor for porosity_n
   * @return Accessor
   */
  arrayView2d< real64 > const getPorosity_n() { return m_porosity_n; }


  /**
   * @brief Const/non-mutable accessor for referencePorosity.
   * @return Accessor
   */
  arrayView1d< real64 const > const  getReferencePorosity() const { return m_referencePorosity; }


  /**
   * @brief Const/non-mutable accessor for dPorosity_dPressure
   * @return Accessor
   */
  arrayView2d< real64 const > const  dPorosity_dPressure() const { return m_dPorosity_dPressure; }

  /**
   * @brief Const/non-mutable accessor for dPorosity_dTemperature
   * @return Accessor
   */
  arrayView2d< real64 const > const  dPorosity_dTemperature() const { return m_dPorosity_dTemperature; }

  /**
   * @brief Utility function to scale the reference porosity (for instance, by net-to-gross)
   * @param[in] scalingFactors the vector of scaling factors (one value per cell) for the reference porosity
   */
  void scaleReferencePorosity( arrayView1d< real64 const > scalingFactors ) const;

  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const override;

  /**
   * @brief Initialize newPorosity and porosity_n.
   */
  virtual void initializeState() const;

  virtual arrayView1d< real64 const > const getBiotCoefficient() const
  {
    GEOS_ERROR( "getBiotPorosity() not implemented for this model" );

    array1d< real64 > out;
    return out.toViewConst();
  }

  using KernelWrapper = PorosityBaseUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_newPorosity,
                          m_porosity_n,
                          m_dPorosity_dPressure,
                          m_dPorosity_dTemperature,
                          m_initialPorosity,
                          m_referencePorosity );
  }


protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_newPorosity;

  array2d< real64 > m_porosity_n;

  array2d< real64 > m_dPorosity_dPressure;

  array2d< real64 > m_dPorosity_dTemperature;

  array2d< real64 > m_initialPorosity;

  array1d< real64 > m_referencePorosity;

  real64 m_defaultReferencePorosity;

};

} /* namespace constitutive */

} /* namespace geos */


#endif //GEOS_CONSTITUTIVE_POROSITY_POROSITYBASE_HPP_
