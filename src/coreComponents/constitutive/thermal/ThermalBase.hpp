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
 * @file ContinuumBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_THERMALBASE_HPP_
#define GEOS_CONSTITUTIVE_THERMALBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Base class for all thermal constitutive kernel wrapper classes.
 *
 * The responsibility of this base is to:
 *
 * 1) Contain views to state and parameter data for solid models.
 * 2) Specify an interface for state update functions.
 *
 * In general, the ArrayView data in the wrapper is specified to be of type
 * "arrayView<T> const" or "arrayView<T const> const". The "const-ness"
 * of the data indicates whether it is a parameter" or a state variable,
 * with the parameters being "T const" and state variables being "T".
 *
 * @note
 * If an allocation occurs on the underlying Array after a KernelWrapper is created,
 * then the ArrayView members of that KernelWrapper are silently invalid.
 */
class ThermalBaseUpdates
{
protected:
  /**
   * @brief constructor
   * @param[in] newTemperature The new temperature data from the constitutive model class.
   * @param[in] oldTemperature The old temperature data from the constitutive model class.
   * @param[in] heatCapacity The heat capacity data for the constitutive mode class.
   */
  ThermalBaseUpdates( arrayView2d< real64 > const & newTemperature,
                      arrayView2d< real64 > const & oldTemperature,
                      arrayView2d< real64 > const & heatCapacity ):
    m_newTemperature( newTemperature ),
    m_oldTemperature( oldTemperature ),
    m_heatCapacity( heatCapacity ),
  {}

  /// Deleted default constructor
  ThermalBaseUpdates() = delete;

  /**
   * @brief Copy Constructor
   * @param source Object to copy
   */
  ThermalBaseUpdates( ThermalBaseUpdates const & source ) = default;

  /**
   * @brief Move Constructor
   * @param source Object to move resources from
   */
  ThermalBaseUpdates( ThermalBaseUpdates && source ) = default;

  /// Deleted copy assignment operator
  ThermalBaseUpdates & operator=( ThermalBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  ThermalBaseUpdates & operator=( ThermalBaseUpdates && ) =  delete;

public:

  /// A reference the current material temperature at quadrature points.
  arrayView3d< real64 > const m_newTemperature;

  /// A reference the previous material temperature at quadrature points.
  arrayView3d< real64 > const m_oldTemperature;

  /// A reference the current material heat capacity at quadrature points.
  arrayView2d< real64 > m_heatCapacity;

  /**
   * @name Update Interfaces: Stress and Stiffness
   *
   * We define a variety of interfaces for constitutive models using different
   * strain theories.  Derived classes only need to implement the subset of interfaces
   * most relevant to them.
   *
   * This group of interfaces returns stress and stiffness simultaneously, and
   * are most useful for implicit finite element formulations.
   */
  ///@{

  /**
   * @brief Incremental thermal update.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] timeIncrement time increment for rate-dependent models.
   * @param[out] temperature New temperature value
   */
  GEOS_HOST_DEVICE
  /**
   * this function is not virtual to avoid a compilation warning with nvcc.
   */
  void incrementalThermalUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement
                          real64 & temperature ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( timeIncrement );
    GEOS_UNUSED_VAR( temperature );
    GEOS_ERROR( "incrementalThermalUpdate() not implemented for this model" );
  }


  /**
   * @brief Helper to save point stress back to m_newStress array
   *
   * This is mostly defined for improving code readability.
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   * @param[in] stress Stress to be save to m_newStress[k][q]
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void saveTemperature( localIndex const k,
                        localIndex const q,
                        real64 const  & temperature ) const
  {
    LvArray::tensorOps::copy< 6 >( m_newTemperature[k][q], temperature );
  }

  /**
   * @brief Save converged state data at index (k,q)
   *
   * @param[in] k Element index.
   * @param[in] q Quadrature point index.
   */
  GEOS_HOST_DEVICE
  inline
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const
  {
    LvArray::tensorOps::copy< 6 >( m_oldTemperature[k][q], m_newTemperature[k][q] );
  }
};


/**
 * @class ThermalBase
 * This class serves as the base class for thermal constitutive models.
 */
class ThermalBase : public constitutive::ConstitutiveBase
{
public:
  /**
   * @brief Constructor
   * @param name Name of the ThermalBase object in the repository.
   * @param parent The parent group of the ThermalBase object.
   */
  ThermalBase( string const & name,
                 Group * const parent );

  /**
   * Destructor
   */
  virtual ~ThermalBase() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "ThermalBase";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// Keys for data in this class
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * temperatureString() { return "temperature"; }                  ///< New temperature key
    static constexpr char const * oldTemperatureString() { return "oldTemperature"; }            ///< Old temperature key
    static constexpr char const * heatCapacityString() { return "heatCapacity"; }                ///< Heat capacity key
    static constexpr char const * defaultHeatCapacityString() { return "defaultHeatCapacity"; }  ///< Default heat capacity key
  };

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (element subregion)
   * @param numConstitutivePointsPerParentIndex Number of quadrature points per element
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /// Save state data in preparation for next timestep
  virtual void saveConvergedState() const override;

  /**
   * @brief Number of elements storing solid data
   * @return Number of elements
   */
  localIndex numElem() const
  {
    return m_oldTemperature.size( 0 );
  }

  /**
   * @brief Number of quadrature points storing solid data
   * @return Number of quadrature points
   */
  localIndex numQuad() const
  {
    return m_oldTemperature.size( 1 );
  }

  /**
   * @name Accessors
   */
  ///@{

  /**
   * @brief Non-const/mutable accessor for temperature
   * @return Accessor
   */
  arrayView3d< real64 > const getTemperature()
  {
    return m_newTemperature;
  }

  /**
   * @brief Const/non-mutable accessor for temperature
   * @return Accessor
   */
  arrayView3d< real64 const > const getTemperature() const
  {
    return m_newTemperature;
  }

  /**
   * @brief Non-const/mutable accessor for old temperature
   * @return Accessor
   */
  arrayView3d< real64 > const getOldTemperature()
  {
    return m_oldTemperature;
  }

  /**
   * @brief Const/non-mutable accessor for old temperature
   * @return Accessor
   */
  arrayView3d< real64 const> const getOldTemperature() const
  {
    return m_oldTemperature;
  }

  /**
   * @brief Non-const/Mutable accessor for heat capacity.
   * @return Accessor
   */
  arrayView2d< real64 > const getHeatCapacity()
  {
    return m_heatCapacity;
  }

  /**
   * @brief Const/non-mutable accessor for heat capacity
   * @return Accessor
   */
  arrayView2d< real64 const > const getgetHeatCapacity() const
  {
    return m_heatCapacity;
  }

  ///@}
  //

protected:

  /// Post-process XML input
  virtual void postInputInitialization() override;

  /// The current temperature at a quadrature point (i.e. at timestep n, global newton iteration k)
  array3d< real64, solid::STRESS_PERMUTATION > m_newTemperature;

  /// The previous temperature at a quadrature point (i.e. at timestep (n-1))
  array3d< real64, solid::STRESS_PERMUTATION > m_oldTemperature;

  /// The material heat capacity at a quadrature point.
  array2d< real64 > m_heatCapacity;

  /// The default heat capacity for new allocations.
  real64 m_defaultHeatCapacity;

  /// band-aid fix...going to have to remove this after we clean up
  /// initialization for constitutive models.
  bool m_postProcessed = false;
};

} // namespace constitutive
} // namespace geos

#endif /* GEOS_CONSTITUTIVE_THERMALBASE_HPP_ */
