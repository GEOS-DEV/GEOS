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
 * @file SolidThermalExpansion.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDTHERMALEXPANSION_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDTHERMALEXPANSION_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class SolidThermalExpansionUpdates
 *
 * @brief Device-side evaluator of the drained linear thermal expansion coefficient alpha(T).
 *
 * The evaluation strategy and all the data it needs (scalars or array views) are fixed on the
 * host when the object is constructed. On device, compute() only branches on an enum: there is
 * no string, no pointer test, no virtual call and no table object.
 *
 * Every supported law is piecewise linear in temperature, so the mean value of alpha over an
 * interval is computed exactly (the mean of a linear function over an interval is its value at
 * the midpoint).
 */
class SolidThermalExpansionUpdates
{
public:

  /// Evaluation strategy, chosen on the host.
  enum class Strategy : integer
  {
    FromSolidBase, ///< per-element constant alpha read from the solid model (legacy input)
    Linear,        ///< alpha(T) = alpha0 + dalpha/dT * ( T - T0 )
    Tabulated      ///< alpha(T) piecewise linear between table knots, constant outside
  };

  /**
   * @brief Build a wrapper reading a per-element constant coefficient (legacy input).
   * @param[in] thermalExpansionCoefficient per-element thermal expansion coefficient [1/temperature]
   */
  explicit SolidThermalExpansionUpdates( arrayView1d< real64 const > const & thermalExpansionCoefficient ):
    m_strategy( Strategy::FromSolidBase ),
    m_thermalExpansionCoefficient( thermalExpansionCoefficient )
  {}

  /**
   * @brief Build a wrapper for the linear law.
   * @param[in] referenceThermalExpansionCoefficient coefficient at the reference temperature [1/temperature]
   * @param[in] dThermalExpansionCoefficient_dTemperature derivative of the coefficient with respect to temperature [1/temperature^2]
   * @param[in] referenceTemperature temperature at which the coefficient equals @p referenceThermalExpansionCoefficient [temperature]
   */
  SolidThermalExpansionUpdates( real64 const referenceThermalExpansionCoefficient,
                                real64 const dThermalExpansionCoefficient_dTemperature,
                                real64 const referenceTemperature ):
    m_strategy( Strategy::Linear ),
    m_referenceThermalExpansionCoefficient( referenceThermalExpansionCoefficient ),
    m_dThermalExpansionCoefficient_dTemperature( dThermalExpansionCoefficient_dTemperature ),
    m_referenceTemperature( referenceTemperature )
  {}

  /**
   * @brief Build a wrapper for the tabulated law.
   * @param[in] knots strictly increasing table temperatures [temperature]
   * @param[in] values coefficient at each knot [1/temperature]
   * @param[in] cumulativeIntegral integral of the coefficient from the first knot to each knot [dimensionless]
   */
  SolidThermalExpansionUpdates( arrayView1d< real64 const > const & knots,
                                arrayView1d< real64 const > const & values,
                                arrayView1d< real64 const > const & cumulativeIntegral ):
    m_strategy( Strategy::Tabulated ),
    m_knots( knots ),
    m_values( values ),
    m_cumulativeIntegral( cumulativeIntegral )
  {}

  /**
   * @brief Evaluate the thermal expansion coefficient over a temperature interval.
   *
   * The thermal strain increment over [@p temperatureLower, @p temperatureUpper] is
   * meanThermalExpansionCoefficient * ( temperatureUpper - temperatureLower ), i.e. the exact integral of alpha.
   * Its derivative with respect to @p temperatureUpper is alpha( temperatureUpper ).
   *
   * @param[in] k element index (used by the FromSolidBase strategy only)
   * @param[in] temperatureLower lower bound of the interval (typically the temperature at the previous time step) [temperature]
   * @param[in] temperatureUpper upper bound of the interval (typically the current temperature) [temperature]
   * @param[out] meanThermalExpansionCoefficient mean value of alpha over the interval [1/temperature];
   *             equal to alpha when the bounds coincide
   * @param[out] thermalExpansionCoefficient value of alpha at @p temperatureUpper [1/temperature]
   */
  GEOS_HOST_DEVICE
  inline
  void compute( localIndex const k,
                real64 const temperatureLower,
                real64 const temperatureUpper,
                real64 & meanThermalExpansionCoefficient,
                real64 & thermalExpansionCoefficient ) const
  {
    switch( m_strategy )
    {
      case Strategy::FromSolidBase:
      {
        thermalExpansionCoefficient = m_thermalExpansionCoefficient[k];
        meanThermalExpansionCoefficient = thermalExpansionCoefficient;
        break;
      }
      case Strategy::Linear:
      {
        thermalExpansionCoefficient = m_referenceThermalExpansionCoefficient + m_dThermalExpansionCoefficient_dTemperature * ( temperatureUpper - m_referenceTemperature );
        meanThermalExpansionCoefficient = m_referenceThermalExpansionCoefficient + m_dThermalExpansionCoefficient_dTemperature *
                                          ( 0.5 * ( temperatureLower + temperatureUpper ) - m_referenceTemperature );
        break;
      }
      case Strategy::Tabulated:
      {
        computeTabulated( temperatureLower, temperatureUpper, meanThermalExpansionCoefficient, thermalExpansionCoefficient );
        break;
      }
    }
  }

  /**
   * @return the evaluation strategy
   */
  GEOS_HOST_DEVICE
  inline
  Strategy getStrategy() const { return m_strategy; }

private:

  /**
   * @brief Index of the table segment containing a temperature.
   * @param[in] temperature the temperature
   * @return -1 below the first knot, numKnots-1 at or above the last knot, otherwise i such that
   *         knots[i] <= temperature < knots[i+1]
   */
  GEOS_HOST_DEVICE
  inline
  localIndex findSegment( real64 const temperature ) const
  {
    localIndex const numKnots = m_knots.size();
    if( temperature < m_knots[0] )
    {
      return -1;
    }
    if( temperature >= m_knots[numKnots - 1] )
    {
      return numKnots - 1;
    }
    localIndex lower = 0;
    localIndex upper = numKnots - 1;
    while( upper - lower > 1 )
    {
      localIndex const middle = ( lower + upper ) / 2;
      if( m_knots[middle] <= temperature )
      {
        lower = middle;
      }
      else
      {
        upper = middle;
      }
    }
    return lower;
  }

  /**
   * @brief Value of alpha at a temperature located in a given segment.
   * @param[in] segment segment index returned by findSegment()
   * @param[in] temperature the temperature
   * @return alpha( temperature )
   */
  GEOS_HOST_DEVICE
  inline
  real64 coefficientInSegment( localIndex const segment, real64 const temperature ) const
  {
    localIndex const numKnots = m_knots.size();
    if( segment < 0 )
    {
      return m_values[0];
    }
    if( segment >= numKnots - 1 )
    {
      return m_values[numKnots - 1];
    }
    real64 const slope = ( m_values[segment + 1] - m_values[segment] ) / ( m_knots[segment + 1] - m_knots[segment] );
    return m_values[segment] + slope * ( temperature - m_knots[segment] );
  }

  /**
   * @brief Exact mean value of the tabulated alpha over an interval, and its value at the upper bound.
   * @param[in] temperatureLower lower bound
   * @param[in] temperatureUpper upper bound
   * @param[out] meanThermalExpansionCoefficient mean value of alpha over the interval
   * @param[out] thermalExpansionCoefficient value of alpha at @p temperatureUpper
   */
  GEOS_HOST_DEVICE
  inline
  void computeTabulated( real64 const temperatureLower,
                         real64 const temperatureUpper,
                         real64 & meanThermalExpansionCoefficient,
                         real64 & thermalExpansionCoefficient ) const
  {
    thermalExpansionCoefficient = coefficientInSegment( findSegment( temperatureUpper ), temperatureUpper );

    // The mean value does not depend on the orientation of the interval.
    real64 const a = ( temperatureLower < temperatureUpper ) ? temperatureLower : temperatureUpper;
    real64 const b = ( temperatureLower < temperatureUpper ) ? temperatureUpper : temperatureLower;
    localIndex const segmentA = findSegment( a );
    localIndex const segmentB = findSegment( b );

    if( segmentA == segmentB )
    {
      // alpha is linear on [a, b]: its mean value is its value at the midpoint.
      // This branch also covers a == b without any division.
      meanThermalExpansionCoefficient = coefficientInSegment( segmentA, 0.5 * ( a + b ) );
      return;
    }

    // [a, b] spans at least one knot: partial first segment + whole segments + partial last segment.
    real64 const firstKnot = m_knots[segmentA + 1];
    real64 const lastKnot = m_knots[segmentB];
    real64 const integral = ( firstKnot - a ) * coefficientInSegment( segmentA, 0.5 * ( a + firstKnot ) )
                            + ( m_cumulativeIntegral[segmentB] - m_cumulativeIntegral[segmentA + 1] )
                            + ( b - lastKnot ) * coefficientInSegment( segmentB, 0.5 * ( lastKnot + b ) );
    meanThermalExpansionCoefficient = integral / ( b - a );
  }

  /// Evaluation strategy
  Strategy m_strategy = Strategy::Linear;

  /// Coefficient at the reference temperature (Linear)
  real64 m_referenceThermalExpansionCoefficient = 0.0;

  /// Derivative of the coefficient with respect to temperature (Linear)
  real64 m_dThermalExpansionCoefficient_dTemperature = 0.0;

  /// Temperature at which the coefficient equals m_referenceThermalExpansionCoefficient (Linear)
  real64 m_referenceTemperature = 0.0;

  /// Per-element coefficient (FromSolidBase)
  arrayView1d< real64 const > m_thermalExpansionCoefficient;

  /// Table temperatures (Tabulated)
  arrayView1d< real64 const > m_knots;

  /// Coefficient at each knot (Tabulated)
  arrayView1d< real64 const > m_values;

  /// Integral of the coefficient from the first knot to each knot (Tabulated)
  arrayView1d< real64 const > m_cumulativeIntegral;
};


/**
 * @class SolidThermalExpansion
 *
 * @brief Drained linear thermal expansion coefficient alpha(T), either linear in temperature or
 *        given by a one-dimensional TableFunction of temperature.
 *
 * The model is attached to a coupled solid (e.g. PorousElasticIsotropic) through the optional
 * attribute solidThermalExpansionModelName. When it is not attached, the coupled solid keeps
 * using the constant per-element coefficient of the solid model.
 */
class SolidThermalExpansion : public ConstitutiveBase
{
public:

  /// Kernel wrapper type
  using KernelWrapper = SolidThermalExpansionUpdates;

  /**
   * @brief Constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  SolidThermalExpansion( string const & name, dataRepository::Group * const parent );

  /**
   * @return the catalog name
   */
  static string catalogName() { return "SolidThermalExpansion"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create the kernel wrapper. Must be called on the host.
   * @return the kernel wrapper
   */
  KernelWrapper createKernelUpdates() const;

  /**
   * @return true if the coefficient depends on temperature (table, or non-zero derivative)
   */
  bool isTemperatureDependent() const;

  /// Keys of the data repository wrappers
  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    /// key of the coefficient at the reference temperature
    static constexpr char const * referenceDrainedLinearTECString() { return "referenceDrainedLinearTEC"; }
    /// key of the derivative of the coefficient with respect to temperature
    static constexpr char const * dDrainedLinearTEC_dTemperatureString() { return "dDrainedLinearTEC_dTemperature"; }
    /// key of the reference temperature
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
    /// key of the table name
    static constexpr char const * drainedLinearTECTableNameString() { return "drainedLinearTECTableName"; }
    /// key of the table temperatures copied from the table
    static constexpr char const * tableKnotsString() { return "tableKnots"; }
    /// key of the table values copied from the table
    static constexpr char const * tableValuesString() { return "tableValues"; }
    /// key of the cumulative integral of the table
    static constexpr char const * tableCumulativeIntegralString() { return "tableCumulativeIntegral"; }
  };

protected:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

private:

  /// Coefficient at the reference temperature
  real64 m_referenceDrainedLinearTEC;

  /// Derivative of the coefficient with respect to temperature
  real64 m_dDrainedLinearTEC_dTemperature;

  /// Temperature at which the coefficient equals m_referenceDrainedLinearTEC
  real64 m_referenceTemperature;

  /// Name of the TableFunction giving the coefficient as a function of temperature
  string m_drainedLinearTECTableName;

  /// Table temperatures, copied from the TableFunction during initialization
  array1d< real64 > m_tableKnots;

  /// Coefficient at each knot, copied from the TableFunction during initialization
  array1d< real64 > m_tableValues;

  /// Integral of the coefficient from the first knot to each knot
  array1d< real64 > m_tableCumulativeIntegral;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDTHERMALEXPANSION_HPP_ */
