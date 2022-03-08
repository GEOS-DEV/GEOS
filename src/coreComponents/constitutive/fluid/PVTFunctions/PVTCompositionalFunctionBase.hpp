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
 * @file PVTCompositionalFunctionBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTCOMPOSITIONALFUNCTIONBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTCOMPOSITIONALFUNCTIONBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

enum class PVTCompositionalFunctionType { UNKNOWN, DENSITY, VISCOSITY, ENTHALPY, INTERNAL_ENERGY};

class PVTCompositionalFunctionBaseUpdate
{
public:

  PVTCompositionalFunctionBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                                      arrayView1d< real64 const > const & componentCriticalPressure,
                                      arrayView1d< real64 const > const & componentCriticalTemperature,
                                      arrayView1d< real64 const > const & componentCriticalVolume,
                                      arrayView1d< real64 const > const & componentAcentricFactor,
                                      arrayView1d< real64 const > const & componentVolumeShift,
                                      arrayView2d< real64 const > const & componentBinaryCoeff )
    : m_componentMolarWeight( componentMolarWeight ),
      m_componentCriticalPressure ( componentCriticalPressure ),
      m_componentCriticalTemperature ( componentCriticalTemperature ),
      m_componentCriticalVolume ( componentCriticalVolume ),
      m_componentAcentricFactor ( componentAcentricFactor ),
      m_componentVolumeShift ( componentVolumeShift ),
      m_componentBinaryCoeff ( componentBinaryCoeff )
  {}

  /**
   * @brief Move the KernelWrapper to the given execution space, optionally touching it.
   * @param space the space to move the KernelWrapper to
   * @param touch whether the KernelWrapper should be touched in the new space or not
   * @note This function exists to enable holding KernelWrapper objects in an ArrayView
   *       and have their contents properly moved between memory spaces.
   */
  virtual void move( LvArray::MemorySpace const space, bool const touch )
  {
    m_componentMolarWeight.move( space, touch );
    m_componentCriticalPressure.move( space, touch );
    m_componentCriticalTemperature.move( space, touch );
    m_componentCriticalVolume.move( space, touch );
    m_componentAcentricFactor.move( space, touch );
    m_componentVolumeShift.move( space, touch );
    m_componentBinaryCoeff.move( space, touch );
  }

protected:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using FluidProp = MultiFluidVar< real64, 2, multifluid::LAYOUT_FLUID, multifluid::LAYOUT_FLUID_DC >;

  /// ArrayView storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

  /// ArrayView storing the component critical pressures
  arrayView1d< real64 const > m_componentCriticalPressure;

  /// ArrayView storing the component critical temperatures
  arrayView1d< real64 const > m_componentCriticalTemperature;

  /// ArrayView storing the component critical volumes
  arrayView1d< real64 const > m_componentCriticalVolume;

  /// ArrayView storing the component acentric factor
  arrayView1d< real64 const > m_componentAcentricFactor;

  /// ArrayView storing the component volume shift
  arrayView1d< real64 const > m_componentVolumeShift;
  
  /// ArrayView storing the component binary coefficient
  arrayView2d< real64 const > m_componentBinaryCoeff;
};

class PVTCompositionalFunctionBase
{

public:

  PVTCompositionalFunctionBase( string const & name,
                                string_array const & componentNames,
                                array1d< real64 > const & componentMolarWeight,
                                array1d< real64 > const & componentCriticalPressure,
                                array1d< real64 > const & componentCriticalTemperature,
                                array1d< real64 > const & componentCriticalVolume,
                                array1d< real64 > const & componentAcentricFactor,
                                array1d< real64 > const & componentVolumeShift,
                                array2d< real64 > const & componentBinaryCoeff )
    :
    m_functionName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight ),
    m_componentCriticalPressure( componentCriticalPressure ),
    m_componentCriticalTemperature( componentCriticalTemperature ),
    m_componentCriticalVolume( componentCriticalVolume ),
    m_componentAcentricFactor( componentAcentricFactor ),
    m_componentVolumeShift( componentVolumeShift ),
    m_componentBinaryCoeff( componentBinaryCoeff )
  {}

  virtual ~PVTCompositionalFunctionBase() = default;

  using CatalogInterface = dataRepository::CatalogInterface< PVTCompositionalFunctionBase,
                                                             string const &,
                                                             string_array const &,
                                                             array1d< real64 > const &,
                                                             array1d< real64 > const &,
                                                             array1d< real64 > const &,
                                                             array1d< real64 > const &,
                                                             array1d< real64 > const &,
                                                             array1d< real64 > const &,
                                                             array2d< real64 > const & >;

  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() const = 0;

  string const & functionName() const { return m_functionName; }

  virtual PVTCompositionalFunctionType functionType() const = 0;

protected:

  /// Name of the PVT function
  string m_functionName;

  /// Array storing the name of the components
  array1d< string > m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

  /// Array storing the component critical pressures
  array1d< real64 > m_componentCriticalPressure;

  /// Array storing the component critical temperatures
  array1d< real64 > m_componentCriticalTemperature;

  /// Array storing the component critical temperatures
  array1d< real64 > m_componentCriticalVolume;

  /// Array storing the component acentric factor
  array1d< real64 > m_componentAcentricFactor;

  /// Array storing the component volume shift
  array1d< real64 > m_componentVolumeShift;
  
  /// Array storing the component binary coefficient
  array2d< real64 > m_componentBinaryCoeff;
};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTCOMPOSITIONALFUNCTIONBASE_HPP_
