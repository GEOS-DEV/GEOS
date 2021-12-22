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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "constitutive/fluid/PVTFunctions/LohrenzBrayClarkViscosity.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

LohrenzBrayClarkViscosity::LohrenzBrayClarkViscosity( string const & name,
                                                      string_array const & componentNames,
                                                      array1d< real64 > const & componentMolarWeight,
                                                      array1d< real64 > const & componentCriticalPressure,
                                                      array1d< real64 > const & componentCriticalTemperature,
                                                      array1d< real64 > const & componentCriticalVolume,
                                                      array1d< real64 > const & componentAcentricFactor,
                                                      array1d< real64 > const & componentVolumeShift,
                                                      array2d< real64 > const & componentBinaryCoeff ):
  PVTCompositionalFunctionBase( name,
                                componentNames,
                                componentMolarWeight,
                                componentCriticalPressure,
                                componentCriticalTemperature,
                                componentCriticalVolume,
                                componentAcentricFactor,
                                componentVolumeShift,
                                componentBinaryCoeff )
{}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        m_componentCriticalPressure,
                        m_componentCriticalTemperature,
                        m_componentCriticalVolume,
                        m_componentAcentricFactor,
                        m_componentVolumeShift,
                        m_componentBinaryCoeff );
}

REGISTER_CATALOG_ENTRY( PVTCompositionalFunctionBase,
                        LohrenzBrayClarkViscosity, 
                        string const &, 
                        string_array const &, 
                        array1d< real64 > const &, 
                        array1d< real64 > const &, 
                        array1d< real64 > const &, 
                        array1d< real64 > const &, 
                        array1d< real64 > const &, 
                        array1d< real64 > const &, 
                        array2d< real64 > const & )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geosx

