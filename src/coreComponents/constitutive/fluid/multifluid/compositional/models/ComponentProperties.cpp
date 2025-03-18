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
 * @file ComponentProperties.cpp
 */

#include "ComponentProperties.hpp"
#include "ComponentType.hpp"

#include "common/format/StringUtilities.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

ComponentProperties::ComponentProperties( string_array const & componentNames,
                                          array1d< real64 > const & componentMolarWeight ):
  m_componentNames ( componentNames ),
  m_componentMolarWeight ( componentMolarWeight )
{
  classifyComponents( componentNames, m_componentType );
}

ComponentProperties::KernelWrapper::KernelWrapper(
  arrayView1d< integer const > const & componentType,
  arrayView1d< real64 const > const & componentMolarWeight,
  arrayView1d< real64 const > const & componentCriticalPressure,
  arrayView1d< real64 const > const & componentCriticalTemperature,
  arrayView1d< real64 const > const & componentAcentricFactor,
  arrayView1d< real64 const > const & componentVolumeShift,
  arrayView2d< real64 const > const & componentBinaryCoeff ):
  m_componentType ( componentType ),
  m_componentMolarWeight ( componentMolarWeight ),
  m_componentCriticalPressure ( componentCriticalPressure ),
  m_componentCriticalTemperature( componentCriticalTemperature ),
  m_componentAcentricFactor( componentAcentricFactor ),
  m_componentVolumeShift( componentVolumeShift ),
  m_componentBinaryCoeff( componentBinaryCoeff )
{}

ComponentProperties::KernelWrapper
ComponentProperties::createKernelWrapper() const
{
  return KernelWrapper( m_componentType,
                        m_componentMolarWeight,
                        m_componentCriticalPressure,
                        m_componentCriticalTemperature,
                        m_componentAcentricFactor,
                        m_componentVolumeShift,
                        m_componentBinaryCoeff );
}

void ComponentProperties::classifyComponents( string_array const & componentNames, array1d< integer > & componentType )
{
  integer const numComps = componentNames.size();
  componentType.resize( numComps );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    componentType[ic] = static_cast< integer >( getComponentTypeFromName( componentNames[ic] ) );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
