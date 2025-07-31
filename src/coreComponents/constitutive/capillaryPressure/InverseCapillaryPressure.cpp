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
 * @file InverseCapillaryPressure.cpp
 */

#include "InverseCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"

namespace geos
{

namespace constitutive
{

template< typename CAP_PRESSURE >
InverseCapillaryPressureUpdate< CAP_PRESSURE >::InverseCapillaryPressureUpdate( CAP_PRESSURE & capPressure,
                                                                                arrayView1d< real64 const > const & phaseMinVolumeFraction )
  : m_capPressureWrapper( capPressure.createKernelWrapper() ),
  m_phaseMinVolumeFraction( phaseMinVolumeFraction )
{
  m_sumMinVolumeFraction = 0.0;
  for( real64 const saturation : m_phaseMinVolumeFraction )
  {
    m_sumMinVolumeFraction += saturation;
  }
}

template< typename CAP_PRESSURE >
InverseCapillaryPressure< CAP_PRESSURE >::InverseCapillaryPressure( CAP_PRESSURE & capPressure )
  : m_capPressure( capPressure )
{}

template< typename CAP_PRESSURE >
typename InverseCapillaryPressure< CAP_PRESSURE >::KernelWrapper
InverseCapillaryPressure< CAP_PRESSURE >::createKernelWrapper()
{
  string const mivVolumeKey = CapillaryPressureBase::viewKeyStruct::phaseMinVolumeFractionString();
  array1d< real64 > & phaseMinVolumeFraction = m_capPressure.template getReference< array1d< real64 > >( mivVolumeKey );
  return KernelWrapper( m_capPressure, phaseMinVolumeFraction.toViewConst() );
}

template class InverseCapillaryPressure< BrooksCoreyCapillaryPressure >;
template class InverseCapillaryPressure< TableCapillaryPressure >;

} // namespace constitutive

} // namespace geos
