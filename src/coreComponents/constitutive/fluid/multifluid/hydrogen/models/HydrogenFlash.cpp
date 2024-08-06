/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrogenFlash.cpp
 */

#include "HydrogenFlash.hpp"

namespace geos
{
namespace constitutive
{

HydrogenFlashUpdate::HydrogenFlashUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                                          integer const h2ComponentIndex,
                                          integer const h2oComponentIndex,
                                          integer const gasPhaseIndex,
                                          integer const watPhaseIndex ):
  m_numComps( componentMolarWeight.size() ),
  m_h2Index( h2ComponentIndex ),
  m_h2oIndex( h2oComponentIndex ),
  m_gasPhaseIndex( gasPhaseIndex ),
  m_watPhaseIndex( watPhaseIndex )
{
  GEOS_UNUSED_VAR( componentMolarWeight );
}

void HydrogenFlashUpdate::move( LvArray::MemorySpace const space, bool const touch )
{
  GEOS_UNUSED_VAR( space );
  GEOS_UNUSED_VAR( touch );
}

HydrogenFlash::HydrogenFlash( string const & name,
                              arrayView1d< real64 > const & componentMolarWeight,
                              integer const h2ComponentIndex,
                              integer const h2oComponentIndex,
                              integer const gasPhaseIndex,
                              integer const watPhaseIndex,
                              bool const printTable ):
  m_componentMolarWeight( componentMolarWeight ),
  m_h2Index( h2ComponentIndex ),
  m_h2oIndex( h2oComponentIndex ),
  m_gasPhaseIndex( gasPhaseIndex ),
  m_watPhaseIndex( watPhaseIndex )
{
  GEOS_UNUSED_VAR( name );
  GEOS_UNUSED_VAR( printTable );
}

HydrogenFlash::KernelWrapper HydrogenFlash::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight.toViewConst(),
                        m_h2Index,
                        m_h2oIndex,
                        m_gasPhaseIndex,
                        m_watPhaseIndex );
}

} // namespace constitutive

} // end namespace geos
