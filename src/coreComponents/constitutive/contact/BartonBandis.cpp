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
 * @file BartonBandis.cpp
 */


#include "BartonBandis.hpp"


namespace geos
{

namespace constitutive
{

BartonBandis::BartonBandis( string const & name,
                                                Group * const parent ):
  HydraulicApertureBase( name, parent ),
  m_referenceNormalStress(0.0),
  m_aperture0(0.0)
{
  registerWrapper( viewKeyStruct::referenceNormalStressString(), &m_referenceNormalStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "" );

  registerWrapper( viewKeyStruct::apertureZeroString(), &m_aperture0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "" );
}  

BartonBandis::~BartonBandis()
{}



} /* namespace constitutive */

} /* namespace geos */