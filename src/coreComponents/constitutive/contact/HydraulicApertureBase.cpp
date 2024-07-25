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
 * @file HydraulicApertureBase.cpp
 */


#include "HydraulicApertureBase.hpp"


namespace geos
{

namespace constitutive
{

HydraulicApertureBase::HydraulicApertureBase( string const & name,
                                              Group * const parent ):
  ConstitutiveBase( name, parent )
{}

HydraulicApertureBase::~HydraulicApertureBase()
{}



} /* namespace constitutive */

} /* namespace geos */
