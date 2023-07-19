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
 * @file SeismicityRateBase.cpp
 */

#include "SeismicityRateBase.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;

//START_SPHINX_INCLUDE_CONSTRUCTOR
SeismicityRateBase::SeismicityRateBase( const string & name,
                              Group * const parent ):
  SolverBase( name, parent ) {}
//END_SPHINX_INCLUDE_CONSTRUCTOR

SeismicityRateBase::~SeismicityRateBase()
{
  // TODO Auto-generated destructor stub
}

} // namespace geos
