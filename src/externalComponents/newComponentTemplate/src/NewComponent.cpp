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
 * @file NewComponent.cpp
 */

#include "NewComponent.hpp"

namespace geos
{

NewComponent::NewComponent( string const & name,
                            Group * const parent ):
    SolverBase(name,parent)
{

}

NewComponent::~NewComponent()
{}





real64 NewComponent::solverStep( real64 const & /*time_n*/,
               real64 const & /*dt*/,
               integer const /*cycleNumber*/,
               DomainPartition & /*domain*/ )
{
  return 0;
}

REGISTER_CATALOG_ENTRY( SolverBase, NewComponent, string const &, dataRepository::Group * const )

} /* namespace geos */
