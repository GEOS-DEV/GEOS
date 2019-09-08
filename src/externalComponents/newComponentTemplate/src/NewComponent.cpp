/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * NewComponent.cpp
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#include "NewComponent.hpp"

namespace geosx
{

NewComponent::NewComponent( std::string const & name,
                            Group * const parent ):
    SolverBase(name,parent)
{

}

NewComponent::~NewComponent()
{}





real64 NewComponent::SolverStep( real64 const & /*time_n*/,
               real64 const & /*dt*/,
               integer const /*cycleNumber*/,
               DomainPartition * /*domain*/ )
{
  return 0;
}

REGISTER_CATALOG_ENTRY( SolverBase, NewComponent, std::string const &, dataRepository::Group * const )

} /* namespace geosx */
