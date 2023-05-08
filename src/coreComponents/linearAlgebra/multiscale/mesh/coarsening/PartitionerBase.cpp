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

/**
 * @file PartitionerBase.cpp
 */

#include "PartitionerBase.hpp"

#include "linearAlgebra/multiscale/mesh/coarsening/CartesianPartitioner.hpp"
#include "linearAlgebra/multiscale/mesh/coarsening/GraphPartitioner.hpp"
#include "linearAlgebra/multiscale/mesh/coarsening/SemistructuredPartitioner.hpp"

namespace geos
{
namespace multiscale
{

std::unique_ptr< PartitionerBase >
PartitionerBase::create( LinearSolverParameters::Multiscale::Coarsening params )
{
  using PartitionType = LinearSolverParameters::Multiscale::Coarsening::PartitionType;
  switch( params.partitionType )
  {
    case PartitionType::graph:          return std::make_unique< GraphPartitioner >( std::move( params ) );
    case PartitionType::cartesian:      return std::make_unique< CartesianPartitioner >( std::move( params ) );
    case PartitionType::semistructured: return std::make_unique< SemistructuredPartitioner >( std::move( params ) );
    default:
    {
      GEOS_THROW( "Multiscale partitioning not supported yet: " << params.partitionType, std::runtime_error );
    }
  }
}

} // namespace multiscale
} // namespace geos
