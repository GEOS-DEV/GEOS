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
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Copyright (c) 2019, Lawrence Livermore National Security, LLC.
*
* Produced at the Lawrence Livermore National Laboratory
*
* LLNL-CODE-746361
*
* All rights reserved. See COPYRIGHT for details.
*
* This file is part of the GEOSX Simulation Framework.
*
* GEOSX is a free software; you can redistribute it and/or modify it under
* the terms of the GNU Lesser General Public License (as published by the
* Free Software Foundation) version 2.1 dated February 1999.
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
#include "AggregateElementSubRegion.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
AggregateElementSubRegion::AggregateElementSubRegion( string const & name,
dataRepository::Group * const parent):
ElementSubRegionBase( name, parent )
{
}

AggregateElementSubRegion::~AggregateElementSubRegion()
{
}

void AggregateElementSubRegion::CreateFromFineToCoarseMap( localIndex nbAggregates,
array1d< localIndex > const & fineToCoarse,
array1d< R1Tensor > const & barycenters)
{
m_elementCenter = barycenters;
m_nbFineCellsPerCoarseCell.resize( nbAggregates + 1 );
m_fineToCoarse.resize( fineToCoarse.size() );

/// First loop to count the number of fine cells per coarse cell
for( localIndex fineCell = 0; fineCell < fineToCoarse.size(); fineCell++ )
{
m_nbFineCellsPerCoarseCell[1 + fineToCoarse[fineCell]]++;
}

/// Second loop to cumulate the number of fine cells
for( localIndex coarseCell = 1; coarseCell <= nbAggregates; coarseCell++ )
{
m_nbFineCellsPerCoarseCell[coarseCell] += m_nbFineCellsPerCoarseCell[coarseCell-1];
}

/// Third loop to order the index of the fine cells
array1d< localIndex > offset( nbAggregates );
for( localIndex fineCell = 0; fineCell < fineToCoarse.size(); fineCell++ )
{
localIndex coarseCell = fineToCoarse[fineCell];
m_fineToCoarse[m_nbFineCellsPerCoarseCell[coarseCell] + offset[coarseCell]++] = fineCell;

}
}
}
