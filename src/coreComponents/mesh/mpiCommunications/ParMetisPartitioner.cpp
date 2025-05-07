
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
 * @file ParMetisPartitioner.cpp
 */

#include "ParMetisPartitioner.hpp"

void ParMetisPartitioner::partition() {
    // Implement partitioning logic using ParMETIS
}

std::vector<int> ParMetisPartitioner::getNeighbors() const {
    // Return neighbor partitions
    std::vector<int> u;
    return u;
}

int ParMetisPartitioner::getColor() const {
    // Return the color of the partition
    return -1;
}
