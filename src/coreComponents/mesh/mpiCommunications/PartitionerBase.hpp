
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
 * @file PartitionerBase.hpp
 */

#ifndef GEOS_PARTITIONER_PARTITIONERBASE_HPP_
#define GEOS_PARTITIONER_PARTITIONERBASE_HPP_

#include<vector>


class PartitionerBase {
public:
    virtual void partition() = 0;

    virtual std::vector<int> getNeighbors() const = 0;
    virtual std::vector<int> getNeighbors() = 0;

    virtual int getColor() const = 0;
    virtual int getColor() = 0;

    virtual ~PartitionerBase() = default;

    PartitionerBase() : m_numColors(-1), m_size(-1) {}


private:
    int m_numColors;
    int m_size;
};

#endif // GEOS_PARTITIONER_PARTITIONERBASE_HPP_
