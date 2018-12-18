/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/**
 * @file GeometricalAggregation.hpp
 */

#pragma once

#include "AggregationManager.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * The Aggregation Manager stands for the aggregation of a group of fine cells into one coarse cell
 * It is based on the SparsityPattern given by the DofManager
 */
class GeometricalAggregation : public AggregationManager {
  public:
    static string CatalogName() { return "GeometricalAggragation"; }

    /**
     * @brief Performs the aggregrations
     * @details Implements the geometrical aggregation : The aggregation is done without taking
     * into account any geological or flow bases consideration. It is done simply using METIS within the domain.
     * @param[in]     fineSparsityPattern the SparsityPattern corresponding to the fine mesh
     * @param[in,out] coarseSparsityPattern the SparsityPattern corresponding to the fine mesh
     * @param[in]     totalNumberOfAggregates the total number of coarse cells.
     */
    virtual void Aggregate(SparsityPattern const & fineSparsityPattern,
                           SparsityPattern & coarseSparsityPattern,
                           int totalNumberOfAggregates) final;

};
}
