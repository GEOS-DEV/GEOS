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
 * @file CommonMeshOutputUtilities.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_COMMONMESHOUTPUTUTILITIES_HPP_
#define GEOS_FILEIO_OUTPUTS_COMMONMESHOUTPUTUTILITIES_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/RestartFlags.hpp"

namespace conduit
{
class Node;
}

namespace geos
{

class MeshLevel;

namespace dataRepository
{
class Group;
}

namespace outputUtilities
{

/**
 * @brief Common options for building a Conduit mesh tree for post-processing metadata.
 */
struct CommonMeshOptions
{
  dataRepository::PlotLevel plotLevel = dataRepository::PlotLevel::LEVEL_1;
  int outputFullQuadratureData = 0;
  int writeGhostObjects = 1;
  int writeFieldData = 1;
  int writeGhostFieldData = 1;
};

/**
 * @brief Common options for writing the shared HDF5 heavy data.
 */
struct CommonHDFWriteOptions
{
  string outputSubdirectory;
  int writeParallelFiles = 1;
};

/**
 * @brief Output details returned by the shared HDF5 writer.
 */
struct CommonHDFWriteResult
{
  string cyclePath;
  string filePathForRank;
  string hdfFileName;
  integer numberOfFiles = 0;
  integer fileRank = 0;
  int wroteFileOnThisRank = 0;
};

/**
 * @brief Build a Conduit mesh tree with common GEOS mesh/field data.
 * @param meshLevel The mesh level to serialize.
 * @param time_n Current simulation time.
 * @param cycleNumber Current cycle number.
 * @param options Shared mesh output options.
 * @param parentGroup Group used for temporary averaged constitutive fields.
 * @param meshRoot Root node to populate. Expected layout: /mesh/...
 */
void buildCommonMeshTree( MeshLevel const & meshLevel,
                          real64 const time_n,
                          integer const cycleNumber,
                          CommonMeshOptions const & options,
                          dataRepository::Group & parentGroup,
                          conduit::Node & meshRoot );

/**
 * @brief Write shared HDF5 heavy data and an associated root file.
 * @param meshRoot Root Conduit node to write.
 * @param outputDirectory Base output directory.
 * @param cycleNumber Current cycle number.
 * @param options Shared HDF write options.
 * @param rootNode Optional pre-populated root node to write (.root file payload).
 * @return Information about files written and file naming.
 */
CommonHDFWriteResult writeCommonMeshHDF5( conduit::Node & meshRoot,
                                          string const & outputDirectory,
                                          integer const cycleNumber,
                                          CommonHDFWriteOptions const & options,
                                          conduit::Node * rootNode = nullptr );

} /* namespace outputUtilities */

} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_COMMONMESHOUTPUTUTILITIES_HPP_ */
