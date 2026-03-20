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
 * @file ArchiveInputDeck.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_ARCHIVEINPUTDECK_HPP_
#define GEOS_FILEIO_OUTPUTS_ARCHIVEINPUTDECK_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace archiveInputDeck
{

/**
 * @brief Copy XML input files into the output directory, preserving the
 *        folder structure
 * @param inputFileNames Container of XML file names to start the copy from
 * @param outputDirectory The output directory to copy files into
 *
 * Copy XML input files and every included files they contain (specified in
 * the <Included> tag. This function creates a somewhat similar folder
 * structure to the actual structure in the disk.
 *
 * Note: XML files that are located "behind" the callpoint (the path to
 *       the first input file given as the -i paramater) will be prefixed
 *       with "__" for every "../" in the relative path from the callpoint.
 */
void archiveInputDeck( string_array const & inputFileNames,
                       string const & outputDirectory );

} /* namespace archiveInputDeck */

} /* namespace geos */


#endif // GEOS_FILEIO_OUTPUTS_ARCHIVEINPUTDECK_HPP_
