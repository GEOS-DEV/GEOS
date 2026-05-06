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
 * @brief Archive the XML input deck (and optionally the XSD schema) into the
 *        output directory.
 * @param inputFileNames Container of XML file names to start the copy from
 * @param outputDirectory The output directory to copy files into
 * @param xmlTagOrder The order of the XML tags in the XML archive file
 * @param level Archiving strategy level:
 *              - 0: no archiving (returns immediately)
 *              - 1: XML inputs only (flattened into a single file)
 *              - 2: XML inputs + the XSD schema
 *
 * Copy XML input files and every included files they contain (specified in
 * the Included tag) into a single flat file. When @p level is at least 2, the
 * XSD schema is also copied next to the flattened input.
 */
void archiveInputDeck( string_array const & inputFileNames,
                       string const & outputDirectory,
                       string_array const & xmlTagOrder,
                       integer level );

} /* namespace archiveInputDeck */

} /* namespace geos */


#endif // GEOS_FILEIO_OUTPUTS_ARCHIVEINPUTDECK_HPP_
