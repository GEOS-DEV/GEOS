/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_HDFFILE_HPP
#define GEOSX_HDFFILE_HPP

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @class HDFFile
 * A class used to control access to an HDF file target.
 */
class HDFFile
{
public:
  /**
   * @brief Constructor -- this creates/opens the target file for read/write.
   * @param fnm The filename.
   * @param deleteExisting Whether to remove/recreate if a file with the same name exists.
   * @param parallelAccess Whether to access one file in parallel or one file per rank in the comm.
   * @param comm An MPI communicator where each rank in the communicator will be accesing the target file.
   */
  HDFFile( string const & fnm, bool deleteExisting, bool parallelAccess, struct ompi_communicator_t * comm );

  /**
   * Destructor -- Close the file and acccessors.
   */
  ~HDFFile();

  /**
 * @brief Whether a dataset/group with the specified name exists in the target.
 * @param name The dataset/group name to check for.
 * @return Whether the dataset/group exists in the target.
 */
  bool checkInTarget( const string & name ) const;

  /**
   * @brief Get the HDF hid_t file identifier.
   * @return the HDF hid_t file id.
   */
  operator int64_t() const { return m_fileId; }
private:
  /// The filename
  string m_filename;
  /// The hdf file id
  int64_t m_fileId;
  /// The hdf file access properties list id
  int64_t m_faplId;
  /// Whether to open the same file in parallel, or open a file per process
  bool m_mpioFapl;
  /// The comminator to operate on the file collectively over
  struct ompi_communicator_t * m_comm;
};

}

#endif //GEOSX_HDFFILE_HPP
