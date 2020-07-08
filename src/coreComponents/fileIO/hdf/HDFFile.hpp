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
 * @file HDFFile.hpp
 */

#ifndef GEOSX_HDFFILE_HPP_
#define GEOSX_HDFFILE_HPP_

#include "LvArray/src/Array.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

#include "managers/TimeHistory/HistoryDataSpec.hpp"
#include "managers/TimeHistory/HistoryIO.hpp"
#include <hdf5.h>

#include <string>

namespace geosx
{

/**
 * @class HDFTarget
 * @brief An abstract class representing an HDF output target.
 */
class HDFTarget
{
public:
  /**
   * @brief Get the HDF hid_t of the target.
   * @return The hid_t of the target.
   */
  virtual operator hid_t() { return 0; }

  /**
   * @brief Whether a dataset/group with the specified name exists in the target.
   * @param name The dataset/group name to check for.
   * @return Whether the dataset/group exists in the target.
   */
  virtual bool CheckInTarget( const string & name )
  {
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Gget_objinfo(this->operator hid_t(), name.c_str(), 0, NULL);
    } H5E_END_TRY
    return (exists == 0);
  }
};

/**
 * @class HDFile
 * @brief A class used to control access to an HDF file target.
 */
class HDFFile : public HDFTarget
{
public:
  /**
   * @brief Constructor -- this creates/opens the target file for read/write.
   * @param fnm The filename.
   * @param delete_existing Whether to remove/recreate if a file with the same name exists.
   * @param comm An MPI communicator where each rank in the communicator will be accesing the target file.
   */
  HDFFile(string const & fnm, bool delete_existing = false, MPI_Comm comm = MPI_COMM_GEOSX);

  /**
   * Destructor -- Close the file and acccessors.
   */
  ~HDFFile();

  /**
   * @brief Get the HDF hid_t file identifier.
   * @return the HDF hid_t file id.
   */
  virtual operator hid_t() final { return m_file_id; }
private:
  string m_filename;
  hid_t m_file_id;
  hid_t m_fapl_id;
  hid_t m_dxpl_id;
  MPI_Comm m_comm;
};

/**
 * @class HDFHistIO
 * @brief Perform buffered history I/O for a single type(really just output) on using HDF5.
 */
class HDFHistIO : public BufferedHistoryIO
{
  public:
  /**
   * @brief Constructor
   * @param filename The filename to perform history output to.
   * @paran rank The rank of the history data being collected.
   * @param dims The dimensional extent for each dimension of the history data being collected.
   * @param name The name to use to create/modify the dataset for the history data.
   * @param type_id The std::type_index(typeid(T)) of the underlying data type.
   * @param write_head How many time history states have been written to the file (used on restart and to compress data on exit).
   * @param init_alloc How many states to preallocate the internal buffer to hold.
   * @param comm A communicator where every rank will participate in writting to the output file.
   */
  HDFHistIO( string const & filename,
             localIndex rank,
             const localIndex * dims,
             string const & name,
             std::type_index type_id,
             localIndex write_head = 0,
             localIndex init_alloc = 4,
             MPI_Comm comm = MPI_COMM_GEOSX);

  /**
   * @brief Constructor
   * @param filename The filename to perform history output to.
   * @param spec HistoryMetadata to use to call the other constructor.
   * @param write_head How many time states have been written to the file (used on restart and to compress data on exit).
   * @param init_alloc How many states to preallocate the internal buffer to hold.
   * @param comm A communicator where every rank will participate in writing to the output file.
   */
  HDFHistIO( string const & filename, const HistoryMetadata & spec, localIndex write_head = 0, localIndex init_alloc = 4, MPI_Comm comm = MPI_COMM_GEOSX ) :
    HDFHistIO( filename, spec.getRank(), spec.getDims(), spec.getName(), spec.getType(), write_head, init_alloc, comm )
  { }

  /// Destructor
  virtual ~HDFHistIO() { }

  /// @copydoc BufferedHistoryIO::Init( )
  virtual void Init( bool exists_okay ) override;

  /// @copydoc HistoryMetadataIO::Write( )
  virtual void Write( ) override;

  /// @copydoc BufferedHistoryIO::CompressInFile( )
  virtual void CompressInFile( ) override;

  /**
   * @brief Resize the dataspace in the target file if needed to perform the current write of buffered states.
   * @param buffered_count The number of buffered states to use to determine if the file needs to be resized.
   */
  inline void ResizeFileIfNeeded( localIndex buffered_count );

  /// @copydoc BufferedHistoryIO::GetRankOffset( )
  virtual globalIndex GetRankOffset( ) override;

  protected:
  virtual void resizeBuffer( ) override;

  private:
    // file io params
    string m_filename;
    const localIndex m_overalloc_multiple;
    globalIndex m_global_idx_offset;
    globalIndex m_global_idx_count;
    localIndex m_write_limit;
    localIndex m_write_head;
    // history metadata
    hsize_t m_hdf_type;
    size_t m_type_size;
    hsize_t m_type_count; // prod(dims[0:n])
    hsize_t m_rank;
    std::vector< hsize_t > m_dims;
    string m_name;
    MPI_Comm m_comm;
    MPI_Comm m_subcomm;
};

}

#endif
