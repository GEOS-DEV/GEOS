/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_FILEIO_TIMEHISTORY_HDFFILE_HPP_
#define GEOS_FILEIO_TIMEHISTORY_HDFFILE_HPP_


#include "dataRepository/HistoryDataSpec.hpp"
#include "BufferedHistoryIO.hpp"
#include "common/DataTypes.hpp"

#include <hdf5.h>

namespace geos
{

/**
 * @class HDFHistoryIO
 * @brief Perform buffered history I/O for a single type(really just output) on using HDF5.
 */
class HDFHistoryIO : public BufferedHistoryIO
{
public:
  /**
   * @brief Constructor
   * @param[in] filename The filename to perform history output to.
   * @param[in] rank The rank of the history data being collected.
   * @param[in] dims The dimensional extent for each dimension of the history data being collected.
   * @param[in] name The name to use to create/modify the dataset for the history data.
   * @param[in] typeId The std::type_index(typeid(T)) of the underlying data type.
   * @param[in] writeHead How many time history states have been written to the file (used on restart and to compress data on exit).
   * @param[in] initAlloc How many states to preallocate the internal buffer to hold.
   * @param[in] overallocMultiple Integer to scale the internal buffer when we fill the existing space.
   * @param[in] comm A communicator where every rank will participate in writting to the output file.
   */
  HDFHistoryIO( string const & filename,
                localIndex rank,
                std::vector< localIndex > const & dims,
                string const & name,
                std::type_index typeId,
                localIndex writeHead = 0,
                MPI_Comm comm = MPI_COMM_GEOSX );

  /**
   * @brief Constructor
   * @param[in] filename The filename to perform history output to.
   * @param[in] spec HistoryMetadata to use to call the other constructor.
   * @param[in] writeHead How many time states have been written to the file (used on restart and to compress data on exit).
   * @param[in] initAlloc How many states to preallocate the internal buffer to hold.
   * @param[in] overallocMultiple Integer to scale the internal buffer when we fill the existing space.
   * @param[in] comm A communicator where every rank will participate in writing to the output file.
   */
  HDFHistoryIO( string const & filename,
                const HistoryMetadata & spec,
                localIndex writeHead = 0,
                MPI_Comm comm = MPI_COMM_GEOSX ):
    HDFHistoryIO( filename,
                  spec.getRank(),
                  spec.getDims(),
                  spec.getName(),
                  spec.getType(),
                  writeHead,
                  comm )
  { }

  /// Destructor
  virtual ~HDFHistoryIO() { }

  /// @copydoc geos::BufferedHistoryIO::init
  virtual void init( bool existsOkay ) override;

  /// @copydoc geos::BufferedHistoryIO::write
  virtual void write( ) override;

  /// @copydoc geos::BufferedHistoryIO::compressInFile
  virtual void finalize( ) override;

  /// @copydoc geos::BufferedHistoryIO::shareable
  virtual bool shareable( ) const override { return true; }


private:

  /**
   * @brief Setup the parallel 'partitioning' of the data to allow dynamically sized output over time
   * @param[in] localIdxCount The number of pieces of data associated with the local rank
   * @note This is collective over the communicator provided in the constructor
   */
  void setupPartition( globalIndex localIdxCount );

  /**
   * @brief Update the extent of the dataset in the target file.
   * @param[in] rowLimit The new discrete 'row' count (the first output dim) to extend the extent to.
   * @note The second dimension is set to the global index highwater ( the largest number of output
   *       pieces of data encountered during execution ).
   */
  void updateDatasetExtent( hsize_t rowLimit );

  /**
   * @brief Resize the dataspace in the target file if needed to perform the current write of buffered states.
   * @param[in] bufferedCount The number of buffered states to use to determine if the file needs to be resized.
   */
  void resizeFileIfNeeded( localIndex bufferedCount );

  // file io params
  /// The filename to write to
  string m_filename;
  /// The global index offset for this mpi rank for this data set
  globalIndex m_globalIdxOffset;
  /// The global index count for this mpi rank for this data set
  globalIndex m_globalIdxCount;
  /// The largest index count encountered (globally) during execution (this determines the 2nd dimension of
  ///   the output data array)
  globalIndex m_globalIdxHighwater;
  /// The current chunk size for writing to file. This is the smallest (necessarily nonzero) local index
  ///   count from the last partition setup.
  hsize_t m_chunkSize;
  /// The current limit in discrete history counts for this data set in the file
  localIndex m_writeLimit;
  /// The current history count for this data set in the file
  localIndex m_writeHead;

  // history metadata
  /// The underlying data type for this history data set
  hsize_t m_hdfType;

  /// The name of the data set
  string m_name;
  /// The communicator across which the data set is distributed
  MPI_Comm m_comm;
  /// The communicator with only members of the m_comm comm which have nonzero ammounts of local data (required for chunking output ->
  /// growing the data size in the file)
  MPI_Comm m_subcomm;
};

}

#endif
