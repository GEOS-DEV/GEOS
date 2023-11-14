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

#ifndef GEOS_FILEIO_TIMEHISTORY_CSVHISTORYIO_HPP_
#define GEOS_FILEIO_TIMEHISTORY_CSVHISTORYIO_HPP_


#include "dataRepository/HistoryDataSpec.hpp"
#include "BufferedHistoryIO.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @class CSVHistoryIO
 * @brief Perform buffered history I/O for a single type(really just output) to CSV files.
 */
class CSVHistoryIO : public BufferedHistoryIO
{
public:
  /**
   * @brief Constructor
   * @param[in] filename The filename to perform history output to.
   * @param[in] rank The rank of the history data being collected.
   * @param[in] dims The dimensional extent for each dimension of the history data being collected.
   * @param[in] typeId The std::type_index(typeid(T)) of the underlying data type.
   * @param[in] writeHead How many time history states have been written to the file (used on restart and to compress data on exit).
   * @param[in] initAlloc How many states to preallocate the internal buffer to hold.
   * @param[in] overallocMultiple Integer to scale the internal buffer when we fill the existing space.
   * @param[in] comm A communicator where every rank will participate in writting to the output file.
   */
  CSVHistoryIO( string const & filename,
                localIndex rank,
                std::vector< localIndex > const & dims,
                std::type_index typeId,
                localIndex writeHead = 0,
                MPI_Comm comm = MPI_COMM_GEOSX );

  /**
   * @brief Constructor
   * @param[in] filename The filename to perform history output to.
   * @param[in] spec HistoryMetadata to use to call the other constructor.
   * @param[in] writeHead How many time states have been written to the file (used on restart and to compress data on exit).
   * @param[in] comm A communicator where every rank will participate in writing to the output file.
   */
  CSVHistoryIO( string const & filename,
                const HistoryMetadata & spec,
                localIndex writeHead = 0,
                MPI_Comm comm = MPI_COMM_GEOSX ):
    CSVHistoryIO( filename,
                  spec.getRank(),
                  spec.getDims(),
                  spec.getType(),
                  writeHead,
                  comm )
  { }

  /// Destructor
  virtual ~CSVHistoryIO() { }

  /// @copydoc geos::BufferedHistoryIO::init
  virtual void init( bool existsOkay ) override;

  /// @copydoc geos::BufferedHistoryIO::write
  virtual void write( ) override;

  /// @copydoc geos::BufferedHistoryIO::finalize
  virtual void finalize( ) override {};

  /// @copydoc geos::BufferedHistoryIO::shareable
  virtual bool shareable( ) const override { return false; }

private:
  // file io params
  /// The filename to write to
  string m_filename;
  /// The current history count for this data set in the output file
  localIndex m_writeHead;

  /// The communicator across which the data set is distributed
  MPI_Comm m_comm;
};

}

#endif
