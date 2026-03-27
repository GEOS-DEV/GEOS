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
 * @file TableMpiComponents.hpp
 * @brief contains variation of tables components for MPI communication.
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLEMPICOMPONENTS_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLEMPICOMPONENTS_HPP

#include "TableFormatter.hpp"

namespace geos
{

/**
 * @struct TableMpiLayout
 * @brief Layout information specific to MPI distributed tables, completing those in TableLayout.
 */
struct TableMpiLayout
{
  /// Enable a separating line between ranks in the table output.
  bool m_separatorBetweenRanks = false;
  /// Title for each rank's section visible in the separating line.
  string m_rankTitle;
};

/**
 * @brief class to format data in a formatted text format, allowing contributions from multiple
 *        MPI ranks.
 */
class TableTextMpiFormatter : public TableTextFormatter
{
public:
  /// base class
  using Base = TableTextFormatter;
  /// Callable comparison function object used for std::sort for a TableData
  using SortingFunc = std::function< bool (stdVector< TableData::CellData >, stdVector< TableData::CellData >) >;
  /**
   * @brief Construct a default Table Formatter without layout specification (to only insert data in it,
   * without any column / title). Feature is not tested.
   * @param mpiLayout MPI-specific layout information (default is having contiguous ranks data).
   */
  TableTextMpiFormatter( TableMpiLayout mpiLayout = TableMpiLayout() );

  /**
   * @brief Construct a new TableTextMpiOutput from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   * @param mpiLayout MPI-specific layout information (default is having contiguous ranks data).
   */
  TableTextMpiFormatter( TableLayout const & tableLayout,
                         TableMpiLayout mpiLayout = TableMpiLayout() );

  /**
   * @brief Convert a data source to a table string.
   * @param tableData The data source to convert.
   * @param outputStream The target output stream for rank 0, to output the table string representation
   *                     of the TableData. Each rank contributing to the common rank 0 output stream
   *                     with their local data. It may be the log or a file stream.
   * @note This method must be called by all MPI ranks.
   */
  template< typename DATASOURCE >
  void toStream( std::ostream & outputStream, DATASOURCE const & tableData ) const;

  /**
   * @brief Set the Sorting Func object
   * @param func The callable comparison function object
   */
  void setSortingFunc( SortingFunc && func )
  { m_sortingFunctor = std::make_unique< SortingFunc >( std::move( func )); }


private:

  // hiding toString() methods as they are not implemented with MPI support.
  using Base::toString;

  struct Status
  {
    bool const m_isMasterRank;
    bool const m_isContributing;
    bool m_hasContent;
    string m_sepLine;
  };

  TableMpiLayout m_mpiLayout;

  /// The custom comparison function object for std::sort
  std::unique_ptr< SortingFunc > m_sortingFunctor;

  /**
   * @brief Expend the columns width to accomodate with the content of all MPI ranks.
   *        As it is based on MPI communications, every ranks must call this method.
   * @param columnsWidth The array to store the resulting columns width in.
   * @param tableGrid The grid of cells containing content.
   */
  void stretchColumnsByRanks( stdVector< size_t > & columnsWidth,
                              Status const & status ) const;

  /**
   * @brief Parse a string row to a TablaData cells.
   * @param rowString The string row string to parse.
   * @return The parsed row as a vector of CellData.
   */
  stdVector< TableData::CellData > parseStringRow( string_view rowString ) const;

  /**
   * @brief Gather all the TableData to the rank 0.
   * @param localTableData The local TableData to send to rank 0;
   */
  TableData gatherTableDataRank0( TableData const & localTableData ) const;


  /**
   * @brief Gather data cell rows across all MPI ranks and output them to rank 0 in rank order.
   * @param tableOutput The output stream to display the resulting table
   * @param dataCellsLayout  The layout for the data cells
   * @param tableLayout The layout of the table
   * @param status The TableMpi status for the current rank
   */
  void gatherAndOutputTableDataInRankOrder( std::ostream & tableOutput,
                                            CellLayoutRows const & dataCellsLayout,
                                            PreparedTableLayout const & tableLayout,
                                            TableTextMpiFormatter::Status & status )const;
};

}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLEMPICOMPONENTS_HPP */
