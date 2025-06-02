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
 * @brief class to format data in a formatted text format, allowing contributions from multiple
 *        MPI ranks.
 */
class TableTextMpiOutput : public TableTextFormatter
{
public:

  using Base = TableTextFormatter;

  enum class ParallelOutputMode
  {
    /// Each ranks output cannot be mixed with the content of other ranks. Each rank flushes after full output.
    InsecableRanks,
    /// The rows of every ranks are output randomly. Each rank flushes after each row formatting.
    MixedRanksRows,
  };

  /**
   * @brief Construct a default Table Formatter without layout specification (to only insert data in it,
   * without any column / title). Feature is not tested.
   */
  TableTextMpiOutput( ParallelOutputMode parallelOutputMode = ParallelOutputMode::MixedRanksRows );

  /**
   * @brief Construct a new TableTextMpiOutput from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableTextMpiOutput( TableLayout const & tableLayout,
                      ParallelOutputMode parallelOutputMode = ParallelOutputMode::MixedRanksRows );

  /**
   * @brief Convert a data source to a table string.
   * @param tableData The data source to convert.
   * @param outputStream The same target output stream for all ranks, to output the table string
   *                     representation of the TableData. The output is partial, each rank
   *                     contributing to common output stream with their local data. It may be the
   *                     log or a file stream.
   */
  template< typename DATASOURCE >
  void toStream( std::ostream & outputStream, DATASOURCE const & tableData ) const;

private:

  // hiding toString() methods as they are not implemented with MPI support.
  using Base::toString;

  ParallelOutputMode m_parallelOutputMode;

  /**
   * @brief Expend the columns width to accomodate with the content of all MPI ranks.
   *        As it is based on MPI communications, every ranks must call this method.
   * @param columnsWidth The array to store the resulting columns width in.
   * @param tableGrid The grid of cells containing content.
   */
  void stretchColumnsByRanks( stdVector< size_t > & columnsWidth ) const;

};

}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLEMPICOMPONENTS_HPP */
