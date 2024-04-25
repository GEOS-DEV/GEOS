/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StencilOutput.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_STENCILOUTPUT_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_STENCILOUTPUT_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geos
{

/**
 * @class StencilOutput
 * Task class allowing the ouput of communicating elements data.
 * Only current usage is to output the cell-cell transmissibility (= "CellToCellOutput"), but we could add more
 * implementations to output the transmissibility of boundary or surface connections (or everything that is managed
 * by the stencils).
 */
class StencilOutput : public FieldStatisticsBase< FlowSolverBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  StencilOutput( const string & name,
                 Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "CellToCellOutput"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

private:

  using Base = FieldStatisticsBase< FlowSolverBase >;
  using LocalToGlobalMap = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  class StencilOutputKernel;

  /**
   * @brief Element-element connection data extracted from the kernel.
   */
  struct KernelConnectionData
  {
    real64 transmissibility[2];
    int regionId[2];
    int subRegionId[2];
    int elementId[2];
  };

  /**
   * @brief Element-element connection data to be output by the class.
   */
  struct ConnectionData
  {
    real64 transmissibility[2];
    int globalId[2];

    /**
     * @brief minus operator to be able to sort instances before outputing them.
     * The sorting is done by the first global id, then by the second global id.
     * @param other compared instance
     * @return true if this instance is to sort before the other instance.
     */
    bool operator<( ConnectionData const & other ) const
    {
      return globalId[0] != other.globalId[0] ? globalId[0] < other.globalId[0] :
             globalId[1] < other.globalId[1];
    }
  };


  /**
   * @brief the list of the current registered output target filepaths.
   */
  std::set< string > outputFiles;


  /**
   * @brief gather the element-element connection data of the current timestep using a given StencilWrapper.
   * @tparam STENCILWRAPPER_T the type of the StencilWrapper
   * @param mesh the mesh for which we want the data
   * @param stencilWrapper the StencilWrapper to use to compute the element-element connection data
   * @return Return the gathered data in an LvArray
   */
  template< typename STENCILWRAPPER_T >
  array1d< KernelConnectionData > gatherTimeStepData( MeshLevel const & mesh,
                                                      STENCILWRAPPER_T const & stencilWrapper ) const;

  /**
   * @return a converted ConnectionData from the given kernel data.
   * @param kernelData the kernel extracted data to be converted.
   * @param localToGlobalMap the local to global id map, useful to convert the kernel extracted ids of the elements.
   */
  static ConnectionData ToConnectionData( KernelConnectionData const & kernelData,
                                          LocalToGlobalMap const & localToGlobalMap );

  /**
   * @brief Output the element-element connection data of the current timestep.
   * @param mesh the specific mesh for which we output the data. We will also need it to convert the ids to global ids.
   * @param stencil the specific mesh for which we output the data.
   * @param outputTime the time for when we gathered the data
   * @param kernelData the connection data, gathered by a kernel.
   */
  void outputTimeStepData( MeshLevel const & mesh, string_view stencilName, real64 outputTime,
                           arrayView1d< KernelConnectionData > const & kernelData );

  /**
   * @return the output file path for the given stencil, under the m_outputDir folder.
   * @param meshName the name of the mesh for which we want to output the data.
   * @param stencilName the name of the stencil for which we want to output the data.
   * @param extension the file extension
   */
  string getOutputFileName( string_view stencilName, string_view extension ) const;

  /**
   * @brief Write the CSV header.
   * @param outputFile the filestream where the header will be writen.
   */
  void writeCsvHeader( std::ofstream & outputFile ) const;

  /**
   * @brief Write the data in a CSV file.
   * @tparam CONN_DATA_CONT Type of the iterable container containing ConnectionData instances.
   * @param outputFile the filestream where the data will be writen.
   * @param data an iterable container of ConnectionData instances
   * @param time data timestep.
   */
  template< typename CONN_DATA_CONT >
  void writeCsvData( std::ofstream & outputFile, CONN_DATA_CONT const & data, real64 const time ) const;
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_ */
