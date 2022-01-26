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

#ifndef GEOSX_FILEIO_VTK_VTKPOLYDATAWRITERINTERFACE_HPP_
#define GEOSX_FILEIO_VTK_VTKPOLYDATAWRITERINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "dataRepository/Wrapper.hpp"
#include "fileIO/vtk/VTKPVDWriter.hpp"
#include "fileIO/vtk/VTKVTMWriter.hpp"

class vtkUnstructuredGrid;
class vtkPointData;
class vtkCellData;

namespace geosx
{

class DomainPartition;
class ElementRegionBase;
class EmbeddedSurfaceNodeManager;
class ElementRegionManager;
class NodeManager;

namespace vtk
{

enum struct VTKOutputMode
{
  BINARY,
  ASCII
};

/**
 * @brief Encapsulate output methods for vtk
 */
class VTKPolyDataWriterInterface
{
public:
  /**
   * @brief Constructor
   * @param[in] outputName folder name in which all the files will be written
   */
  explicit VTKPolyDataWriterInterface( string outputName );

  /**
   * @brief Sets the plot level
   * @details All fields have an associated plot level. If it is <= to \p plotLevel,
   * the field will be output.
   * @param[in] plotLevel the limit plotlevel
   */
  void setPlotLevel( integer plotLevel )
  {
    m_plotLevel = dataRepository::toPlotLevel( plotLevel );
  }

  /**
   * @brief Set the binary mode
   * @param[in] mode output mode to be set
   */
  void setOutputMode( VTKOutputMode mode )
  {
    m_outputMode = mode;
  }

  /**
   * @brief Set the output directory name
   * @param[in] outputDir global output directory location
   * @param[in] outputName name of the VTK output subdirectory and corresponding PVD file
   */
  void setOutputLocation( string outputDir, string outputName )
  {
    m_outputDir = std::move( outputDir );
    m_outputName = std::move( outputName );
    m_pvd.setFileName( joinPath( m_outputDir, m_outputName ) + ".pvd" );
  }

  /**
   * @brief Main method of this class. Write all the files for one time step.
   * @details This method writes a .pvd file (if a previous one was created from a precedent time step,
   * it is overwritten). The .pvd file contains relative path to every .vtm files (one vtm file per time step).
   * This method triggers also the writing of a .vtm file. A .vtm file containts relative paths to blocks
   * with the following hierarchy :
   *  - CellElementRegion
   *    - CellElementRegion1
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - CellElementRegion2
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - ...
   *  -WellElementRegion
   *    - Well1
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - Well2
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   * @param[in] time the time step to be written
   * @param[in] cycle the current cycle of event
   * @param[in] domain the computation domain of this rank
   */
  void write( real64 time, integer cycle, DomainPartition const & domain );

private:

  /**
   * @brief Given a time-step \p time, returns the relative path
   * to the subfolder containing the files concerning this time-step
   * @param[in] cycle the current cycle number
   * @return the relative path to the folder of the time step
   */
  string getCycleSubFolder( integer const cycle ) const;

  /**
   * @brief Writes the files for all the CellElementRegions.
   * @details There will be one file written per CellElementRegion and per rank.
   * @param[in] time the time-step
   * @param[in] cycle the current cycle number
   * @param[in] elemManager the ElementRegionManager containing the CellElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void writeCellElementRegions( real64 time,
                                integer const cycle,
                                ElementRegionManager const & elemManager,
                                NodeManager const & nodeManager ) const;

  /**
   * @brief Writes the files containing the well representation
   * @details There will be one file written per WellElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] cycle the current cycle number
   * @param[in] elemManager the ElementRegionManager containing the WellElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void writeWellElementRegions( real64 time,
                                integer const cycle,
                                ElementRegionManager const & elemManager,
                                NodeManager const & nodeManager ) const;

  /**
   * @brief Writes the files containing the faces elements
   * @details There will be one file written per FaceElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] cycle the current cycle number
   * @param[in] elemManager the ElementRegionManager containing the FaceElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void writeSurfaceElementRegions( real64 time,
                                   integer const cycle,
                                   ElementRegionManager const & elemManager,
                                   NodeManager const & nodeManager,
                                   EmbeddedSurfaceNodeManager const & embSurfNodeManager ) const;

  /**
   * @brief Write all the fields associated to the nodes of \p nodeManager if their plotlevel is <= m_plotLevel
   * @param[in] pointData a VTK object containing all the fields associated with the nodes
   * @param[in] nodeManager the NodeManager associated with the domain being written
   */
  void writeNodeFields( NodeManager const & nodeManager,
                        vtkPointData & pointData ) const;

  /**
   * @brief Writes all the fields associated to the elements of \p er if their plotlevel is <= m_plotLevel
   * @param[in] subRegion ElementRegion being written
   * @param[in] cellData a VTK object containing all the fields associated with the elements
   */
  template< class SUBREGION >
  void writeElementFields( ElementRegionBase const & subRegion,
                           vtkCellData & cellData ) const;

  /**
   * @brief Writes an unstructured grid
   * @details The unstructured grid is the last element in the hierarchy of the output,
   * it contains the cells connectivities and the vertices coordinates as long as the
   * data fields associated with it
   * @param[in] ug a VTK SmartPointer to the VTK unstructured grid.
   * @param[in] cycle the current cycle number
   * @param[in] name the name of the ElementRegionBase to be written
   */
  void writeUnstructuredGrid( integer const cycle,
                              string const & name,
                              vtkUnstructuredGrid & ug ) const;

private:

  /// Output directory name
  string m_outputDir;

  /// Name of the PVD file and associated directory containing the outputs
  string m_outputName;

  /// A writter specialized for PVD files. There is one PVD file per simulation. It is the root
  /// file containing all the paths to the VTM files.
  VTKPVDWriter m_pvd;

  /// Maximum plot level to be written.
  dataRepository::PlotLevel m_plotLevel;

  /// The previousCycle
  integer m_previousCycle;

  /// Output mode, could be ASCII or BINARAY
  VTKOutputMode m_outputMode;
};

} // namespace vtk
} // namespace geosx

#endif
