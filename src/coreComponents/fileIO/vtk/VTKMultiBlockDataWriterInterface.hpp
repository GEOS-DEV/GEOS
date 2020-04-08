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

#ifndef GEOSX_FILEIO_VTK_VTKMULTIBLOCKWRITERINTERFACE_HPP_
#define GEOSX_FILEIO_VTK_VTKMULTIBLOCKWRITERINTERFACE_HPP_

#include "common/DataTypes.hpp"

#include "managers/DomainPartition.hpp"

#include "mesh/CellElementSubRegion.hpp"

#include "dataRepository/WrapperBase.hpp"
#include "dataRepository/Wrapper.hpp"

#include "VTKPVDWriter.hpp"
#include "VTKVTMWriter.hpp"
#include "VTKGEOSXData.hpp"

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>

namespace geosx
{
using namespace dataRepository;
namespace vtk
{
/*!
 * @brief Encapsulate IO methods from vtk
 */
class VTKPolyDataWriterInterface
{
  public:
  VTKPolyDataWriterInterface( string const & outputName );


  /*!
   * @brief Sets the plot level
   * @details All fields have an associated plot level. If it is <= to \p plotLevel,
   * the field will be output.
   * @param[in] plotLevel the limit plotlevel
   */
  void SetPlotLevel(integer plotLevel )
  {
    m_plotLevel = dataRepository::IntToPlotLevel(plotLevel);
  }

  /*!
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
  void Write( real64 time, integer cycle, DomainPartition const & domain  );

  private:

  /*!
   * @brief Create a folder at the given time-step \p time
   * @details the name of the folder will be the time-step. This folder
   * will contains every files concerning the time-step \p time.
   * (aka one file per ElementRegion and per rank).
   * @param[in] time the time-step
   */
  void CreateTimeStepSubFolder( real64 time ) const;

  /*!
   * @brief Given a time-step \p time, returns the relative path
   * to the subfolder containing the files concerning this time-step
   * @param[in] time the time-step
   * @return the relative path to the folder of the time step
   */
  string GetTimeStepSubFolder( real64 time ) const;

  /*!
   * @brief Writes the files for all the CellElementRegions.
   * @details There will be one file written per CellElementRegion and per rank.
   * @param[in] time the time-step
   * @param[in] domain the computation domain for this rank
   */
  void WriteCellElementRegions( real64 time, DomainPartition const & domain ) const;

  /*!
   * @brief Gets the VTK Object points encapsulating
   * the cells connectivities of \p er
   * @param[in] er the CellElementRegion to be written
   * @return a tuple, first value is a VTK object containing the connectivity information,
   * the second value is a table with the same size than the total number of element in the CellElementRegion
   * containg the type of the cells.
   */
  std::tuple< vtkSmartPointer< vtkCellArray >, std::vector< int> > GetVTKCells( CellElementRegion const & er ) const;

  /*!
   * @brief Gets the VTK Object points encapsulating
   * the vertices coordinates of \p nodeManager
   * @param[in] nodeManager the NodeManager associated with the domain being written
   */
  vtkSmartPointer< vtkPoints > GetVTKPoints( NodeManager const & nodeManager ) const;

  /*!
   * @brief Writes the files containing the well representation
   * @details There will be one file written per WellElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] domain the computation domain for this rank
   */
  void WriteWellElementRegions( real64 time, DomainPartition const & domain ) const;

  /*!
   * @brief Gets the VTK Object points encapsulating
   * the cells connectivities of \p es
   */
  std::tuple< vtkSmartPointer< vtkPoints >,  vtkSmartPointer< vtkCellArray > >GetWell( WellElementSubRegion  const & esr , NodeManager const & nodeManager) const;

  /*!
   * @brief Writes the files containing the faces elements
   * @details There will be one file written per FaceElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] domain the computation domain for this rank
   */
  void WriteFaceElementRegions( real64 time, DomainPartition const & domain ) const;

  std::tuple< vtkSmartPointer< vtkPoints >,  vtkSmartPointer< vtkCellArray > >GetSurface( FaceElementSubRegion  const & esr , NodeManager const & nodeManager) const;

  /*!
   * @brief Writes a VTM file for the time-step \p time.
   * @details a VTM file is a VTK Multiblock file. It contains reltive path to different files organized in blocks.
   * @param[in] time the time-step
   * @param[in] domain the computation domain for this rank
   * @param[in] vtmWrite a writer specialized for the VTM file format
   */
  void WriteVTMFile( real64 time, DomainPartition const & domain, VTKVTMWriter const & vtmWriter ) const;

  /*!
   * @brief Write all the fields associated to the nodes of \p nodeManager if their plotlevel is <= m_plotLevel
   * @param[in] pointdata a VTK object containing all the fields associated with the nodes
   * @param[in] nodeManager the NodeManager associated with the domain being written
   */
  void WriteNodeFields( vtkSmartPointer< vtkPointData > const pointdata, NodeManager const & nodeManager) const;

  /*!
   * @brief Writes all the fields associated to the cellss of \p er if their plotlevel is <= m_plotLevel
   * @param[in] celldata a VTK object containing all the fields associated with the cells
   * @param[in] er CellElementRegion being written
   */
  void WriteCellFields( vtkSmartPointer< vtkCellData > const celldata, CellElementRegion const & er ) const;

  /*!
   * @brief Writes a field from \p wrapperBase
   * @details Sets the number of components, the number of value and fill the VTK data structure using
   * a wrapper around a field.
   * @param[in] wrapperBase a wrapper around the field to be written
   * @param[in,out] data a VTK data container derived to be suitable for some GEOSX types.
   * @param[in] size the number of values in the field
   * @param[in,out] a counter that is incremented each time a value is written. This is useful
   * for CellElementSubRegion.
   */
  void WriteField( WrapperBase const & wrapperBase, vtkSmartPointer < VTKGEOSXData > data, localIndex size, localIndex & count ) const;

  private:

  /// Folder name in which all the files will be written
  string const m_outputFolder;
  
  /// A writter specialized for PVD files. There is one PVD file per simulation. It is the root
  /// file containing all the paths to the VTM files.
  VTKPVDWriter m_pvd;
  
  /// Maximum plot level to be written.
  dataRepository::PlotLevel m_plotLevel;

  /// The previousCycle
  integer m_previousCycle;
};

} // namespace vtk
} // namespace geosx

#endif
