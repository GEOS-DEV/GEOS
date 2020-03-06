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

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>

namespace geosx
{
namespace vtk
{
  class VTKPolyDataWriterInterface
  {
    public:
    VTKPolyDataWriterInterface( string const & outputName );


    /*!
     * @brief Writes the .pvd file
     */
    void Write( real64 time, DomainPartition const * const domain  ) const;

    private:
    /*!
     * @brief Gets the VTK Object points encapsulating
     * the vertices coordinates of \p nodeManager
     */
    vtkSmartPointer< vtkPoints > GetVTKPoints( NodeManager const * const nodeManager ) const;

    /*!
     * @brief Gets the VTK Object points encapsulating
     * the cells connectivities of \p es
     */
    vtkSmartPointer< vtkCellArray > GetVTKCells( CellElementSubRegion const * const esr ) const;

    /*!
     * @brief Sets the name of the file to be output.
     */
    void SetFileName( double time ) const;

    void LinkMesh( DomainPartition const * const domain ) const;

    private:
    vtkSmartPointer< vtkXMLUnstructuredGridWriter > m_writer;

    string const m_outputFolder;
  };

} // namespace vtk
} // namespace geosx

#endif
