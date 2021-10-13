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
 * @file VTKOutput.hpp
 */

#ifndef GEOSX_FILEIO_OUTPUTS_VTKOUTPUT_HPP_
#define GEOSX_FILEIO_OUTPUTS_VTKOUTPUT_HPP_

#include "OutputBase.hpp"
#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"


namespace geosx
{

/**
 * @brief A class for creating vtk outputs
 */
class VTKOutput : public OutputBase
{
public:

  /// @copydoc geosx::dataRepository::Group::Group(string const & name, Group * const parent)
  VTKOutput( string const & name, Group * const parent );

  /// Destructor
  virtual ~VTKOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "VTK"; }

  virtual void postProcessInput() override;

  /**
   * @brief Writes out a set of vtk files.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Write one final set of vtk files as the code exits
   * @copydoc ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain );
  }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct : OutputBase::viewKeysStruct
  {
    static constexpr auto plotFileRoot = "plotFileRoot";
    static constexpr auto writeFEMFaces = "writeFEMFaces";
    static constexpr auto plotLevel = "plotLevel";
    static constexpr auto binaryString = "writeBinaryData";

  } vtkOutputViewKeys;
  /// @endcond

private:

  string m_plotFileRoot;
  integer m_writeFaceMesh;
  integer m_plotLevel;
  integer m_writeBinaryData;

  vtk::VTKPolyDataWriterInterface m_writer;

};


} /* namespace geosx */

#endif /* GEOSX_FILEIO_OUTPUTS_VTKOUTPUT_HPP_ */
