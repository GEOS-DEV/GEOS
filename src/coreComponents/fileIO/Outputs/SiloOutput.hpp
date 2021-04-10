/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SiloOutput.hpp
 */

#ifndef GEOSX_FILEIO_OUTPUTS_SILOOUTPUT_HPP_
#define GEOSX_FILEIO_OUTPUTS_SILOOUTPUT_HPP_

#include "OutputBase.hpp"


namespace geosx
{

/**
 * @class SiloOutput
 *
 * A class for creating silo-based outputs
 */
class SiloOutput : public OutputBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(string const & name, Group * const parent)
  SiloOutput( string const & name,
              Group * const parent );

  /// Destructor
  virtual ~SiloOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "Silo"; }

  /**
   * @brief Writes out a Silo plot file.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Writes out a Silo plot file at the end of the simulation.
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
    static constexpr auto writeEdgeMesh = "writeEdgeMesh";
    static constexpr auto writeFaceMesh = "writeFEMFaces";
    static constexpr auto writeCellElementMesh = "writeCellElementMesh";
    static constexpr auto writeFaceElementMesh = "writeFaceElementMesh";
    static constexpr auto plotLevel = "plotLevel";

  } siloOutputViewKeys;
  /// @endcond

private:
  string m_plotFileRoot;
  integer m_writeEdgeMesh;
  integer m_writeFaceMesh;
  integer m_writeCellElementMesh;
  integer m_writeFaceElementMesh;
  integer m_plotLevel;

};


} /* namespace geosx */

#endif /* GEOSX_FILEIO_OUTPUTS_SILOOUTPUT_HPP_ */
