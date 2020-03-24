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

/**
 * @file SiloOutput.hpp
 */

#ifndef GEOSX_MANAGERS_OUTPUTS_SILOOUTPUT_HPP_
#define GEOSX_MANAGERS_OUTPUTS_SILOOUTPUT_HPP_

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
  /// Main constructor
  SiloOutput( std::string const & name,
              Group * const parent );

  /// Destructor
  virtual ~SiloOutput() override;

  /// Catalog name interface
  static string CatalogName() { return "Silo"; }

  /// This method will be called by the event manager if triggered
  virtual void Execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override;

  /// Write one final output as the code exits
  virtual void Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        dataRepository::Group * domain ) override
  {
    Execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain );
  }

  struct viewKeysStruct : OutputBase::viewKeysStruct
  {
    static constexpr auto plotFileRoot = "plotFileRoot";
    static constexpr auto writeEdgeMesh = "writeEdgeMesh";
    static constexpr auto writeFaceMesh = "writeFEMFaces";
    static constexpr auto writeCellElementMesh = "writeCellElementMesh";
    static constexpr auto writeFaceElementMesh = "writeFaceElementMesh";
    static constexpr auto plotLevel = "plotLevel";

  } siloOutputViewKeys;

private:
  string m_plotFileRoot;
  integer m_writeEdgeMesh;
  integer m_writeFaceMesh;
  integer m_writeCellElementMesh;
  integer m_writeFaceElementMesh;
  integer m_plotLevel;

};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_SILOOUTPUT_HPP_ */
