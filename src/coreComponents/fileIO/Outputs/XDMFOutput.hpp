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
 * @file XDMFOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_XDMFOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_XDMFOUTPUT_HPP_

#include "fileIO/Outputs/OutputBase.hpp"

namespace conduit
{
class Node;
}

namespace geos
{

/**
 * @class XDMFOutput
 * @brief A class for creating XDMF outputs with heavy data stored in HDF5.
 */
class XDMFOutput : public OutputBase
{
public:

  XDMFOutput( string const & name,
              Group * const parent );

  virtual ~XDMFOutput() override
  {}

  static string catalogName() { return "XDMF"; }

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  { execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain ); }

private:

  static string buildXdmfDocument( conduit::Node const & mesh,
                                   string const & hdfFileName,
                                   real64 time );

  dataRepository::PlotLevel m_plotLevel = dataRepository::PlotLevel::LEVEL_1;

  int m_outputFullQuadratureData = 0;

  int m_writeParallelFiles = 1;

  int m_writeGhostObjects = 1;

  int m_writeFieldData = 1;

  int m_writeGhostFieldData = 1;
};

} // namespace geos

#endif // GEOS_FILEIO_OUTPUTS_XDMFOUTPUT_HPP_
