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
 * @file LoggingObject.hpp
 */

#ifndef GEOS_COMMON_LOGGINGOUTPUT_HPP
#define GEOS_COMMON_LOGGINGOUTPUT_HPP

#include "common/LogMsg.hpp"

namespace geos
{ // TODO document
namespace logging
{

class LogOutput {
public:

  LogOutput() = default;
  LogOutput( LogOutput && other ) = default;
  LogOutput( LogOutput const & other ) = delete;
  LogOutput & operator=( LogOutput const & other ) = delete;
  LogOutput & operator=( LogOutput && other ) = delete;

  virtual void log( LogMsg message ) = 0;

};

/**
 * @brief LogOutput implementations
 */
///@{

class LogHDF5Output final : public LogOutput {
public:

  LogHDF5Output( string_view filePath );

  void log( LogMsg message ) override;

private:

  class HDF5Resources;

  std::unique_ptr< HDF5Resources > m_resources;

};


class LogTextStreamOutput final : public LogOutput {
public:

  LogTextOutput( string_view streamToOutputTo );

  void log( LogMsg message ) override;

private:

  std::ostream * m_outStream;
  
};


class LogTextFileOutput final : public LogOutput {
public:

  LogTextOutput( string_view rankFileOutputDir );

  void log( LogMsg message ) override;

private:

  std::ofstream m_rankFileStream;
  
};

///@}

} /* namespace logging */
} /* namespace geos */


#endif /* GEOS_COMMON_LOGGINGOUTPUT_HPP */