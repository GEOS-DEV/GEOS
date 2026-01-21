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
 * @file LogHistory.hpp
 */

#ifndef GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP
#define GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP

#include "common/StdContainerWrappers.hpp"
#include "common/format/LogPart.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "DiagnosticMessage.hpp"


namespace geos
{

struct NumMsg
{
  string fileLocation;
  integer codeLocation;
  integer count;
};

class LogHistory
{
public:

  void notifyMsg( LogPart::Type logPartName, DiagnosticMsg const & diagMsg, integer threadCount );

  auto const & get() const
  { return m_messageCounts; }

private:
  stdMap< LogPart::Type, stdMap< MsgType, stdVector< NumMsg > > > m_messageCounts;

};
template<>
string TableTextFormatter::toString< LogHistory >( LogHistory const & ) const;

}

#endif
