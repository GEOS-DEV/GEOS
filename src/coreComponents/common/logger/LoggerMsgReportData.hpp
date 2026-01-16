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
 * @file LoggerMsgReportData.hpp
 */

#ifndef GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP
#define GEOS_COMMON_LOGGER_MSG_REPORT_DATA_HPP

#include "common/format/LogPart.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "MsgType.hpp"

namespace geos
{

struct NumMsg;

struct LoggerMsgReportData
{
  stdMap< LogPart::Type, NumMsg > numMsgByPart;

  void increment( LogPart::Type logPartName, MsgType );
};

struct NumMsg
{
  stdMap< MsgType, int > numMsg;
  stdMap< MsgType, int > numMsgLoc;
};

template<>
string TableTextFormatter::toString< LoggerMsgReportData >( LoggerMsgReportData const & ) const;

}

#endif
