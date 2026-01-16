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
 * @file LoggerMsgReportData.cpp
 */

#include "LoggerMsgReportData.hpp"
#include "common/format/EnumStrings.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"
#include "common/logger/MsgType.hpp"
#include <sstream>

namespace geos
{

void LoggerMsgReportData::increment( LogPart::Type logPartName, MsgType msgType )
{
  if( numMsgByPart.count( logPartName ) == 0 )
  {
    NumMsg numMsg{ { {msgType, 1}}, {{msgType, 1} }};
    numMsgByPart[logPartName] = numMsg;
  }
  else
  {
    NumMsg & numMsg = numMsgByPart.at( logPartName );
    numMsg.numMsg[msgType]++;
    numMsg.numMsgLoc[msgType]++;
  }
}

template<>
string TableTextFormatter::toString< LoggerMsgReportData >( LoggerMsgReportData const & report ) const
{
  TableLayout tableLayoutPerSection;
  TableData dataPerSection;

  std::ostringstream oss;
  tableLayoutPerSection.addColumn( " Types " );

  for( auto const & [ logPartName, numMsg ] : report.numMsgByPart )
  {
    tableLayoutPerSection.addColumn( EnumStrings< LogPart::Type >::toString( logPartName ) );
  }


  for( size_t i  = (size_t) MsgType::Error; i  != (size_t)MsgType::Undefined; i++ )
  {
    MsgType const currentType = (MsgType) i;
    stdVector< TableData::CellData > row;

    row.push_back( {CellType::Value, EnumStrings< MsgType >::toString( (MsgType) i ) } );
    for( auto const & [ _, msgTypes ] : report.numMsgByPart )
    {
      auto it =  msgTypes.numMsg.find( currentType );
      int const count = ( it != msgTypes.numMsg.end() ) ? it->second : 0;
      row.push_back( {CellType::Value, std::to_string( count )} );

    }
    dataPerSection.addRow( row );
  }

  TableTextFormatter textFormatter( tableLayoutPerSection );
  oss << textFormatter.toString( dataPerSection ) << "\n";

  return oss.str();
}

};
