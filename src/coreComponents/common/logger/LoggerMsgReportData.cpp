#include "LoggerMsgReportData.hpp"
#include "common/format/EnumStringsCore.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/logger/MsgType.hpp"
#include <sstream>

namespace geos
{

void LoggerMsgReportData::increment( string const & logPartName, MsgType msgType )
{
  if( numMsgByPart.count( logPartName ) ==0 )
  {
    NumMsg numMsg{ { {msgType, 1}}, {{msgType, 1} }};
    numMsgByPart[logPartName] = numMsg;
  }
  else
  {
    NumMsg numMsg = numMsgByPart[logPartName];
    numMsg.numMsg[msgType]++;
    numMsg.numMsgLoc[msgType]++;
  }
}

template<>
string TableTextFormatter::toString< LoggerMsgReportData >( LoggerMsgReportData const & report ) const
{
  stdMap< string, TableLayout > tableLayoutPerSection;
  stdMap< string, TableData > dataPerSection;

  std::ostringstream oss;

  for( auto const & [ logPartName, numMsg ] : report.numMsgByPart )
  {

    TableLayout layout;
    TableData tableData;
    stdVector< TableData::CellData > cells;
    layout.setTitle( logPartName );
    for( auto const & [ msgType, count ] : numMsg.numMsg )
    {
      layout.addColumn( EnumStringsCore< MsgType >::toRawString( msgType ));
      cells.push_back( {CellType::Value, std::to_string( count )} );
    }
    tableData.addRow( cells );
    tableLayoutPerSection.get_inserted( logPartName ) = layout;
    dataPerSection.get_inserted( logPartName ) = tableData;
  }

  for( auto const & [ logPartName, layout ] : tableLayoutPerSection )
  {
    auto const & data = dataPerSection[logPartName];
    TableTextFormatter textFormatter( layout );
    oss << textFormatter.toString( data ) << "\n"; 
  }

  return oss.str();
}

};
