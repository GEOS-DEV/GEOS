#include "LoggerMsgReportData.hpp"
#include "common/format/EnumStringsCore.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"
#include "common/logger/MsgType.hpp"
#include <sstream>

namespace geos
{

void LoggerMsgReportData::increment( string const & logPartName, MsgType msgType )
{
  std::cout << "DEBUG-- "<< logPartName<< " " << EnumStringsCore< MsgType >::toRawString( msgType ) << std::endl;
  if( numMsgByPart.count( logPartName ) ==0 )
  {
    NumMsg numMsg{ { {msgType, 1}}, {{msgType, 1} }};
    numMsgByPart[logPartName] = numMsg;
  }
  else
  {
    std::cout << "check "<< std::endl;
    NumMsg & numMsg = numMsgByPart.at(logPartName);
    numMsg.numMsg[msgType]++;
    numMsg.numMsgLoc[msgType]++;
  }
}

template<>
string TableTextFormatter::toString< LoggerMsgReportData >( LoggerMsgReportData const & report ) const
{
  TableLayout tableLayoutPerSection;
  TableData dataPerSection;
  /// {Numerical Methods : { Warning : 4 }, { Exception : 2 } }
  stdMap< string, stdMap< string, int > > cells;

  std::ostringstream oss;
  tableLayoutPerSection.addColumn( " Types " );

  for( auto const & [ logPartName, numMsg ] : report.numMsgByPart )
  {
    tableLayoutPerSection.addColumn( logPartName );
    for( auto const & [ msgType, count ] : numMsg.numMsg )
    {
      string const msgTypeStr = EnumStringsCore< MsgType >::toRawString( msgType );
      cells.get_inserted( logPartName ).get_inserted( msgTypeStr ) = count;
    }

    for( size_t fooInt = (size_t) MsgType::Error; fooInt != (size_t)MsgType::Undefined; fooInt++ )
    {
      string type =  EnumStringsCore< MsgType >::toRawString( (MsgType) fooInt );
      if( cells.get_inserted( logPartName ).count( type ) == 0 )
        cells.get_inserted( logPartName ).get_inserted( type ) = 0;
    }
  }

  for( size_t fooInt = (size_t) MsgType::Error; fooInt != (size_t)MsgType::Undefined; fooInt++ )
  {
    stdVector< TableData::CellData > countPerLogPart;
    string const type=  EnumStringsCore< MsgType >::toRawString( (MsgType) fooInt );
    countPerLogPart.push_back( {CellType::Value, type } );
    for( auto const & [ logPartName, msgTypes ] : cells )
    {
      if( msgTypes.count( type ) != 0 )
      {
        countPerLogPart.push_back( {CellType::Value, std::to_string( msgTypes.at( type ) )} );
      }
      else
      {
        countPerLogPart.push_back( {CellType::Value, 0} );
      }
    }
    dataPerSection.addRow( countPerLogPart );
  }

  TableTextFormatter textFormatter( tableLayoutPerSection );
  oss << textFormatter.toString( dataPerSection ) << "\n";

  return oss.str();
}

};
