/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "LASFile.hpp"

namespace geosx
{

void LASLine::ParseLine( string const & line )
{
  m_keywordname.reserve( 8 );
  m_unit.reserve( 4 );
  m_data.reserve( 60 );
  m_description.reserve( 100 );

  TYPE curType = TYPE::MNEM;
  for( char const & c : line )
  {
    if( c == '.' && curType == TYPE::MNEM )
    {
      curType = TYPE::UNITS;
      continue;
    }
    else if( c == ' ' && curType == TYPE::UNITS )
    {
      curType = TYPE::DATA;
      continue;
    }
    else if( c == ':' && curType == TYPE::DATA )
    {
      curType = TYPE::DESCRIPTION;
      continue;
    }

    if( curType == TYPE::MNEM )
    {
      if( c == ' ') continue;
      m_keywordname.push_back(c);
    }
    else if( curType == TYPE::UNITS )
    {
      m_unit.push_back(c);
    }
    else if( curType == TYPE::DATA )
    {
      if( c == ' ') continue;
      m_data.push_back(c);
    }
    else if( curType == TYPE::DESCRIPTION )
    {
      m_data.push_back(c);
    }
  }
}

std::streampos LASSection::ParseSection( std::ifstream & file )
{
   string curLine;
   std::streampos pos = file.tellg();
   while ( std::getline(file, curLine) )
   {
     stringutilities::TrimLeft( curLine );
     if( curLine[0] == '#' ) continue;
     if( curLine[0] == '~' )                // We reach a new section
     {
       break;
     }
     pos = file.tellg();
     ParseLine( curLine );
   }
   return pos;
}

void LASInformationSection::WriteSection( std::ofstream & file ) const
{
  file << "~" << GetName() << " Section\n";
  this->forLines([&]( auto & line )
  {
    file << line.GetLine() << "\n";
  });
}

LASLine const & LASInformationSection::GetLine( string const & keyword ) const
{
  for( auto & lasLine : m_lines )
  {
    if( lasLine.GetKeyword() == keyword )
    {
      return lasLine;
    }
  }
  GEOSX_ERROR( "Keyword " << keyword << " not found in section " << GetName() );
  return m_lines[0]; // Should never be reached;
}

bool LASInformationSection::HasKeyword( string const & keyword ) const
{
  for( auto & lasLine : m_lines )
  {
    if( lasLine.GetKeyword() == keyword )
    {
      return true;
    }
  }
  return false;
}
localIndex LASWellInformationSection::GetNumberOfLogEntries() const
{
  real64 start = GetLine("STRT").GetDataAsReal64();
  real64 stop = GetLine("STOP").GetDataAsReal64();
  real64 step = GetLine("STEP").GetDataAsReal64();

  real64 length = stop - start;
  return std::round( length / step ) +1;
}

void LASInformationSection::ParseLine( string const & line )
{
  LASLine curLine( line );
  GEOSX_ERROR_IF( HasKeyword( curLine.GetKeyword() ) != 0, "Keyword " << curLine.GetKeyword()
                 << " was already defined in "<< GetName() );
  m_lines.push_back( curLine );
}

localIndex LASCurveInformationSection::FindLogIndex( string const & logName ) const
{
  localIndex logIndex = -1;
  for( auto & line : m_lines )
  {
    logIndex++;
    if( line.GetKeyword() == logName )
    {
      return logIndex;
    }
  }
  return logIndex;
}
void LASASCIILogDataSection::WriteSection( std::ofstream & file ) const
{
  file << "~" << GetName() << " Section\n";
  for( localIndex i = 0; i < m_nbLogEntries; i++ )
  {
    for( localIndex j = 0; j < m_nbCurves; j++ )
    {
      file << m_logs[j][i] << " ";
    }
    file << "\n";
  }
}

void LASASCIILogDataSection::ParseLine( string const & line )
{
  string_array splitLine = stringutilities::Tokenize( line, " \t\n\r" );
  for( integer i = 0; i < m_nbCurves; i++ )
  {
    m_logs[i][m_count] = std::stold( splitLine[i] );
  }
  m_count++;
}

void LASFile::Load( string const& fileName )
{
  std::ifstream file( fileName );
  GEOSX_ERROR_IF( !file.is_open(), "Can't open " << fileName );
  string curLine;
  while ( std::getline(file, curLine) )
  {
    stringutilities::TrimLeft( curLine );
    if( curLine[0] == '#' ) continue;  // Comment

    if( curLine[0] == '~' )            // Section
    {
      if( curLine[1] != 'A' )
      {
        std::unique_ptr< LASInformationSection > curLASInformationSection = LASInformationSection::CreateLASInformationSection( curLine[1] );
        std::streampos curPos = curLASInformationSection->ParseSection( file );
        curLASInformationSection->CheckKeywords();
        file.seekg( curPos );
        m_lasInformationSections.push_back( std::move( curLASInformationSection ) );
      }
      else
      {
        LASWellInformationSection * lastWellInformationSection = GetLastSection<LASWellInformationSection>();
        LASCurveInformationSection * lastCurveInformationSection = GetLastSection<LASCurveInformationSection>();
        LASASCIILogDataSection curLASASCIISection( lastWellInformationSection->GetNumberOfLogEntries(),
                                                   lastCurveInformationSection->GetNumberOfCurves() );
        std::streampos curPos = curLASASCIISection.ParseSection( file );
        m_lasASCIILogDataSection.push_back( curLASASCIISection );
        file.seekg( curPos );
      }
    }
  }
  file.close();
}

void LASFile::Save( string const& fileName ) const
{
  std::ofstream file( fileName );
  file << "# LAS Log file written by GEOSX" << "\n";
  localIndex countLog = 0;
  set< string > sectionsOutputed;

  this->forInformationSections( [&]( auto & informationSection )
  {
    if( sectionsOutputed.count( informationSection->GetName() ) )
    {
      m_lasASCIILogDataSection[countLog++].WriteSection( file );
      sectionsOutputed.clear();
    }
    informationSection->WriteSection( file );
    sectionsOutputed.insert( informationSection->GetName() );
  });
  if( countLog == 0 )
  {
    GEOSX_ASSERT( m_lasASCIILogDataSection.size() == 1 );
    m_lasASCIILogDataSection[0].WriteSection( file );
  }
  else
  {
    m_lasASCIILogDataSection[countLog].WriteSection( file );
  }
  file.close();
}

arraySlice1d< real64> LASFile::GetLog( string const & logName ) const
{
  if( !HasLog( logName ) )
  {
    GEOSX_ERROR( logName << " not found in LAS file");
  }
  localIndex logSectionIndex = 0;
  for( auto const & section :  m_lasInformationSections )
  {
    if( section->GetName() == LASCurveInformationSection::GetNameStatic() )
    {
      LASCurveInformationSection * lasCurve = dynamic_cast< LASCurveInformationSection* > ( section.get() );
      if( lasCurve->HasKeyword( logName ) )
      {
        localIndex logIndex = lasCurve->FindLogIndex( logName );
        GEOSX_ASSERT( logIndex > -1 );
        return m_lasASCIILogDataSection[logSectionIndex].GetLog( logIndex );
      }
      else
      {
        logSectionIndex++;
      }
    }
  }
  return m_lasASCIILogDataSection[0].GetLog(0); // should never be reached
}


localIndex LASFile::LogSize( string const& logName ) const
{
  if( !HasLog( logName ) )
  {
    GEOSX_ERROR( logName << " not found in LAS file");
  }
  localIndex logSectionIndex = 0;
  for( auto const & section : m_lasInformationSections )
  {
    if( section->GetName() == LASCurveInformationSection::GetNameStatic() )
    {
      LASCurveInformationSection * lasCurve = dynamic_cast< LASCurveInformationSection* > ( section.get() );
      if( lasCurve->HasKeyword( logName ) )
      {
        GEOSX_ASSERT( lasCurve->FindLogIndex( logName ) > -1 );
        return m_lasASCIILogDataSection[logSectionIndex].LogSize();
      }
      else
      {
        logSectionIndex++;
      }
    }
  }
  return 0; // should never be reached
}

bool LASFile::HasLog( string const & logName ) const
{
  for( auto const & section : m_lasInformationSections )
  {
    if( section->GetName() == LASCurveInformationSection::GetNameStatic() )
    {
      if ( section->HasKeyword( logName ) )
      {
        return true;
      }
    }
  }
  return false;
}


template< typename T >
T * LASFile::GetLastSection()
{
  for( auto const & section : m_lasInformationSections )
  {
    if( section->GetName() == T::GetNameStatic() )
    {
      return dynamic_cast< T * >( section.get() );
    }
  }
  GEOSX_ERROR(" LAS Log file  is not valid: Log Data Section was " <<
             " declared before the " << T::GetNameStatic() << " section");
  return nullptr;
}

std::unique_ptr< LASInformationSection > LASInformationSection::CreateLASInformationSection( char const & name )
{
  if( name == 'V' )          // Version Information
  {
    return std::unique_ptr < LASVersionInformationSection > ( new LASVersionInformationSection() );
  }
  else if( name == 'W' )     // Well Information
  {
    return std::unique_ptr< LASWellInformationSection > ( new LASWellInformationSection() );
  }
  else if( name == 'C' )     // Curve Information
  {
    return std::unique_ptr< LASCurveInformationSection > ( new LASCurveInformationSection() );
  }
  else if( name == 'P' )     // Parameter Information
  {
    return std::unique_ptr< LASParameterInformationSection > ( new LASParameterInformationSection() );
  }
  else if( name == 'O' )     // Other Information
  {
    return std::unique_ptr< LASOtherInformationSection > ( new LASOtherInformationSection() );
  }
  else
  {
    GEOSX_ERROR( name << " is not a valid section for LAS files" );
    return nullptr;
  }
}

}
