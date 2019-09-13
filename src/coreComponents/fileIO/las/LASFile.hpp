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

#ifndef GEOSX_SRC_CORECOMPONENTS_FILEIO_LAS_LASFILE_
#define GEOSX_SRC_CORECOMPONENTS_FILEIO_LAS_LASFILE_

#include "common/DataTypes.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

class LASVersionInformationSection;

/*!
 * @brief This class contains a parsed line of a LAS file
 * @details A LAS line have always this structure
 * MNEM.UNIT    DATA  :  DESCRIPTION
 * With '.' as a first delimiter, the first space as the second
 * delimiter and ':' as the third delimiter
 */
class LASLine
{
  public:
    LASLine() = default;

    LASLine( string const& line )
    {
      ParseLine( line );
    }

    string const & GetKeyword() const
    {
      return m_keywordname;
    }

    string const & GetData() const
    {
      return m_data;
    }

    bool GetDataAsBool() const
    {
      GEOS_ASSERT( m_data == "YES" || m_data == "NO");
      return ( m_data == "YES" ? true : false );
    }

    localIndex GetDataAsLocalIndex() const
    {
      return stringutilities::fromString< localIndex >( m_data );
    }

    real64 GetDataAsReal64() const
    {
      return stringutilities::fromString< real64 >( m_data );
    }

    string const & GetDescription() const
    {
      return m_description;
    }

    string const & GetUnit() const
    {
      return m_unit;
    }

    bool HasUnit() const
    {
      return !m_unit.empty();
    }

    /*
     * @brief Build a string containing the line
     */
    string GetLine() const
    {
      return m_keywordname + "." + m_unit + "    " + m_data + "    :    " + m_description;
    }
  private:
    void ParseLine( string const & line )
    {
      // First get the keyword and the rest of the line
      string_array keywordAndRest = stringutilities::Split( line, "." );
      stringutilities::RemoveSpaces( keywordAndRest[0] );
      m_keywordname = keywordAndRest[0];

      // Second get the unit and the rest of the line
      string_array unitsAndRest = stringutilities::Split( keywordAndRest[1], " ");
      stringutilities::RemoveSpaces( unitsAndRest[0] );
      m_unit = unitsAndRest[0];

      // Third get the value and the rest of te line
      string_array valueAndRest = stringutilities::Split( unitsAndRest[1], ":" );
      stringutilities::Trim( valueAndRest[0] );
      m_data = valueAndRest[0];

      // Finally, get the description
      stringutilities::Trim( valueAndRest[1] );
      m_description = valueAndRest[1];
    }
  private:
    /// Name of the keyword
    string m_keywordname;

    /// Unit (if applicable);
    string m_unit;

    /// Data
    string m_data;

    /// Description
    string m_description;
};

class LASSection
{
  public:
    LASSection() = default;
    LASSection( LASSection const &  ) = default;
    virtual ~LASSection()
    {
    }

    virtual std::streampos ParseSection( std::ifstream & file )
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

    virtual void WriteSection( std::ofstream & file ) const = 0;

    virtual string const GetName() const = 0;

  protected:
    virtual void ParseLine( string const & line ) = 0;
};
/*!
 * @brief Basis class for a LAS Information Section
 * @details a LAS Information Section contains an ensemble of LASLine,
 * referenced by their keywords
 */
class LASInformationSection : public LASSection
{
  public:
    LASInformationSection() = default;

    /*!
     * Check if all the mandatory keywords are registered
     */
    virtual void CheckKeywords()
    {
      for( string & keyword : m_mandatoryKeyword )
      {
        GEOS_ERROR_IF( HasKeyword( keyword ) == 0, "Mandatory keyword " << keyword << " not found in "
                                                      << GetName() );
      }
    }

    static LASInformationSection * CreateLASInformationSection( char const & name );

    virtual void WriteSection( std::ofstream & file ) const override
    {
      file << "~" << GetName() << " Section\n";
      for( auto & line: m_lines )
      {
        file << line.GetLine() << "\n";
      }
    }

    LASLine const & GetLine( string const & keyword ) const
    {
      for( auto & lasLine : m_lines )
      {
        if( lasLine.GetKeyword() == keyword )
        {
          return lasLine;
        }
      }
      GEOS_ERROR( "Keyword " << keyword << " not found in section " << GetName() );
      return m_lines[0]; // Should never be reached;
    }

    bool HasKeyword( string const & keyword ) const
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

  protected:
    virtual void ParseLine( string const & line ) override
    {
      LASLine curLine( line );
      GEOS_ERROR_IF( HasKeyword( curLine.GetKeyword() ) != 0, "Keyword " << curLine.GetKeyword()
                     << " was already defined in "<< GetName() );
      m_lines.push_back( curLine );
    }


  protected:
    /// Contains the mandatory keyword for the section
    array1d< string > m_mandatoryKeyword;

    /// Contains all the lines (A line = keyword + value(s))
    std::vector< LASLine > m_lines;
};

/*!
 * @brief This section identify the LAS format ans whether the wrap format is used
 */
class LASVersionInformationSection : public LASInformationSection
{
  public:
    LASVersionInformationSection() :
      LASInformationSection()
    {
      m_mandatoryKeyword.push_back( "VERS" );
      m_mandatoryKeyword.push_back( "WRAP" );
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Version Information";
    }

};

/*!
 * @brief This section identify the well
 */
class LASWellInformationSection : public LASInformationSection
{
  public:
    LASWellInformationSection() :
      LASInformationSection()
    {
      m_mandatoryKeyword.push_back( "STRT" );
      m_mandatoryKeyword.push_back( "STOP" );
      m_mandatoryKeyword.push_back( "STEP" );
      m_mandatoryKeyword.push_back( "NULL" );
      m_mandatoryKeyword.push_back( "COMP" );
      m_mandatoryKeyword.push_back( "WELL" );
      m_mandatoryKeyword.push_back( "FLD" );
      m_mandatoryKeyword.push_back( "LOC" );
      m_mandatoryKeyword.push_back( "PROV" );
      m_mandatoryKeyword.push_back( "SRVC" );
      m_mandatoryKeyword.push_back( "DATE" );
      m_mandatoryKeyword.push_back( "UWI" );
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Well Information";
    }

    localIndex GetNumberOfLogEntries() const
    {
      real64 start = GetLine("STRT").GetDataAsReal64();
      real64 stop = GetLine("STOP").GetDataAsReal64();
      real64 step = GetLine("STEP").GetDataAsReal64();

      real64 length = stop - start;
      return static_cast< localIndex > ( length / step ) + 1;
    }
};

/*!
 * @brief This section describes the curves and their units
 */
class LASCurveInformationSection : public LASInformationSection
{
  public:
    LASCurveInformationSection() :
      LASInformationSection()
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Curve Information";
    }

    /*!
     * @brief Returns the number of curves
     * @details DEPTH, DEPT, TIME or INDEX is one the curve
     * A curve is a log
     */
    localIndex GetNumberOfCurves() const
    {
      return m_lines.size();
    }
  private:
    /*!
     * @brief For this specific section, one keyword must be
     * DEPTH, DEPT, TIME or INDEX"
     * */
    virtual void CheckKeywords() override
    {
      GEOS_ERROR_IF( !HasKeyword( "DEPTH" ) && !HasKeyword( "DEPT" ) && !HasKeyword( "TIME" )  && !HasKeyword( "INDEX "),
                     "Invalid " << GetName() << " section. It musts contains at least one those keyword \n"
                     << "DEPTH DEPT TIME INDEX");
    }
};

/*!
 * @brief This optional section describes optional parameters
 */
class LASParameterInformationSection : public LASInformationSection
{
  public:
    LASParameterInformationSection() :
      LASInformationSection()
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Parameter Information";
    }
};

/*!
 * @brief This optional section describes some comments
 */
class LASOtherInformationSection : public LASInformationSection
{
  public:
    LASOtherInformationSection() :
      LASInformationSection()
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "Other Information";
    }
  private:
    virtual void ParseLine( string const & line ) override
    {
      m_comments += line;
    }

    /*!
     * @brief For this specific section there is no keyword
     * */
    virtual void CheckKeywords() override
    {
      GEOS_ERROR_IF( !m_lines.empty(),
                     "Invalid " << GetName() << " section. No keyword should be defined, only data.");
    }

    virtual void WriteSection( std::ofstream & file ) const override
    {
      file << "~" <<GetName() << " Section\n";
      file << m_comments << "\n";
    }

  private:
    /// Content of the section
    string m_comments;

};

/*!
 * @brief This optional section contains all the log data
 */
class LASASCIILogDataSection : public LASSection
{
  public:
    LASASCIILogDataSection( localIndex nbEntries, localIndex nbCurves ) :
      LASSection(),
      m_logs(nbCurves, nbEntries),
      m_nbCurves( nbCurves ),
      m_nbLogEntries( nbEntries ),
      m_count(0)
    {
    }

    virtual string const GetName() const override
    {
      return GetNameStatic();
    }

    static string GetNameStatic() 
    {
      return "ASCII Log Data";
    }

    virtual void WriteSection( std::ofstream & file ) const override
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
  private:
    virtual void ParseLine( string const & line ) override
    {
      string_array splitLine = stringutilities::Tokenize( line, " \t\n\r" );
      GEOS_ASSERT( splitLine.size() == m_nbCurves );
      for( integer i = 0; i < splitLine.size(); i++ )
      {
        m_logs[i][m_count] = stringutilities::fromString< real64 >( splitLine[i] );
      }
      m_count++;
    }
  private:
    /// Contains the well logs
    array2d< real64 > m_logs;

    localIndex m_nbCurves;

    localIndex m_nbLogEntries;

    localIndex m_count;
};

class LASFile
{
  public:
    LASFile() = default;

    ~LASFile()
    {
      for( auto lasInformationSection : m_lasInformationSections )
      {
        delete lasInformationSection;
      }
      m_lasInformationSections.clear();
    }
    
    void Load( string const& fileName)
    {
      std::ifstream file( fileName );
      GEOS_ERROR_IF( !file.is_open(), "Can't open " << fileName );
      string curLine;
      while ( std::getline(file, curLine) )
      {
        stringutilities::TrimLeft( curLine );
        if( curLine[0] == '#' ) continue;  // Comment

        if( curLine[0] == '~' )            // Section
        {
          if( curLine[1] != 'A' )
          {
            LASInformationSection * curLASInformationSection = LASInformationSection::CreateLASInformationSection( curLine[1] );
            std::streampos curPos = curLASInformationSection->ParseSection( file );
            curLASInformationSection->CheckKeywords();
            file.seekg( curPos );
            m_lasInformationSections.push_back( std::move( curLASInformationSection ) );
          }
          else
          {
            std:: cout << "AA" << std::endl;
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
      std::cout << m_lasASCIILogDataSection.size() << std::endl;
      file.close();
    }
    
    void Save( string const& fileName ) const
    {
      std::ofstream file( fileName );
      file << "# LAS Log file written by GEOSX" << "\n";
      localIndex countLog = 0;
      string toto;
      set< string > sectionsOutputed;
      for( auto informationSection : m_lasInformationSections )
      {
        if( sectionsOutputed.count( informationSection->GetName() ) )
        {
          m_lasASCIILogDataSection[countLog++].WriteSection( file );
          sectionsOutputed.clear();
        }
        informationSection->WriteSection( file );
        sectionsOutputed.insert( informationSection->GetName() );
      }
      std::cout << "count log " << countLog << std::endl;
      if( countLog == 0 )
      {
        GEOS_ASSERT( m_lasASCIILogDataSection.size() == 1 );
        m_lasASCIILogDataSection[0].WriteSection( file );
      }
      else
      {
        m_lasASCIILogDataSection[countLog].WriteSection( file );
      }
      file.close();
    }
  private:
    template< typename T >
    T * GetLastSection()
    {
      for( auto itr = m_lasInformationSections.rbegin(); itr != m_lasInformationSections.rend(); itr++)
      {
        if( (*itr)->GetName() == T::GetNameStatic() )
        {
          return dynamic_cast< T * >( *itr);
        }
      }
      GEOS_ERROR(" LAS Log file  is not valid: Log Data Section was " <<
                 " declared before the " << T::GetNameStatic() << " section");
      return nullptr;
    }
  private:
    /// All information sections of the LAS file
    std::vector< LASInformationSection * > m_lasInformationSections;

    /// Log data sections
    std::vector< LASASCIILogDataSection > m_lasASCIILogDataSection;
};
}

#endif
