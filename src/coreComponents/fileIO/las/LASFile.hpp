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
      GEOSX_ASSERT( m_data == "YES" || m_data == "NO");
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
    void ParseLine( string const & line );

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

    virtual std::streampos ParseSection( std::ifstream & file );

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
        GEOSX_ERROR_IF( HasKeyword( keyword ) == 0, "Mandatory keyword " << keyword << " not found in "
                                                      << GetName() );
      }
    }

    static std::unique_ptr< LASInformationSection > CreateLASInformationSection( char const & name );

template< typename LAMBDA >
void forLines( LAMBDA&& lambda )
{
  for( auto & line: m_lines )
  {
    lambda( line );
  }
}

template< typename LAMBDA >
void forLines( LAMBDA&& lambda ) const
{
  for( auto & line: m_lines )
  {
    lambda( line );
  }
}

    virtual void WriteSection( std::ofstream & file ) const override;

    LASLine const & GetLine( string const & keyword ) const;

    bool HasKeyword( string const & keyword ) const;

  protected:
    virtual void ParseLine( string const & line ) override;

  protected:
    /// Contains the mandatory keyword for the section
    array1d< string > m_mandatoryKeyword;

    /// Contains all the lines (A line = keyword + value(s))
    array1d< LASLine > m_lines;
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

    localIndex GetNumberOfLogEntries() const;
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

    /*!
     * @brief returns the index of the log referenced in this InformationSection
     * @returns -1 if the log is not found, > -1 value otherwise.
     */
    localIndex FindLogIndex( string const & logName ) const;
  private:
    /*!
     * @brief For this specific section, one keyword must be
     * DEPTH, DEPT, TIME or INDEX"
     * */
    virtual void CheckKeywords() override
    {
      GEOSX_ERROR_IF( !HasKeyword( "DEPTH" ) && !HasKeyword( "DEPT" ) && !HasKeyword( "TIME" )  && !HasKeyword( "INDEX "),
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
      GEOSX_ERROR_IF( !m_lines.empty(),
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
      m_logs(nbCurves, nbEntries ),
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

    virtual void WriteSection( std::ofstream & file ) const override;

    /*!
     * @brief returns the size of the stored logs
     */
    localIndex LogSize() const
    {
      return m_nbLogEntries;
    }

    /*!
     * @brief returns the number of logs
     */
    localIndex NbLogs() const
    {
      return m_nbCurves;
    }

    /*!
     * @brief returns the ith log
     * @param[in] i  index of the log
     */
    arraySlice1d< real64 > GetLog( localIndex i ) const
    {
      return m_logs[i];
    }
    
  private:
    virtual void ParseLine( string const & line ) override;

  private:
    /// Contains the well logs
    array2d< real64 > m_logs;

    localIndex m_nbCurves;

    localIndex m_nbLogEntries;

    /// Use for internal counting while parsing the file
    localIndex m_count;
};

class LASFile
{
  public:
    LASFile() = default;

    ~LASFile(){}
    
    void Load( string const& fileName );
    
    void Save( string const& fileName ) const;

    arraySlice1d< real64> GetLog( string const & logName ) const;

    template< typename LAMBDA >
    void forInformationSections( LAMBDA&& lambda )
    {
      for( auto& informationSection : m_lasInformationSections )
      {
        lambda( informationSection );
      }
    }

    template< typename LAMBDA >
    void forInformationSections( LAMBDA&& lambda ) const
    {
      for( auto& informationSection : m_lasInformationSections )
      {
         lambda( informationSection );
      }
    }

    template< typename LAMBDA >
    void forLogSections( LAMBDA&& lambda )
    {
      for( auto& logSection : m_lasASCIILogDataSection )
      {
        lambda( logSection );
      }
    }

     template< typename LAMBDA >
     void forLogSections( LAMBDA&& lambda ) const
     {
      for( auto& logSection : m_lasASCIILogDataSection )
      {
        lambda( logSection );
      }
    }

    LASInformationSection const & GetInformationSection( integer i ) const
    {
      return *m_lasInformationSections[i];
    }

    LASASCIILogDataSection const &GetLogSection( integer i ) const
    {
      return m_lasASCIILogDataSection[i];
    }

    /*!
     * @brief returns the number of log entries for a property
     * @param[in] logName the name of the property
     */
    localIndex LogSize( string const& logName ) const;

    /*!
     * @brief Tells wether or not the LAS file contains a log.
     * @param[in] logName the name the log
     */
    bool HasLog( string const & logName ) const;

    /*!
     * @brief Tells wether or not the LAS section contains a line
     * @param[in] keyword the keyword of the line
     */
    template< typename SECTION >
    bool HasLine( string const & keyword ) const
    {
      for( auto const & section : m_lasInformationSections )
      {
        if( section->GetName() == SECTION::GetNameStatic() )
        {
          if( section->HasKeyword( keyword ) )
          {
            return true;
          }
        }
      }
       return false;
     }

    /*!
     * @brief returns all the lines corresponding to a keyword in specific sections(s)
     * @details a vector is returned as it is possible to have plenty of same sections in some LAS files
     */
     template< typename T >
     std::vector< LASLine const * > const GetLASLines( string const& keyword ) const
     {
      std::vector< LASLine const * > result;
      for( auto const & section : m_lasInformationSections )
      {
        if( section->GetName() == T::GetNameStatic() )
        {
          result.push_back( &section->GetLine( keyword ) );
        }
      }
      GEOSX_ASSERT( result.size() > 0 );
      return result;
    }

  private:
    template< typename T >
    T * GetLastSection();
  private:
    /// All information sections of the LAS file
    array1d< std::unique_ptr< LASInformationSection > > m_lasInformationSections;

    /// Log data sections
    array1d< LASASCIILogDataSection > m_lasASCIILogDataSection;
};
}

#endif
