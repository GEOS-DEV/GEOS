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
#include <string>

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

  /*!
   * @brief default constructor
   */
  LASLine() = default;

  /*!
   * @brief default destructor
   */
  ~LASLine() = default;

  /*!
   * @brief instantiate the class using the string containing tha actual line
   * @param[in] line the input line
   */
  LASLine( string const & line )
  {
    ParseLine( line );
  }

  /*!
   * @brief Get the keyword of the line
   * @return the keyword of the line
   */
  string const & GetKeyword() const
  {
    return m_keywordname;
  }

  /*!
   * @brief Get the data associated with the line as a string
   * @return the data associatied with the line
   */
  string const & GetData() const
  {
    return m_data;
  }

  /*!
   * @brief Get the data associated with the line as a bool
   * @return the data associatied with the line as a bool
   */
  bool GetDataAsBool() const
  {
    GEOSX_ASSERT( m_data == "YES" || m_data == "NO" );
    return ( m_data == "YES" ? true : false );
  }

  /*!
   * @brief Get the data associated with the line as a localIndex
   * @return the data associatied with the line as a localIndex
   */
  localIndex GetDataAsLocalIndex() const
  {
    return std::stoi( m_data );
  }

  /*!
   * @brief Get the data associated with the line as a real64
   * @return the data associatied with the line as a real64
   */
  real64 GetDataAsReal64() const
  {
    return std::stold( m_data );
  }

  /*!
   * @brief Get the description of the line
   * @return the description of the line
   */
  string const & GetDescription() const
  {
    return m_description;
  }

  /*!
   * @brief Get the unit of the data of the line as a string
   * @return the unit of the data of the line as a string
   */
  string const & GetUnit() const
  {
    return m_unit;
  }

  /*!
   * @brief Tell wether or not the line is associated with an unit
   * @return true if the line is associated with an unit, false otherwise
   */
  bool HasUnit() const
  {
    return !m_unit.empty();
  }

  /*!
   * @brief Build a string containing the line
   * @return a string with the line content
   */
  string GetLine() const
  {
    return m_keywordname + "." + m_unit + "    " + m_data + "    :    " + m_description;
  }
private:
  /*!
   * @brief Parse the line an decompose it into a keyword, an unit, data and description
   */
  void ParseLine( string const & line );

private:
  /*!
   * @brief enum descriping the different type of information stored in a line
   */
  enum struct TYPE
  {
    MNEM,
    UNITS,
    DATA,
    DESCRIPTION
  };
  /// Name of the keyword
  string m_keywordname;

  /// Unit (if applicable);
  string m_unit;

  /// Data
  string m_data;

  /// Description
  string m_description;
};

/*!
 * @brief base class defining a LAS Section
 */
class LASSection
{
public:

  /*!
   * @brief default constructor
   */
  LASSection() = default;

  /*!
   * @brief default copy constructor
   */
  LASSection( LASSection const & ) = default;

  /*!
   * @brief default destructor
   */
  virtual ~LASSection()
  {}

  /*!
   * @brief Parse a whole LAS Section
   * @param[in] file the stream to the input file
   * @return the position of the cursor while reading the LAS File
   */
  virtual std::streampos ParseSection( std::ifstream & file );

  /*!
   * @brief Write the section to a file
   * @param[in] file the stream to the output file
   */
  virtual void WriteSection( std::ofstream & file ) const = 0;

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const = 0;

protected:
  /*!
   * @brief Parse the line of a section
   */
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
  /*!
   * @brief default constructor
   */
  LASInformationSection() = default;

  /*!
   * @brief Check if all the mandatory keywords are registered according to the LAS standard
   */
  virtual void CheckKeywords()
  {
    for( string & keyword : m_mandatoryKeyword )
    {
      GEOSX_ERROR_IF( HasKeyword( keyword ) == 0, "Mandatory keyword " << keyword << " not found in "
                                                                       << GetName() );
    }
  }

  /*!
   * @brief Create an instance of a LASInformationSection
   * @param[in] name first character of the name f the name of the section
   * @return an unique_ptr to the LASInformationSection
   */
  static std::unique_ptr< LASInformationSection > CreateLASInformationSection( char const & name );

  /*!
   * @brief Loop over the lines of the information section and apply a \p lambda funtion
   * @param[in] lambda the lambda function
   */
  template< typename LAMBDA >
  void forLines( LAMBDA && lambda )
  {
    for( auto & line: m_lines )
    {
      lambda( line );
    }
  }

  template< typename LAMBDA >

  /*!
   * @copydoc forLines( LAMBDA && )
   */
  void forLines( LAMBDA && lambda ) const
  {
    for( auto & line: m_lines )
    {
      lambda( line );
    }
  }

  /*!
   * @brief Writes the whole section
   * @param[in] file the file to be output
   */
  virtual void WriteSection( std::ofstream & file ) const override;

  /*!
   * @brief Get a line associated with the \p keyword
   * @param[in] keyword the keyword associated to the line to get
   * @return the LASLine associated with the keyword
   */
  LASLine const & GetLine( string const & keyword ) const;

  /*!
   * @brief check if the \p keyword is associated with a line for this section
   * @return true if the keyword is associated with a line in this section. false otherwise
   */
  bool HasKeyword( string const & keyword ) const;

protected:
  /*!
   * @brief Parse a line for the information sections
   * @param[in] line the line to be parsed
   */
  virtual void ParseLine( string const & line ) override;

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
  /*!
   * @brief default constructor
   */
  LASVersionInformationSection():
    LASInformationSection()
  {
    m_mandatoryKeyword.push_back( "VERS" );
    m_mandatoryKeyword.push_back( "WRAP" );
  }

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const override
  {
    return GetNameStatic();
  }

  /*!
   * @copydoc GetName()
   */
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

  /*!
   * @brief default constructor
   */
  LASWellInformationSection():
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

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const override
  {
    return GetNameStatic();
  }

  /*!
   * @copydoc GetName()
   */
  static string GetNameStatic()
  {
    return "Well Information";
  }

  /*!
   * @brief Get the number of log entries
   * @details the number of log entries is computed using the information
   * within START, STOP and STEP keywords
   * @return the number of log entries
   */
  localIndex GetNumberOfLogEntries() const;
};

/*!
 * @brief This section describes the curves and their units
 */
class LASCurveInformationSection : public LASInformationSection
{
public:
  /*!
   * @brief default constructor
   */
  LASCurveInformationSection():
    LASInformationSection()
  {}

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const override
  {
    return GetNameStatic();
  }

  /*!
   * @copydoc GetName()
   */
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
    GEOSX_ERROR_IF( !HasKeyword( "DEPTH" ) && !HasKeyword( "DEPT" ) && !HasKeyword( "TIME" )  && !HasKeyword( "INDEX " ),
                    "Invalid " << GetName() << " section. It musts contains at least one those keyword \n"
                               << "DEPTH DEPT TIME INDEX" );
  }
};

/*!
 * @brief This optional section describes optional parameters
 */
class LASParameterInformationSection : public LASInformationSection
{
public:
  /*!
   * @brief default constructor
   */
  LASParameterInformationSection():
    LASInformationSection()
  {}

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const override
  {
    return GetNameStatic();
  }

  /*!
   * @copydoc GetName()
   */
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
  /*!
   * @brief default constructor
   */
public:
  LASOtherInformationSection():
    LASInformationSection()
  {}

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const override
  {
    return GetNameStatic();
  }

  /*!
   * @copydoc GetName()
   */
  static string GetNameStatic()
  {
    return "Other Information";
  }
private:
  /*!
   * @brief For this specific section, just add the line to the list of lines
   * @param[in] line one line of the other information section
   */
  virtual void ParseLine( string const & line ) override
  {
    m_comments += line;
  }

  /*!
   * @brief For this specific section there is no keyword
   */
  virtual void CheckKeywords() override
  {
    GEOSX_ERROR_IF( !m_lines.empty(),
                    "Invalid " << GetName() << " section. No keyword should be defined, only data." );
  }

  /*!
   * @brief Write the infotmation section to the \p file
   * @param[in] file the file to be output
   */
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
  /*!
   * @brief Constructor of the Log Data Section
   * @param[in] nbEntries the number of entries per logs
   * @param[in] nbCurves the number of data that are logged
   */
  LASASCIILogDataSection( localIndex nbEntries, localIndex nbCurves ):
    LASSection(),
    m_logs( nbCurves, nbEntries ),
    m_nbCurves( nbCurves ),
    m_nbLogEntries( nbEntries ),
    m_count( 0 )
  {}

  /*!
   * @brief Get the name of the section
   * @return a string containing the name of the section
   */
  virtual string const GetName() const override
  {
    return GetNameStatic();
  }

  /*!
   * @copydoc GetName()
   */
  static string GetNameStatic()
  {
    return "ASCII Log Data";
  }

  /*!
   * @brief Write the section to a file
   * @param[in] file the stream to the output file
   */
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
   * @param[in] i index of the log
   */
  arraySlice1d< real64 > GetLog( localIndex i ) const
  {
    return m_logs[i];
  }

private:
  /*!
   * @brief Parse the line an decompose it into a keyword, an unit, data and description
   */
  virtual void ParseLine( string const & line ) override;

private:
  /// Contains the well logs
  array2d< real64 > m_logs;

  /// The number of curves (i.e. the number of data)
  localIndex m_nbCurves;

  /// The number of log entries per curves.
  localIndex m_nbLogEntries;

  /// Use for internal counting while parsing the file
  localIndex m_count;
};

/*!
 * @brief Main class to read and write a LAS file
 */
class LASFile
{
public:
  /*!
   * @brief defaults constructor
   */
  LASFile() = default;

  /*!
   * @brief default destructor
   */
  ~LASFile(){}

  /*!
   * @brief Load a LAS file
   * @param[in] fileName path to the file to be loaded
   */
  void Load( string const & fileName );

  /*!
   * @brief Save a LAS file
   * @param[in] fileName path to the file to be saved
   */
  void Save( string const & fileName ) const;

  /*!
   * @brief Get a log giving its name
   * @param[in] name name of the log.
   * @return a slice containing the data
   */
  arraySlice1d< real64 > GetLog( string const & logName ) const;

  /*!
   * @brief loop through all the information sections
   * @param[in] lambda the lambda function to be applied during the loop
   */
  template< typename LAMBDA >
  void forInformationSections( LAMBDA && lambda )
  {
    for( auto & informationSection : m_lasInformationSections )
    {
      lambda( informationSection );
    }
  }

  /*!
   * @copydoc forInformationSections( LAMBDA && )
   */
  template< typename LAMBDA >
  void forInformationSections( LAMBDA && lambda ) const
  {
    for( auto & informationSection : m_lasInformationSections )
    {
      lambda( informationSection );
    }
  }

  /*!
   * @brief loop through all the log sections
   * @param[in] lambda the lambda function to be applied during the loop
   */
  template< typename LAMBDA >
  void forLogSections( LAMBDA && lambda )
  {
    for( auto & logSection : m_lasASCIILogDataSection )
    {
      lambda( logSection );
    }
  }

  /*!
   * @copydoc forLogSections( LAMBDA && )
   */
  template< typename LAMBDA >
  void forLogSections( LAMBDA && lambda ) const
  {
    for( auto & logSection : m_lasASCIILogDataSection )
    {
      lambda( logSection );
    }
  }

  /*!
   * @brief Get one information section
   * @param[in] i the index of the information section
   */
  LASInformationSection const & GetInformationSection( integer i ) const
  {
    return *m_lasInformationSections[i];
  }

  /*!
   * @brief Get one log section
   * @param[in] i the index of the log section
   */
  LASASCIILogDataSection const & GetLogSection( integer i ) const
  {
    return m_lasASCIILogDataSection[i];
  }

  /*!
   * @brief returns the number of log entries for a property
   * @param[in] logName the name of the property
   */
  localIndex LogSize( string const & logName ) const;

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
   * @param[i] keywork the keyword associated with the line(s) to be returned
   * @tparam Type of the section
   */
  template< typename T >
  std::vector< LASLine const * > const GetLASLines( string const & keyword ) const
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
  /*!
   * @brief Return the last section
   * @tparam Type of the section
   * @return the last section
   */
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
