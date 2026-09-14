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
 * @file WellPropWriter.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPROPWRITER_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPROPWRITER_HPP

#include <fstream>
#include <functional>
#include <utility>
#include <vector>

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementSubRegion.hpp"

namespace geos
{

/******************************** WellPropWriter ********************************/

class WellPropWriter
{
public:
  using ColumnNames = std::vector< string >;

  struct PerforationContext
  {
    localIndex perforationIndex;
    localIndex reservoirRegion;
    localIndex reservoirSubRegion;
    localIndex reservoirElementIndex;
  };

  using SegmentWriter = std::function< void( localIndex, std::ofstream & ) >;
  using PerforationWriter = std::function< void( PerforationContext const &, std::ofstream & ) >;

  WellPropWriter() = default;
  WellPropWriter( WellPropWriter && ) = default;
  WellPropWriter( WellPropWriter const & ) = delete;
  WellPropWriter & operator=( WellPropWriter && ) = default;
  WellPropWriter & operator=( WellPropWriter const & ) = delete;

  ~WellPropWriter()
  {
    closeFile( m_segmentOutputFile );
    closeFile( m_perforationOutputFile );
  }

  void initialize( int const rank,
                   string const & outputDir,
                   string const & wellName,
                   string_array phaseNames,
                   string_array componentNames,
                   WellElementSubRegion & subRegion,
                   PerforationData const & perforationData )
  {
    GEOS_THROW_IF( m_isInitialized, "WellPropWriter may only be initialized once per well subregion.", geos::RuntimeError );

    openFile( m_segmentOutputFile, outputDir, wellName + "_seg_" + std::to_string( rank ) + ".csv" );
    openFile( m_perforationOutputFile, outputDir, wellName + "_perf_" + std::to_string( rank ) + ".csv" );

    m_phaseNames = std::move( phaseNames );
    m_componentNames = std::move( componentNames );
    m_numPhase = m_phaseNames.size();
    m_numComponent = m_componentNames.size();

    m_numSegments = subRegion.size();
    m_elementGhostRank = subRegion.ghostRank();
    m_globalWellElementIndex = subRegion.getGlobalWellElementIndex();

    m_numPerforations = perforationData.size();
    m_reservoirElementRegion = perforationData.getField< fields::perforation::reservoirElementRegion >();
    m_reservoirElementSubRegion = perforationData.getField< fields::perforation::reservoirElementSubRegion >();
    m_reservoirElementIndex = perforationData.getField< fields::perforation::reservoirElementIndex >();
    m_perforationWellElementIndex = perforationData.getField< fields::perforation::wellElementIndex >();
    m_perforationReservoirElementGlobalIndex = perforationData.getField< fields::perforation::reservoirElementGlobalIndex >();

    m_isInitialized = true;
  }

  bool isInitialized() const
  {
    return m_isInitialized;
  }

  template< typename T >
  void registerSegmentScalar( string const & name, T const & values )
  {
    addSegmentColumn( name,
                      [values]( localIndex const segmentIndex, std::ofstream & stream )
    {
      stream << "," << values[segmentIndex];
    } );
  }

  template< typename T >
  void registerSegmentConstitutiveScalar( string const & name, T const & values )
  {
    addSegmentColumn( name,
                      [values]( localIndex const segmentIndex, std::ofstream & stream )
    {
      stream << "," << values[segmentIndex][0];
    } );
  }

  template< typename T >
  void registerSegmentNamedColumns( ColumnNames const & names, T const & values )
  {
    addSegmentColumns(
      names,
      [values, columnCount = integer( names.size() )]( localIndex const segmentIndex, std::ofstream & stream )
      {
        writeColumns( stream, columnCount,
                      [&]( integer const columnIndex )
        {
          return values[segmentIndex][columnIndex];
        } );
      } );
  }

  template< typename T >
  void registerSegmentComponentColumns( string const & name, T const & values )
  {
    addSegmentColumns(
      componentColumnNames( name ),
      [values, componentCount = m_numComponent]( localIndex const segmentIndex, std::ofstream & stream )
      {
        writeColumns( stream, componentCount,
                      [&]( integer const componentIndex )
        {
          return values[segmentIndex][componentIndex];
        } );
      } );
  }

  template< typename T >
  void registerSegmentPhaseColumns( string const & name, T const & values )
  {
    addSegmentColumns(
      phaseColumnNames( name ),
      [values, phaseCount = m_numPhase]( localIndex const segmentIndex, std::ofstream & stream )
      {
        writeColumns( stream, phaseCount,
                      [&]( integer const phaseIndex )
        {
          return values[segmentIndex][phaseIndex];
        } );
      } );
  }

  template< typename T >
  void registerSegmentConstitutivePhaseColumns( string const & name, T const & values )
  {
    addSegmentColumns(
      phaseColumnNames( name ),
      [values, phaseCount = m_numPhase]( localIndex const segmentIndex, std::ofstream & stream )
      {
        writeColumns( stream, phaseCount,
                      [&]( integer const phaseIndex )
        {
          return values[segmentIndex][0][phaseIndex];
        } );
      } );
  }

  template< typename T >
  void registerSegmentPhaseDerivativeColumns( string const & name,
                                              ColumnNames const & derivativeNames,
                                              T const & values )
  {
    addSegmentColumns(
      crossColumnNames( name, m_phaseNames, derivativeNames ),
      [values, phaseCount = m_numPhase, derivativeCount = integer( derivativeNames.size() )](
        localIndex const segmentIndex,
        std::ofstream & stream )
      {
        writeColumnGrid( stream, phaseCount, derivativeCount,
                         [&]( integer const phaseIndex, integer const derivativeIndex )
        {
          return values[segmentIndex][phaseIndex][derivativeIndex];
        } );
      } );
  }

  template< typename T >
  void registerSegmentConstitutivePhaseDerivativeColumns( string const & name,
                                                          ColumnNames const & derivativeNames,
                                                          T const & values )
  {
    addSegmentColumns(
      crossColumnNames( name, m_phaseNames, derivativeNames ),
      [values, phaseCount = m_numPhase, derivativeCount = integer( derivativeNames.size() )](
        localIndex const segmentIndex,
        std::ofstream & stream )
      {
        writeColumnGrid( stream, phaseCount, derivativeCount,
                         [&]( integer const phaseIndex, integer const derivativeIndex )
        {
          return values[segmentIndex][0][phaseIndex][derivativeIndex];
        } );
      } );
  }

  template< typename T >
  void registerSegmentConstitutivePhaseComponentColumns( string const & name, T const & values )
  {
    addSegmentColumns(
      crossColumnNames( name, m_phaseNames, m_componentNames ),
      [values, phaseCount = m_numPhase, componentCount = m_numComponent](
        localIndex const segmentIndex,
        std::ofstream & stream )
      {
        writeColumnGrid( stream, phaseCount, componentCount,
                         [&]( integer const phaseIndex, integer const componentIndex )
        {
          return values[segmentIndex][0][phaseIndex][componentIndex];
        } );
      } );
  }

  template< typename T >
  void registerPerforationScalar( string const & name, T const & values )
  {
    addPerforationColumn( name,
                          [values]( PerforationContext const & context, std::ofstream & stream )
    {
      stream << "," << values[context.perforationIndex];
    } );
  }

  template< typename T >
  void registerPerforationNamedColumns( ColumnNames const & names, T const & values )
  {
    addPerforationColumns(
      names,
      [values, columnCount = integer( names.size() )]( PerforationContext const & context, std::ofstream & stream )
      {
        writeColumns( stream, columnCount,
                      [&]( integer const columnIndex )
        {
          return values[context.perforationIndex][columnIndex];
        } );
      } );
  }

  template< typename T >
  void registerPerforationComponentColumns( string const & name, T const & values )
  {
    addPerforationColumns(
      componentColumnNames( name ),
      [values, componentCount = m_numComponent]( PerforationContext const & context, std::ofstream & stream )
      {
        writeColumns( stream, componentCount,
                      [&]( integer const componentIndex )
        {
          return values[context.perforationIndex][componentIndex];
        } );
      } );
  }

  template< typename T >
  void registerReservoirScalarAtPerforation( string const & name, T const & values )
  {
    addPerforationColumn( name,
                          [values]( PerforationContext const & context, std::ofstream & stream )
    {
      stream << "," << values[context.reservoirRegion][context.reservoirSubRegion][context.reservoirElementIndex];
    } );
  }

  template< typename T >
  void registerReservoirConstitutiveScalarAtPerforation( string const & name, T const & values )
  {
    addPerforationColumn( name,
                          [values]( PerforationContext const & context, std::ofstream & stream )
    {
      stream << "," << values[context.reservoirRegion][context.reservoirSubRegion][0][context.reservoirElementIndex];
    } );
  }

  template< typename T >
  void registerReservoirPhaseColumnsAtPerforation( string const & name, T const & values )
  {
    addPerforationColumns(
      phaseColumnNames( name ),
      [values, phaseCount = m_numPhase]( PerforationContext const & context, std::ofstream & stream )
      {
        writeColumns( stream, phaseCount,
                      [&]( integer const phaseIndex )
        {
          return values[context.reservoirRegion][context.reservoirSubRegion][context.reservoirElementIndex][phaseIndex];
        } );
      } );
  }

  template< typename T >
  void registerReservoirConstitutivePhaseColumnsAtPerforation( string const & name, T const & values )
  {
    addPerforationColumns(
      phaseColumnNames( name ),
      [values, phaseCount = m_numPhase]( PerforationContext const & context, std::ofstream & stream )
      {
        writeColumns( stream, phaseCount,
                      [&]( integer const phaseIndex )
        {
          return values[context.reservoirRegion][context.reservoirSubRegion][context.reservoirElementIndex][0][phaseIndex];
        } );
      } );
  }

  template< typename T >
  void registerReservoirConstitutivePhaseComponentColumnsAtPerforation( string const & name, T const & values )
  {
    addPerforationColumns(
      crossColumnNames( name, m_phaseNames, m_componentNames ),
      [values, phaseCount = m_numPhase, componentCount = m_numComponent](
        PerforationContext const & context,
        std::ofstream & stream )
      {
        writeColumnGrid( stream, phaseCount, componentCount,
                         [&]( integer const phaseIndex, integer const componentIndex )
        {
          return values[context.reservoirRegion][context.reservoirSubRegion][context.reservoirElementIndex][0][phaseIndex][componentIndex];
        } );
      } );
  }

  void writeTimeStep( real64 const time,
                      real64 const dt,
                      integer const cycle,
                      integer const subevent,
                      integer const timeStep,
                      integer const newtonIteration,
                      integer const numTimeStepCuts )
  {
    ensureInitialized();
    writeHeadersIfNeeded();
    validateCurrentRegistration();

    for( localIndex segmentIndex = 0; segmentIndex < m_numSegments; ++segmentIndex )
    {
      if( m_elementGhostRank[segmentIndex] < 0 )
      {
        m_segmentOutputFile << time << "," << dt << "," << cycle << "," << subevent << "," << timeStep << ","
                            << newtonIteration << "," << numTimeStepCuts << "," << m_globalWellElementIndex[segmentIndex];

        for( SegmentWriter const & writer : m_segmentWriters )
        {
          writer( segmentIndex, m_segmentOutputFile );
        }
        m_segmentOutputFile << std::endl;
      }
    }

    for( localIndex perforationIndex = 0; perforationIndex < m_numPerforations; ++perforationIndex )
    {
      localIndex const wellElementIndex = m_perforationWellElementIndex[perforationIndex];
      if( m_elementGhostRank[wellElementIndex] < 0 )
      {
        PerforationContext const context {
          perforationIndex,
          m_reservoirElementRegion[perforationIndex],
          m_reservoirElementSubRegion[perforationIndex],
          m_reservoirElementIndex[perforationIndex]
        };

        m_perforationOutputFile << time << "," << dt << "," << cycle << "," << subevent << "," << timeStep << ","
                                << newtonIteration << "," << numTimeStepCuts << ","
                                << m_perforationReservoirElementGlobalIndex[perforationIndex] << ","
                                << m_globalWellElementIndex[wellElementIndex];

        for( PerforationWriter const & writer : m_perforationWriters )
        {
          writer( context, m_perforationOutputFile );
        }
        m_perforationOutputFile << std::endl;
      }
    }

    m_segmentWriters.clear();
    m_perforationWriters.clear();
    m_numRegisteredSegmentColumns = 0;
    m_numRegisteredPerforationColumns = 0;
  }

private:
  template< typename Writer >
  void addSegmentColumn( string const & name, Writer && writer )
  {
    addSegmentColumns( ColumnNames { name }, std::forward< Writer >( writer ) );
  }

  template< typename Writer >
  void addSegmentColumns( ColumnNames const & names, Writer && writer )
  {
    ensureInitialized();
    trackSegmentColumns( names );
    m_segmentWriters.emplace_back( std::forward< Writer >( writer ) );
  }

  template< typename Writer >
  void addPerforationColumn( string const & name, Writer && writer )
  {
    addPerforationColumns( ColumnNames { name }, std::forward< Writer >( writer ) );
  }

  template< typename Writer >
  void addPerforationColumns( ColumnNames const & names, Writer && writer )
  {
    ensureInitialized();
    trackPerforationColumns( names );
    m_perforationWriters.emplace_back( std::forward< Writer >( writer ) );
  }

  void ensureInitialized() const
  {
    GEOS_THROW_IF( !m_isInitialized, "WellPropWriter must be initialized before registering output columns.", geos::RuntimeError );
  }

  static void closeFile( std::ofstream & stream )
  {
    if( stream.is_open() )
    {
      stream.close();
    }
  }

  static void openFile( std::ofstream & stream,
                        string const & outputDir,
                        string const & fileName )
  {
    closeFile( stream );
    stream.open( outputDir + "/" + fileName );
  }

  template< typename Range >
  static void appendColumns( ColumnNames & columns, Range const & names )
  {
    for( string const & name : names )
    {
      columns.push_back( name );
    }
  }

  void trackSegmentColumns( ColumnNames const & names )
  {
    trackColumns( names, m_segmentHeader, m_numRegisteredSegmentColumns, "segment" );
  }

  void trackPerforationColumns( ColumnNames const & names )
  {
    trackColumns( names, m_perforationHeader, m_numRegisteredPerforationColumns, "perforation" );
  }

  void trackColumns( ColumnNames const & names,
                     ColumnNames & existingHeader,
                     localIndex & registrationOffset,
                     char const * const writerKind )
  {
    localIndex const numNames = localIndex( names.size() );
    if( !m_headersWritten )
    {
      appendColumns( existingHeader, names );
    }
    else
    {
      GEOS_THROW_IF( registrationOffset + numNames > localIndex( existingHeader.size() ),
                     GEOS_FMT( "WellPropWriter {} registration exceeded the number of header columns.", writerKind ),
                     geos::RuntimeError );

      for( localIndex i = 0; i < numNames; ++i )
      {
        GEOS_THROW_IF( existingHeader[registrationOffset + i] != names[i],
                       GEOS_FMT( "WellPropWriter {} registration mismatch at column {}. Expected '{}', received '{}'.",
                                 writerKind,
                                 registrationOffset + i,
                                 existingHeader[registrationOffset + i],
                                 names[i] ),
                       geos::RuntimeError );
      }
    }

    registrationOffset += numNames;
  }

  template< typename Accessor >
  static void writeColumns( std::ofstream & stream,
                            integer const columnCount,
                            Accessor && accessor )
  {
    for( integer columnIndex = 0; columnIndex < columnCount; ++columnIndex )
    {
      stream << "," << accessor( columnIndex );
    }
  }

  template< typename Accessor >
  static void writeColumnGrid( std::ofstream & stream,
                               integer const outerCount,
                               integer const innerCount,
                               Accessor && accessor )
  {
    for( integer outerIndex = 0; outerIndex < outerCount; ++outerIndex )
    {
      for( integer innerIndex = 0; innerIndex < innerCount; ++innerIndex )
      {
        stream << "," << accessor( outerIndex, innerIndex );
      }
    }
  }

  template< typename Range >
  static ColumnNames suffixedColumnNames( string const & prefix, Range const & suffixes )
  {
    ColumnNames columnNames;
    for( string const & suffix : suffixes )
    {
      columnNames.push_back( prefix + "_" + suffix );
    }
    return columnNames;
  }

  template< typename OuterRange, typename InnerRange >
  static ColumnNames crossColumnNames( string const & prefix,
                                       OuterRange const & outerNames,
                                       InnerRange const & innerNames )
  {
    ColumnNames columnNames;
    for( string const & outerName : outerNames )
    {
      for( string const & innerName : innerNames )
      {
        columnNames.push_back( prefix + "_" + outerName + "_" + innerName );
      }
    }
    return columnNames;
  }

  ColumnNames componentColumnNames( string const & prefix ) const
  {
    return suffixedColumnNames( prefix, m_componentNames );
  }

  ColumnNames phaseColumnNames( string const & prefix ) const
  {
    return suffixedColumnNames( prefix, m_phaseNames );
  }

  void writeHeadersIfNeeded()
  {
    if( m_headersWritten )
    {
      return;
    }

    m_segmentOutputFile << "Time,Dt,Cycle,SubEvent,TimeStep,NewtonIteration,TimeStepCuts,Element";
    for( string const & column : m_segmentHeader )
    {
      m_segmentOutputFile << "," << column;
    }
    m_segmentOutputFile << std::endl;

    m_perforationOutputFile << "Time,Dt,Cycle,SubEvent,TimeStep,NewtonIteration,TimeStepCuts,ResElement,WellElement";
    for( string const & column : m_perforationHeader )
    {
      m_perforationOutputFile << "," << column;
    }
    m_perforationOutputFile << std::endl;

    m_headersWritten = true;
  }

  void validateCurrentRegistration() const
  {
    localIndex const numSegmentHeaderColumns = localIndex( m_segmentHeader.size() );
    localIndex const numPerforationHeaderColumns = localIndex( m_perforationHeader.size() );

    GEOS_THROW_IF( m_numRegisteredSegmentColumns != numSegmentHeaderColumns,
                   GEOS_FMT( "WellPropWriter segment registration is incomplete: expected {} columns, received {}.",
                             numSegmentHeaderColumns,
                             m_numRegisteredSegmentColumns ),
                   geos::RuntimeError );

    GEOS_THROW_IF( m_numRegisteredPerforationColumns != numPerforationHeaderColumns,
                   GEOS_FMT( "WellPropWriter perforation registration is incomplete: expected {} columns, received {}.",
                             numPerforationHeaderColumns,
                             m_numRegisteredPerforationColumns ),
                   geos::RuntimeError );
  }

  bool m_isInitialized = false;
  bool m_headersWritten = false;
  localIndex m_numRegisteredSegmentColumns = 0;
  localIndex m_numRegisteredPerforationColumns = 0;

  string_array m_phaseNames;
  string_array m_componentNames;
  localIndex m_numSegments = 0;
  integer m_numComponent = 0;
  integer m_numPhase = 0;

  ColumnNames m_segmentHeader;
  std::ofstream m_segmentOutputFile;
  std::vector< SegmentWriter > m_segmentWriters;

  localIndex m_numPerforations = 0;
  arrayView1d< localIndex const > m_reservoirElementRegion;
  arrayView1d< localIndex const > m_reservoirElementSubRegion;
  arrayView1d< localIndex const > m_reservoirElementIndex;
  arrayView1d< localIndex const > m_perforationWellElementIndex;
  arrayView1d< globalIndex const > m_perforationReservoirElementGlobalIndex;

  ColumnNames m_perforationHeader;
  std::ofstream m_perforationOutputFile;
  std::vector< PerforationWriter > m_perforationWriters;

  /// Global index of local element
  arrayView1d< globalIndex const > m_globalWellElementIndex;
  arrayView1d< integer const > m_elementGhostRank;
};

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLPROPWRITER_HPP
