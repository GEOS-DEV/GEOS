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

/**
 * @file VTKFile.hpp
 */

#ifndef VTKFILE_HPP_
#define VTKFILE_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/RestartFlags.hpp" 
#include "mesh/InterObjectRelation.hpp"
#include "codingUtilities/StringUtilities.hpp"

#include "pugixml.hpp"

#ifdef GEOSX_USE_MPI
#include <mpi.h>
#endif

namespace geosx
{

class DomainPartition;

class VTKFile
{
  public:
  VTKFile() = delete;

  /*!
   * @brief Initialize the VTK file writer
   * @details This constructor will construct the root file for
   * the output (.pvd)
   * @param[in] name the name of the pvd file
   */
  VTKFile( string const & name );

  /*!
   * @brief Set the plot level
   * @param[in] plotLevel the plot level. All fields flagged with
   * a plot level inferior or equal to this prescribed value will
   * be output
   */
  void SetPlotLevel( const int plotLevel )
  {
    m_plotLevel = dataRepository::IntToPlotLevel(plotLevel);
  }

  /*!
   * @brief Set the binary mode
   * @param[in] binary the binary mode
   */
  void SetBinaryMode( const bool binary )
  {
    m_binary = binary;
  }

  /*!
   * @brief Output a file for one time step
   * @param[in] cycle the cycle number
   */
  void Write( double const timeStep, 
              DomainPartition const & domain );

  private:
    /*!
     * @brief Create a XML Node for DataArray
     * @param[in,out] parent the parent XML node
     * @param[in] type a string containing the type of the field
     * @param[in] name the name of the field
     * @param[in] nbComponents dimension of the field
     * @param[in] p is a parallel data array
     * @return the corresponding xml node
     */
  pugi::xml_node CreateDataArray( pugi::xml_node & parent,
                                  string const & type,
                                  string const & name,
                                  int const & nbComponents,
                                  string const & format = "ascii",
                                  bool p = false)
  {
    pugi::xml_node dataArrayNode;
    if( p )
    {
      dataArrayNode = parent.append_child("PDataArray");
    }
    else
    {
      dataArrayNode = parent.append_child("DataArray");
    }
    dataArrayNode.append_attribute("type") = type.c_str();
    dataArrayNode.append_attribute("Name") = name.c_str();
    dataArrayNode.append_attribute("NumberOfComponents") = std::to_string(nbComponents).c_str();
    dataArrayNode.append_attribute("format") = format.c_str();
    return dataArrayNode;
  }

    /*!
     * @brief Create a XML Node for PDataArray
     * @param[in,out] parent the parent XML node
     * @param[in] type a string containing the type of the field
     * @param[in] name the name of the field
     * @param[in] nbComponents dimension of the field
     * @return the corresponding xml node
     */
  pugi::xml_node CreatePDataArray( pugi::xml_node & parent,
                                   string const & type,
                                   string const & name,
                                   int const & nbComponents,
                                   string const & format = "ascii" )
  {
    return CreateDataArray( parent, type, name, nbComponents, format, true);
  }

  private:
    class CustomVTUXMLWriter
    {
      public:
        CustomVTUXMLWriter() = delete;
        CustomVTUXMLWriter( string const & fileName ) :
          m_outFile( fileName, std::ios::binary ),
          m_spaceCount(0)
        {
        }

        void WriteHeader()
        {
          m_outFile << "<?xml version=\"1.0\"?>\n";
        }

        void OpenXMLNode( string const & nodeName, std::initializer_list< std::pair< string, string > > const & args)
        {
          for( int i = 0 ; i < m_spaceCount ; i++)
          {
            m_outFile << " ";
          }
          m_outFile << "<" << nodeName << " ";
          for( auto param  : args )
          {
            m_outFile << param.first << "=\"" << param.second << "\" ";
          }
          m_outFile << ">\n";
          m_spaceCount +=2;
        }

        void WriteVertices( r1_array const & vertices, bool binary )
        {
          if( binary )
          {
            WriteBinaryVertices( vertices );
          }
          else
          {
            WriteAsciiVertices( vertices );
          }
        }

        template< typename NODEMAPTYPE >
        void WriteCellConnectivities( string const& type, NODEMAPTYPE const & connectivities, bool binary )
        {
          if( binary )
          {
            WriteBinaryConnectivities( m_geosxToVTKCellTypeMap.at( type ), connectivities );
          }
          else
          {
            WriteAsciiConnectivities( m_geosxToVTKCellTypeMap.at( type ), connectivities );
          }
        }

        void WriteCellOffsets( localIndex nbNodesPerElement, localIndex nbElements, bool binary )
        {
          if( binary )
          {
            WriteBinaryOffsets( nbNodesPerElement, nbElements );
          }
          else
          {
            WriteAsciiOffsets( nbNodesPerElement, nbElements );
          }
        }


        void WriteCellTypes( string const& type, localIndex nbElements, bool binary )
        {
          if( binary )
          {
            WriteBinaryTypes( m_geosxToVTKCellTypeMap.at( type ), nbElements );
          }
          else
          {
            WriteAsciiTypes( m_geosxToVTKCellTypeMap.at( type ), nbElements );
          }
        }

        template< typename T >
        void WriteData(T const & data, bool binary)
        {
          if( binary )
          {
            WriteBinaryData( data );
          }
          else
          {
            WriteAsciiData( data );
          }
        }

        void CloseXMLNode( string const & nodeName )
        {
          m_spaceCount -=2;
          for( int i = 0 ; i < m_spaceCount ; i++)
          {
            m_outFile << " ";
          }
          m_outFile << "</" << nodeName << ">\n";
        }

      private:

        void WriteBinaryVertices( r1_array const & vertices )
        {
          std::stringstream stream;
          std::uint32_t size = integer_conversion< std::uint32_t >( vertices.size() ) * 3 *  sizeof( real64 );
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ));
          for( auto const & vertex : vertices )
          {
            stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( vertex.Data() ), 3*sizeof( real64 )) ;
          }
          DumpBuffer( stream );
        }

        void WriteAsciiVertices( r1_array const & vertices )
        {
          for( auto vertex : vertices )
          {
            m_outFile << vertex << "\n";
          }
        }

        template< typename NODEMAPTYPE >
        void WriteBinaryConnectivities( integer type, NODEMAPTYPE const & connectivities )
        {
          // TODO hardcoded, waiting for the VTK HEX node ordering in GEOSX
          std::stringstream stream;
          std::uint32_t size = integer_conversion< std::uint32_t > ( connectivities.size() ) * sizeof( localIndex );
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ));
          for( localIndex i = 0; i < connectivities.size() / 8 ; i++ )
          {
            localIndex_array connectivityIndex(8);
            connectivityIndex[0] = connectivities[i][0];
            connectivityIndex[1] = connectivities[i][1];
            connectivityIndex[2] = connectivities[i][3];
            connectivityIndex[3] = connectivities[i][2];
            connectivityIndex[4] = connectivities[i][4];
            connectivityIndex[5] = connectivities[i][5];
            connectivityIndex[6] = connectivities[i][7];
            connectivityIndex[7] = connectivities[i][6];
            stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivityIndex.data() ), sizeof( localIndex ) * 8 );
          }
          DumpBuffer( stream );
        }

        template< typename NODEMAPTYPE >
        void WriteAsciiConnectivities( integer type, NODEMAPTYPE const & connectivities )
        {
          // TODO hardcoded, waiting for the VTK HEX node ordering in GEOSX
          for( localIndex i = 0; i < connectivities.size() / 8 ; i++ )
          {
            m_outFile << connectivities[i][0] << " ";
            m_outFile << connectivities[i][1] << " ";
            m_outFile << connectivities[i][3] << " ";
            m_outFile << connectivities[i][2] << " ";
            m_outFile << connectivities[i][4] << " ";
            m_outFile << connectivities[i][5] << " ";
            m_outFile << connectivities[i][7] << " ";
            m_outFile << connectivities[i][6] << " ";
            m_outFile << "\n";
          }
        }

        void WriteAsciiOffsets( localIndex j, localIndex nb )
        {
          for( localIndex i = 1 ; i < nb + 1 ; i++)
          {
            m_outFile << i*j << "\n";
          }
        }

        void WriteBinaryOffsets( localIndex j, localIndex nb )
        {
          std::stringstream stream;
          std::uint32_t size = integer_conversion< std::uint32_t >( nb )* sizeof( localIndex );
          integer multiplier = FindMultiplier( sizeof( integer ) );
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ));
          localIndex_array offsetFragment( multiplier );
          for( integer i = 0; i < multiplier; i++)
          {
            offsetFragment[i] = ( i + 1 ) * j;
          }
          for( localIndex i = 0 ; i < nb / multiplier ; i++)
          {
            stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( offsetFragment.data() ), sizeof( localIndex ) * multiplier );
            for( integer k = 0; k < multiplier; k++)
            {
              offsetFragment[k] += j * multiplier;
            }
          }
          stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( offsetFragment.data() ), sizeof( localIndex ) * ( nb % multiplier) );
          DumpBuffer( stream );
        }

        void WriteAsciiTypes( integer type, localIndex nb )
        {
          for( localIndex i = 0 ; i < nb ; i++)
          {
            m_outFile << type << "\n";
          }
        }

        void WriteBinaryTypes( integer type, localIndex nb )
        {
          std::stringstream stream;
          std::uint32_t size = integer_conversion< std::uint32_t >( nb ) * sizeof( integer );
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ));
          integer multiplier = FindMultiplier( sizeof( integer ) );
          integer_array typeArray( multiplier );
          for( integer i = 0; i < multiplier; i++ )
          {
            typeArray[i] = type;
          }
          string typeString64 = stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( typeArray.data() ), sizeof( integer ) * multiplier );
          for( localIndex i = 0 ; i < nb / multiplier ; i++)
          {
            stream << typeString64;
          }
          stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( typeArray.data() ), sizeof( integer ) * ( nb % multiplier) );
          DumpBuffer( stream );
        }

        template< typename T >
        void WriteAsciiData( T const & data)
        {
          for( localIndex i = 0; i < data.size() ; i++ )
          {
            m_outFile << data[i] << "\n";
          }
        }

        template< typename T >
        void WriteBinaryData( T const & data )
        {
          std::stringstream stream;
          std::uint32_t size = integer_conversion < std::uint32_t >( data.size() ) * sizeof( real64 );
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ));
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( data.data() ), sizeof( data[0] ) * integer_conversion< integer >( data.size() ) );
          DumpBuffer( stream );
        }

        integer FindMultiplier( integer typeSize )
        {
          integer multiplier = 1;
          while( ( multiplier * typeSize) % 6 )
          {
            multiplier++;
          }
          return multiplier;
        }

      void DumpBuffer( std::stringstream const & stream )
      {
        m_outFile << stream.rdbuf() << '\n';
      }
      private:
        /// vtu output file
        std::ofstream m_outFile;

        /// Space counter to have well indented XML file
        int m_spaceCount;

        /// Map from GEOSX type to VTK cell types
        const unordered_map< string, int > m_geosxToVTKCellTypeMap =
        {
          { "C3D4", 10 },
          { "C3D5", 14 },
          { "C3D6", 13 },
          { "C3D8", 12 },
          { "", 9 } // QUAD ?
        };
    };
    const unordered_map< rtTypes::TypeIDs, string > m_geosxToVTKTypeMap =
    {
      {rtTypes::TypeIDs::integer_id, "Int32"},
      {rtTypes::TypeIDs::localIndex_id, "Int64"},
      {rtTypes::TypeIDs::globalIndex_id, "Int64"},
      {rtTypes::TypeIDs::real32_id, "Float32"},
      {rtTypes::TypeIDs::real64_id, "Float64"},
      {rtTypes::TypeIDs::r1_array_id, "Float64"},
      {rtTypes::TypeIDs::real64_array_id, "Float64"},
      {rtTypes::TypeIDs::real64_array2d_id, "Float64"},
      {rtTypes::TypeIDs::real64_array3d_id, "Float64"},
      {rtTypes::TypeIDs::real32_array_id, "Float32"},
      {rtTypes::TypeIDs::real32_array2d_id, "Float32"},
      {rtTypes::TypeIDs::real32_array3d_id, "Float32"},
      {rtTypes::TypeIDs::integer_array_id, "Int32"},
      {rtTypes::TypeIDs::localIndex_array_id, "Int64"},
      {rtTypes::TypeIDs::localIndex_array2d_id, "Int64"},
      {rtTypes::TypeIDs::localIndex_array3d_id, "Int64"},
      {rtTypes::TypeIDs::globalIndex_array_id, "Int64"},
      {rtTypes::TypeIDs::globalIndex_array2d_id, "Int64"},
      {rtTypes::TypeIDs::globalIndex_array3d_id, "Int64"},
    };
    /// Root file ( .pvd )
    pugi::xml_document m_rootFile;

    /// Unstructured file gathering all vtu files for a time step ( .pvtu )
    pugi::xml_document m_pvtuFile;

    /// Plot level
    dataRepository::PlotLevel m_plotLevel;

    /// Base name of the output
    string m_baseName;

    /// Tells wether or not the output is binary
    bool m_binary;
};

template<>
inline void VTKFile::CustomVTUXMLWriter::WriteBinaryData( r1_array const & data )
{
  std::stringstream stream;
  std::uint32_t size = integer_conversion< std::uint32_t > ( data.size() ) * sizeof( real64 ) * 3;
  stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ));
  for( auto const & elem : data )
  {
    stream << stringutilities::EncodeBase64(  reinterpret_cast< const unsigned char * >( elem.Data() ), sizeof( real64 ) * 3 );
  }
  m_outFile << stream.rdbuf() << '\n';
}
}
#endif /* VTKFILE_H_ */
