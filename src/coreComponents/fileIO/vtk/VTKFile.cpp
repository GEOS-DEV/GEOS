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
 * @file VTKFile.cpp
 */

#include "VTKFile.hpp"
#include <sys/stat.h>

#include "managers/DomainPartition.hpp"

#include "dataRepository/ViewWrapper.hpp"


namespace geosx
{
using namespace dataRepository;

namespace
{
/// Map from GEOSX type to VTK cell types
static std::unordered_map< string, int > geosxToVTKCellTypeMap =
{
  { "C3D4", 10 },
  { "C3D5", 14 },
  { "C3D6", 13 },
  { "C3D8", 12 },
  { "", 9 } // QUAD ?
};

static std::unordered_map< std::type_index, string > geosxToVTKTypeMap =
{
  {std::type_index( typeid( integer ) ), "Int32"},
  {std::type_index( typeid( localIndex ) ), "Int64"},
  {std::type_index( typeid( globalIndex ) ), "Int64"},
  {std::type_index( typeid( real32 ) ), "Float32"},
  {std::type_index( typeid( real64 ) ), "Float64"},
  {std::type_index( typeid( r1_array ) ), "Float64"},
  {std::type_index( typeid( real64_array ) ), "Float64"},
  {std::type_index( typeid( real64_array2d ) ), "Float64"},
  {std::type_index( typeid( real64_array3d ) ), "Float64"},
  {std::type_index( typeid( real32_array ) ), "Float32"},
  {std::type_index( typeid( real32_array2d ) ), "Float32"},
  {std::type_index( typeid( real32_array3d ) ), "Float32"},
  {std::type_index( typeid( integer_array ) ), "Int32"},
  {std::type_index( typeid( localIndex_array ) ), "Int64"},
  {std::type_index( typeid( localIndex_array2d ) ), "Int64"},
  {std::type_index( typeid( localIndex_array3d ) ), "Int64"},
  {std::type_index( typeid( globalIndex_array ) ), "Int64"},
  {std::type_index( typeid( globalIndex_array2d ) ), "Int64"},
  {std::type_index( typeid( globalIndex_array3d ) ), "Int64"}
};
}
class CustomVTUXMLWriter
{
  public:
  CustomVTUXMLWriter() = delete;
  CustomVTUXMLWriter( string const & fileName ) :
      m_outFile( fileName, std::ios::binary ),
      m_spaceCount(0)
  {
  }

  /*!
   * @brief Write the header of the VTU file with the XML version
   */
  void WriteHeader()
  {
    m_outFile << "<?xml version=\"1.0\"?>\n";
  }

  /*!
   * @brief Write an opening XML node
   * @param[in] nodeName the name of the node
   * @param[in[ args list of pair {paramaters,value}
   * @details This function also handle the spacing for good looking file
   */
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

  /*!
   * @brief Write the vertices coordinates
   * @param[in] vertices table of vertice coordinates
   * @param[in] binary tells wether or not the data should be written in binary format
   */
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

  /*!
   * @brief Write the cell connectivities
   * @param[in] type the GEOSX type of the cells
   * @param[in] connectivities the table of connectivities
   * @param[in] binary tells wether or not the data should be written in binary format
   */
  template< typename NODEMAPTYPE >
  void WriteCellConnectivities( string const& type, NODEMAPTYPE const & connectivities, bool binary )
  {
    if( binary )
    {
      WriteBinaryConnectivities( geosxToVTKCellTypeMap.at( type ), connectivities );
    }
    else
    {
      WriteAsciiConnectivities( geosxToVTKCellTypeMap.at( type ), connectivities );
    }
  }

  /*!
   * @brief Write the offsets
   * @details for a full hex mesh : 0, 8, 16, 24....
   * @param[in] nbNodesPerElement the number of nodes by elements
   * @param[in] nbElements the number of elements
   * @param[in] binary tells wether or not the data should be written in binary format
   */
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

  /*!
   * @brief Write the Cell types
   * @param[in] type the GEOSX type of the cells
   * @param[in] nbElements the number of elements
   * @param[in] binary tells wether or not the data should be written in binary format
   */
  void WriteCellTypes( string const& type, localIndex nbElements, bool binary )
  {
    if( binary )
    {
      WriteBinaryTypes( geosxToVTKCellTypeMap.at( type ), nbElements );
    }
    else
    {
      WriteAsciiTypes( geosxToVTKCellTypeMap.at( type ), nbElements );
    }
  }

  template< typename T >
  /*!
   * @brief Write contiguous data
   * @param[in] data table of data
   * @param[in] binary tells wether or not the data should be written in binary format
   */
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

  /*!
   * @brief Method use to write, cells data size
   * @param[in] nbTotalCells number of elements accross all the elementRegion
   * @param[in] factor usually the size of the data type to be written, multiplied by the number of data by cells.
   */
  void WriteSize( localIndex nbTotalCells, localIndex factor )
  {
    std::stringstream stream;
    std::uint32_t size = integer_conversion< std::uint32_t >( nbTotalCells * factor );
    stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ) );
    m_outFile << stream.rdbuf();
  }

  /*!
   * @brief Close a XML node
   * @param[in] nodeName name of the node to be closed
   */
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
      stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ) );
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
        std::stringstream stream;
        std::uint32_t size = integer_conversion< std::uint32_t > ( connectivities.size() ) * sizeof( localIndex );
        //stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), sizeof( std::uint32_t ) );
        if( type == 12 ) // Special case for the hex, because the ordering need to be changed for each cell...
        {
          integer multiplier = FindMultiplier( sizeof( localIndex ) );
          localIndex cellIndex = 0;
          localIndex vertexIndex=0;
          localIndex_array connectivityFragment( multiplier );
          std::cout << multiplier << std::endl;
          for( integer i = 0 ; i < connectivities.size() / multiplier ; i++ )
          {
            for( integer j = 0 ; j < multiplier; j++ )
            {
              if( vertexIndex == 2 )
              {
                connectivityFragment[j] = connectivities[cellIndex][3];
              }
              else if( vertexIndex == 3 )
              {
                connectivityFragment[j] = connectivities[cellIndex][2];
              }
              else if( vertexIndex == 6 )
              {
                connectivityFragment[j] = connectivities[cellIndex][7];
              }
              else if( vertexIndex == 7 )
              {
                connectivityFragment[j] = connectivities[cellIndex][6];
              }
              else
              {
                connectivityFragment[j] = connectivities[cellIndex][vertexIndex];
              }
              vertexIndex++;
              if( vertexIndex == 8)
              {
                vertexIndex = 0;
                cellIndex++;
              }
            }
            stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivityFragment.data() ), sizeof( localIndex ) * multiplier );
          }
          connectivityFragment.resize( 8 - vertexIndex );
          for( integer j = 0 ; j < connectivityFragment.size(); j++ )
          {
            if( vertexIndex == 2 )
            {
              connectivityFragment[j] = connectivities[cellIndex][3];
            }
            else if( vertexIndex == 3 )
            {
              connectivityFragment[j] = connectivities[cellIndex][2];
            }
            else if( vertexIndex == 6 )
            {
              connectivityFragment[j] = connectivities[cellIndex][7];
            }
            else if( vertexIndex == 7 )
            {
              connectivityFragment[j] = connectivities[cellIndex][6];
            }
            else
            {
              connectivityFragment[j] = connectivities[cellIndex][vertexIndex];
            }
            vertexIndex++;
            if( vertexIndex == 8)
            {
              vertexIndex = 0;
              cellIndex++;
            }
          }
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivityFragment.data() ), sizeof( localIndex ) * integer_conversion< integer >( connectivityFragment.size() ) );
        }
        else
        {
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivities.data() ), sizeof( localIndex ) * integer_conversion< integer >( connectivities.size() ) );
        }
        DumpBuffer( stream );
      }

    template< typename NODEMAPTYPE >
    void WriteAsciiConnectivities( integer type, NODEMAPTYPE const & connectivities )
    {
      if( type == 12 ) // Special case for hexahedron because of the internal ordering
      {
        for( localIndex i = 0; i < connectivities.size() / 8  ; i++ )
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
      else
      {
        for( localIndex i = 0; i < connectivities.size()  ; i++ )
        {
            m_outFile << connectivities.data()[i] <<" ";
        }
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
      integer multiplier = FindMultiplier( sizeof( integer ) ); // We do not write all the data at once to avoid creating a big table each time.
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
      integer multiplier = FindMultiplier( sizeof( integer ) );// We do not write all the data at once to avoid creating a big table each time.
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

    /*!
     * @brief This function is used to compute the minimum number of value of a certain type
     * that can be continously encoded into a base64 to be properly written into the VTU file.
     */
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
};

template<>
inline void CustomVTUXMLWriter::WriteBinaryData( r1_array const & data )
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

  VTKFile::VTKFile( string const & name ):
    m_baseName( name ),
    m_binary( false )
  {
    int mpiRank;
    MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
    if( mpiRank == 0 )
    {
      // Declaration of XML version
      auto declarationNode = m_rootFile.append_child(pugi::node_declaration);
      declarationNode.append_attribute("version") = "1.0";

      // Declaration of the node VTKFile
      auto vtkFileNode = m_rootFile.append_child("VTKFile");
      vtkFileNode.append_attribute("type") = "Collection";
      vtkFileNode.append_attribute("version") = "0.1";
      //vtkFileNode.append_attribute("byteOrder") = "LittleEndian";
      //vtkFileNode.append_attribute("compressor") = "vtkZLibDataCompressor";

      // Declaration of the node Collection
      vtkFileNode.append_child("Collection");
      mode_t mode = 0733;
      mkdir( name.c_str(), mode );

      string pvdFileName = name + ".pvd";
      m_rootFile.save_file(pvdFileName.c_str());
    }

  }

  void VTKFile::Write( double const timeStep,
                       DomainPartition const & domain )
  {
    int mpiRank;
    int mpiSize;
    MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
    MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
    ElementRegionManager const * elemManager = domain.getMeshBody(0)->getMeshLevel(0)->getElemManager();
    NodeManager const * nodeManager = domain.getMeshBody(0)->getMeshLevel(0)->getNodeManager();
    string timeStepFolderName = m_baseName + "/" + std::to_string( timeStep );
    string format;
    if( m_binary )
    {
      format = "binary";
    }  
    else
    {
      format = "ascii";
    }
    if( mpiRank == 0 )    
    {                     
      /// Add the new entry to the pvd root file
      auto collectionNode = m_rootFile.child("VTKFile").child("Collection");
      auto dataSetNode = collectionNode.append_child("DataSet");
      dataSetNode.append_attribute("timestep") = std::to_string( timeStep ).c_str();
      dataSetNode.append_attribute("group") = "";
      dataSetNode.append_attribute("part") = "0";
      string pvtuFileName = timeStepFolderName + "/root.pvtu";
      dataSetNode.append_attribute("file") = pvtuFileName.c_str();

      /// Create the pvtu file for this time step
      // Create a directory for this time step
      mode_t mode = 0733;
      mkdir( timeStepFolderName.c_str(), mode );
      pugi::xml_document pvtuFile;

      // Declaration of XML version
      auto declarationNode = pvtuFile.append_child(pugi::node_declaration);
      declarationNode.append_attribute("version") = "1.0";

      // Declaration of the node VTKFile
      auto vtkFileNode = pvtuFile.append_child("VTKFile");
      vtkFileNode.append_attribute("type") = "PUnstructuredGrid";
      vtkFileNode.append_attribute("version") = "0.1";
      if( m_binary )
      {
        vtkFileNode.append_attribute("byteOrder") = "LittleEndian";
      }

      // Declaration of the node PUnstructuredGrid
      auto pUnstructureGridNode = vtkFileNode.append_child("PUnstructuredGrid");
      pUnstructureGridNode.append_attribute("GhostLevel") = "1";

      // Declaration the node PPoints
      auto pPointsNode = pUnstructureGridNode.append_child("PPoints");
      // .... and the data array containg the positions
      CreatePDataArray( pPointsNode, geosxToVTKTypeMap.at( std::type_index( typeid( real64 ) ) ), "Position", 3, format );
      
      // Find all the node fields to output
      auto pPointDataNode = pUnstructureGridNode.append_child("PPointData");
      for( auto const & wrapperIter : nodeManager->wrappers() )
      {
        ViewWrapperBase const * const wrapper = wrapperIter.second;
        if( wrapper->getPlotLevel() < m_plotLevel )
        {
           string const fieldName = wrapper->getName();
           std::type_info const & typeID = wrapper->get_typeid();
           if( !geosxToVTKTypeMap.count( typeID ) )
             continue;
           int dimension = 0;
           rtTypes::TypeIDs fieldType = rtTypes::typeID(wrapper->get_typeid());
           if( fieldType == rtTypes::TypeIDs::r1_array_id )
           {
             dimension = 3;
           }
           else
           {
             dimension = 1;
           }
           CreatePDataArray(pPointDataNode, geosxToVTKTypeMap.at( typeID ), fieldName, dimension, format);
        }
      }

      // Declaration of the node PCells
      auto pCellsNode = pUnstructureGridNode.append_child("PCells");
      // .... and its data array defining the connectivities, types, and offsets
      CreatePDataArray( pCellsNode, geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ), "connectivity", 1, format ); //TODO harcoded for the moment
      CreatePDataArray( pCellsNode, geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ), "offsets", 1, format );
      CreatePDataArray( pCellsNode, geosxToVTKTypeMap.at( std::type_index( typeid( integer ) ) ), "types", 1, format );

      // Find all the cell fields to output
      auto pCellDataNode = pUnstructureGridNode.append_child("PCellData");
      elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                    auto const * const elemRegion )
      {
        elemRegion->forElementSubRegions([&]( auto const * const subRegion )
        {
          for( auto const & wrapperIter : subRegion->wrappers() )
          {
            ViewWrapperBase const * const wrapper = wrapperIter.second;

            if( wrapper->getPlotLevel() < m_plotLevel )
            {
              // the field name is the key to the map
              string const fieldName = wrapper->getName();
              std::type_info const & typeID = wrapper->get_typeid();
              rtTypes::TypeIDs fieldType = rtTypes::typeID(wrapper->get_typeid());
              if( !geosxToVTKTypeMap.count( typeID ) )
                continue;
              int dimension = 0;
              if( fieldType == rtTypes::TypeIDs::r1_array_id )
              {
                dimension = 3;
              }
              else
              {
                dimension = 1;
              }
              CreatePDataArray(pCellDataNode, geosxToVTKTypeMap.at( typeID ), fieldName, dimension, format);
            }
          }
       });
    });
    
    // Declaration of the "Piece" nodes refering to the vtu files
    for( int i = 0 ;  i < mpiSize ; i++ )
    {
      auto curPieceNode = pUnstructureGridNode.append_child("Piece");
      string fileName = std::to_string(i) + ".vtu";
      curPieceNode.append_attribute("Source") = fileName.c_str();
    }

    // Save the files
    string pvdFileName = m_baseName + ".pvd";
    m_rootFile.save_file(pvdFileName.c_str());
    pvtuFile.save_file(pvtuFileName.c_str());
  }
    
  string vtuFileName = timeStepFolderName + "/" + std::to_string(mpiRank) + ".vtu";
  CustomVTUXMLWriter vtuWriter( vtuFileName );
  vtuWriter.WriteHeader();
  vtuWriter.OpenXMLNode( "VTKFile", { {"type", "UnstructuredGrid"},
                                      {"version", "0.1"},
                                      {"byte_order", "LittleEndian"} } );
  vtuWriter.OpenXMLNode( "UnstructuredGrid",{} );
  
  // Declaration of the node Piece and the basic informations of the mesh
  localIndex totalNumberOfCells = 0;
  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    totalNumberOfCells += elemRegion->GetTotalSize();
  });
  vtuWriter.OpenXMLNode( "Piece", { { "NumberOfPoints", std::to_string(nodeManager->size() ) },
                                    { "NumberOfCells", std::to_string( totalNumberOfCells ) } } );

  // Definition of node Points
  vtuWriter.OpenXMLNode( "Points",{} );

  // Definition of the node DataArray that will contain all the node coordinates
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( real64 ) ) ) },
                                        { "Name", "Position" },
                                        { "NumberOfComponents", "3" },
                                        { "format", format } } );
  vtuWriter.WriteVertices( nodeManager->referencePosition(), m_binary );
  vtuWriter.CloseXMLNode( "DataArray" );
  vtuWriter.CloseXMLNode( "Points" );

  // Point data output
  vtuWriter.OpenXMLNode( "PointData", {} );
  for( auto const & wrapperIter : nodeManager->wrappers() )
  {
    ViewWrapperBase const * const wrapper = wrapperIter.second;
    if( wrapper->getPlotLevel() < m_plotLevel )
    {
       string const fieldName = wrapper->getName();
       std::type_info const & typeID = wrapper->get_typeid();
       rtTypes::TypeIDs fieldType = rtTypes::typeID(wrapper->get_typeid());
       if( !geosxToVTKTypeMap.count( typeID ) )
         continue;
       int dimension = 0;
       if( fieldType == rtTypes::TypeIDs::r1_array_id )
       {
         dimension = 3;
       }
       else
       {
         dimension = 1;
       }
       vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( typeID ) },
                                             { "Name", fieldName },
                                             { "NumberOfComponents", std::to_string( dimension ) },
                                             { "format", format } } );
       std::type_index typeIndex = std::type_index( typeID );
       rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ),
                                       [&]( auto type ) -> void
       {
         using cType = decltype(type);
         const ViewWrapper< cType > & view = ViewWrapper<cType>::cast( *wrapper );
         vtuWriter.WriteData( view.reference(), m_binary );
       });
       vtuWriter.CloseXMLNode( "DataArray" );
    }
  }
  vtuWriter.CloseXMLNode( "PointData" );

  // Definition of the node Cells
  vtuWriter.OpenXMLNode( "Cells", {} );

  // Definition of the node DataArray that will contain the connectivities
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ) },
                                        { "Name", "connectivity" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  if( m_binary )
  {
    localIndex totalNumberOfConnectivities = 0;
    elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                  auto const * const elemRegion )
    {
      elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
      {
        totalNumberOfConnectivities = elemSubRegion->size() * elemSubRegion->numNodesPerElement();
      });
    });
    vtuWriter.WriteSize( totalNumberOfConnectivities, sizeof( localIndex ) );
  }

  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellConnectivities( elemSubRegion->GetElementTypeString(), elemSubRegion->nodeList(), m_binary );
    });
  });
  vtuWriter.CloseXMLNode( "DataArray" );

  // Definition of the node DataArray that will contain the offsets
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ) },
                                        { "Name", "offsets" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  vtuWriter.WriteSize( totalNumberOfCells, sizeof( localIndex ) );
  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellOffsets( elemSubRegion->numNodesPerElement(), elemSubRegion->size(), m_binary );
    });
  });
  vtuWriter.CloseXMLNode( "DataArray" );

  // Definition of the node DataArray that will contain the cell types
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( integer ) ) ) },
                                        { "Name", "types" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  vtuWriter.WriteSize( totalNumberOfCells, sizeof( integer ) );
  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellTypes( elemSubRegion->GetElementTypeString(), elemSubRegion->size(), m_binary );
    });
  });

  vtuWriter.CloseXMLNode( "DataArray" );
  vtuWriter.CloseXMLNode( "Cells" );

  // Definition of the CellDataArray node that will contains all the data held by the elements
  vtuWriter.OpenXMLNode( "CellData", {} );
  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    elemRegion->forElementSubRegions([&]( auto const * const subRegion )
    {
      for( auto const & wrapperIter : subRegion->wrappers() )
      {
        ViewWrapperBase const * const wrapper = wrapperIter.second;

        if( wrapper->getPlotLevel() < m_plotLevel )
        {
          string const fieldName = wrapper->getName();
          std::type_info const & typeID = wrapper->get_typeid();
          rtTypes::TypeIDs fieldType = rtTypes::typeID(wrapper->get_typeid());
          if( !geosxToVTKTypeMap.count( typeID ) )
            continue;
          int dimension = 0;
          if( fieldType == rtTypes::TypeIDs::r1_array_id )
          {
            dimension = 3;
          }
          else
          {
            dimension = 1;
          }
          vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( typeID ) },
                                                { "Name", fieldName },
                                                { "NumberOfComponents", std::to_string( dimension ) },
                                                { "format", format } } );
          std::type_index typeIndex = std::type_index( typeID );
          rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ),
                                          [&]( auto type ) -> void
          {
            using cType = decltype(type);
            const ViewWrapper< cType > & view = ViewWrapper<cType>::cast( *wrapper );
            vtuWriter.WriteData( view.reference(), m_binary );
          });
          vtuWriter.CloseXMLNode( "DataArray" );
        }
      }
    });
  });
  vtuWriter.CloseXMLNode( "CellData" );
  vtuWriter.CloseXMLNode( "Piece" );
  vtuWriter.CloseXMLNode( "UnstructuredGrid" );
  vtuWriter.CloseXMLNode( "VTKFile" );
}

}
