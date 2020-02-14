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

/**
 * @file VTKFile.cpp
 */

// Source includes
#include "VTKFile.hpp"
#include "dataRepository/Wrapper.hpp"
#include "managers/DomainPartition.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

// System includes
#include <sys/stat.h>


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
  void WriteVertices( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & vertices, bool binary )
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
   */
  void WriteCellConnectivities( ElementRegionManager const * const elemManager, bool binary )
  {
    if( binary )
    {
      localIndex totalNumberOfConnectivities = 0;
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                        auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          totalNumberOfConnectivities += elemSubRegion->size() * elemSubRegion->numNodesPerElement();
        });
      });
      WriteSize( totalNumberOfConnectivities, sizeof( localIndex ) );
      WriteBinaryConnectivities( elemManager );
    }
    else
    {
      WriteAsciiConnectivities( elemManager );
    }
  }

  /*!
   * @brief Write the offsets
   * @details for a full hex mesh : 0, 8, 16, 24....
   */
  void WriteCellOffsets( ElementRegionManager const * const elemManager, bool binary )
  {
    if( binary )
    {
      WriteSize( elemManager->getNumberOfElements< CellElementSubRegion >(), sizeof( localIndex ) );
      WriteBinaryOffsets( elemManager );
    }
    else
    {
      WriteAsciiOffsets( elemManager );
    }
  }

  /*!
   * @brief Write the Cell types
   */
  void WriteCellTypes( ElementRegionManager const * const elemManager, bool binary )
  {
    if( binary )
    {
      WriteSize( elemManager->getNumberOfElements< CellElementSubRegion >(), sizeof( integer ) );
      WriteBinaryTypes( elemManager );
    }
    else
    {
      WriteAsciiTypes( elemManager );
    }
  }

  template< typename T >
  void WriteCellData(  ElementRegionManager::ElementViewAccessor< T > const & dataView, ElementRegionManager const * const elemManager, bool binary)
  {
    if( binary )
    {
      WriteCellBinaryData( dataView, elemManager );
    }
    else
    {
      WriteCellAsciiData( dataView, elemManager );
    }
  }

  template< typename T >
  void WriteNodeData(  Wrapper< T > const & dataView, bool binary)
  {
    if( binary )
    {
      WriteNodeBinaryData( dataView ); 
    }
    else
    {
      WriteNodeAsciiData( dataView );
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
    string outputString;
    outputString.resize( FindBase64StringLength( sizeof(std::uint32_t ) ) );
    stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), outputString, sizeof( std::uint32_t ) );
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

    void WriteBinaryVertices( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & vertices )
    {
      std::stringstream stream;
      std::uint32_t size = integer_conversion< std::uint32_t >( vertices.size() ) * sizeof( real64 );
      string outputString;
      outputString.resize( FindBase64StringLength( sizeof(std::uint32_t ) ) );
      stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &size ), outputString, sizeof( std::uint32_t ) );
      outputString.resize(FindBase64StringLength( sizeof( real64 ) * 3 ) );
      for ( localIndex i = 0; i < vertices.size( 0 ); ++i )
      {
        for ( localIndex j = 0; j < 3; ++j )
        {
          stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( &vertices( i, j ) ), outputString, sizeof( real64 ) );
        }
      }
      DumpBuffer( stream );
    }

    void WriteAsciiVertices( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & vertices )
    {
      for ( localIndex i = 0; i < vertices.size( 0 ); ++i )
      {
        m_outFile << vertices( i, 0 ) << " " << vertices( i, 1 ) << "" << vertices( i, 2 ) << "\n";
      }
    }

    void WriteBinaryConnectivities( ElementRegionManager const * const elemManager )
    {
      std::stringstream stream;
      integer multiplier = FindMultiplier( sizeof( localIndex ) );
      string outputString;
      outputString.resize( FindBase64StringLength( multiplier * sizeof( localIndex) ) );
      localIndex_array connectivityFragment( multiplier );
      integer countConnectivityFragment = 0;
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                               auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          integer type = geosxToVTKCellTypeMap.at( elemSubRegion->GetElementTypeString() );
          auto & connectivities = elemSubRegion->nodeList();
          if( type == 12 ) // Special case for hexahedron because of the internal ordering
          {
            for( localIndex i = 0 ; i < elemSubRegion->size() ; i++ )
            {
              for( integer j = 0; j < elemSubRegion->numNodesPerElement(); j++ )
              {
                if( j == 2 )
                {
                  connectivityFragment[countConnectivityFragment++] = connectivities[i][3];
                }
                else if( j == 3 )
                {
                  connectivityFragment[countConnectivityFragment++] = connectivities[i][2];
                }
                else if( j == 6 )
                {
                  connectivityFragment[countConnectivityFragment++] = connectivities[i][7];
                }
                else if( j == 7 )
                {
                  connectivityFragment[countConnectivityFragment++] = connectivities[i][6];
                }
                else
                {
                  connectivityFragment[countConnectivityFragment++] = connectivities[i][j];
                }
                if( countConnectivityFragment == multiplier )
                {
                  stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivityFragment.data() ), outputString, sizeof( localIndex ) * multiplier );
                  countConnectivityFragment = 0;
                }
              }
            }
          }
          else
          {
            for( localIndex i = 0 ; i < elemSubRegion->size() ; i++ )
            {
              for( integer j = 0; j < elemSubRegion->numNodesPerElement(); j++ )
              {
                connectivityFragment[countConnectivityFragment++] = connectivities[i][j];
                if( countConnectivityFragment == multiplier )
                {
                  stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivityFragment.data() ), outputString, sizeof( localIndex ) * multiplier );
                  countConnectivityFragment = 0;
                }
              }
            }
          }
        });
      });
      outputString.resize( FindBase64StringLength( sizeof( localIndex ) * ( countConnectivityFragment) ) );
      stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( connectivityFragment.data() ), outputString, sizeof( localIndex ) * ( countConnectivityFragment) );
      DumpBuffer( stream );
    }

    void WriteAsciiConnectivities( ElementRegionManager const * const elemManager )
    {
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          integer type = geosxToVTKCellTypeMap.at( elemSubRegion->GetElementTypeString() );
          auto & connectivities = elemSubRegion->nodeList();
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
        });
      });
    }

    void WriteAsciiOffsets( ElementRegionManager const * const elemManager )
    {
      localIndex curOffset = elemManager->GetRegion(0)->GetSubRegion(0)->numNodesPerElement();
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          localIndex offSetForOneCell = elemSubRegion->numNodesPerElement();
          for( localIndex i =  0; i < elemSubRegion->size(); i++ )
          {
            m_outFile << curOffset << "\n";
            curOffset += offSetForOneCell;
          }
        });
      });
    }

    void WriteBinaryOffsets( ElementRegionManager const * const elemManager )
    {
      std::stringstream stream;
      integer multiplier = FindMultiplier( sizeof( integer ) ); // We do not write all the data at once to avoid creating a big table each time.
      localIndex_array offsetFragment( multiplier );
      string outputString;
      outputString.resize( FindBase64StringLength( sizeof( localIndex ) * multiplier ) );
      integer countOffsetFragmentIndex = 0;
      localIndex curOffset = elemManager->GetRegion(0)->GetSubRegion(0)->numNodesPerElement();
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          localIndex offSetForOneCell = elemSubRegion->numNodesPerElement();
          for( localIndex i =  0; i < elemSubRegion->size(); i++ )
          {
            offsetFragment[countOffsetFragmentIndex++] = curOffset;   
            curOffset += offSetForOneCell;
            if( countOffsetFragmentIndex == multiplier )
            {
              stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( offsetFragment.data() ), outputString, sizeof( localIndex ) * multiplier );
              countOffsetFragmentIndex = 0;
            }
          }
        });
      });
      outputString.resize( FindBase64StringLength( sizeof( localIndex ) * ( countOffsetFragmentIndex) ) );
      stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( offsetFragment.data() ), outputString, sizeof( localIndex ) * ( countOffsetFragmentIndex) );
      DumpBuffer( stream );
    }

    void WriteAsciiTypes( ElementRegionManager const * const elemManager )
    {
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          integer type = geosxToVTKCellTypeMap.at( elemSubRegion->GetElementTypeString() );
          for( localIndex i =  0; i < elemSubRegion->size(); i++ )
          {
            m_outFile <<  type << "\n";
          }
        });
      });
    }

    void WriteBinaryTypes( ElementRegionManager const * const elemManager )
    {
      std::stringstream stream;
      integer multiplier = FindMultiplier( sizeof( integer ) ); // We do not write all the data at once to avoid creating a big table each time.
      integer_array typeFragment( multiplier );
      string outputString;
      outputString.resize( FindBase64StringLength( sizeof( integer ) * multiplier ) );
      integer countTypeFragmentIndex = 0;
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
        {
          integer type = geosxToVTKCellTypeMap.at( elemSubRegion->GetElementTypeString() );
          for( localIndex i =  0; i < elemSubRegion->size(); i++ )
          {
            typeFragment[countTypeFragmentIndex++] = type;   
            if( countTypeFragmentIndex == multiplier )
            {
              stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( typeFragment.data() ), outputString, sizeof( integer ) * multiplier );
              countTypeFragmentIndex = 0;
            }
          }
        });
      });
      outputString.resize( FindBase64StringLength( sizeof( integer ) * ( countTypeFragmentIndex) ) );
      stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( typeFragment.data() ), outputString, sizeof( localIndex ) * ( countTypeFragmentIndex) );
      DumpBuffer( stream );
    }

    template< typename T >
    void WriteCellAsciiData( ElementRegionManager::ElementViewAccessor< T > const & dataView,
                             ElementRegionManager const * const elemManager )
    {
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const er,
                                                                    auto const * const elemRegion )
      {
        elemRegion->template forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                                     auto const * const elemSubRegion )
        {
          for ( localIndex ei = 0; ei < elemSubRegion->size(); ++ei )
          {
            LvArray::forValuesInSlice( dataView[er][esr][ei], [this]( auto const & value ) { m_outFile << value << " "; } );
            m_outFile << "\n";
          }
        });
      });
    }

    template< typename ARRAY_TYPE >
    void WriteCellBinaryData( ElementRegionManager::ElementViewAccessor< ARRAY_TYPE > const & dataView, ElementRegionManager const * const elemManager )
    {
      using VALUE_TYPE = typename ARRAY_TYPE::value_type;

      std::stringstream stream;
      WriteSize( elemManager->getNumberOfElements< CellElementSubRegion >(), sizeof( VALUE_TYPE ) );
      integer multiplier = FindMultiplier( sizeof( VALUE_TYPE ) );// We do not write all the data at once to avoid creating a big table each time.
      string outputString;
      outputString.resize( FindBase64StringLength( sizeof( VALUE_TYPE ) * multiplier ) );
      std::vector< VALUE_TYPE > dataFragment( multiplier );
      integer countDataFragment = 0;
      elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const er,
                                                                        CellElementRegion const * const elemRegion )
      {
        elemRegion->template forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                                     CellElementSubRegion const * const elemSubRegion )
        {
          if ( elemSubRegion->size() > 0 )
          {
            for( localIndex ei = 0; ei < elemSubRegion->size(); ++ei )
            {
              LvArray::forValuesInSlice( dataView[er][esr][ei],
                [&]( VALUE_TYPE const & value )
                {
                  dataFragment[ countDataFragment++ ] = value;
                  if( countDataFragment == multiplier )
                  {
                    stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( dataFragment.data() ), outputString, sizeof( VALUE_TYPE ) * countDataFragment );
                    countDataFragment = 0;
                  }
                }
              );
            }
          }
          else
          {
            real64 nanArray[3] = { std::nan("0"), std::nan("0"), std::nan("0") };
            for( localIndex ei = 0; ei  < elemSubRegion->size(); ei++)
            {
              stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( nanArray ), outputString, sizeof( real64 ) * 3 );
            }
          }
        });
      });

      outputString.resize( FindBase64StringLength( sizeof( VALUE_TYPE ) * ( countDataFragment) ) );
      stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( dataFragment.data() ), outputString, sizeof( VALUE_TYPE ) * ( countDataFragment ) );
      DumpBuffer( stream );
    }

    template< typename ARRAY_TYPE >
    void WriteNodeAsciiData( Wrapper< ARRAY_TYPE > const & dataView )
    {
      ARRAY_TYPE const & array = dataView.reference();
      for( localIndex i = 0; i < array.size( 0 ); i++ )
      {
        LvArray::forValuesInSlice( array[ i ], [this]( auto const & value ) { m_outFile << value << " "; } );
        m_outFile << "\n";
      }
    }

    template< typename ARRAY_TYPE >
    void WriteNodeBinaryData( Wrapper< ARRAY_TYPE > const & dataView )
    {

      std::stringstream stream;
      ARRAY_TYPE const & array = dataView.reference();

      using VALUE_TYPE = typename ARRAY_TYPE::value_type;
      
      integer multiplier = FindMultiplier( sizeof( VALUE_TYPE ) );// We do not write all the data at once to avoid creating a big table each time.
      
      string outputString;
      outputString.resize( FindBase64StringLength( sizeof( VALUE_TYPE ) * multiplier ) );
      
      std::vector< VALUE_TYPE > dataFragment( multiplier );
      
      integer countDataFragment = 0;
      WriteSize( array.size(), sizeof( VALUE_TYPE ) );
      for( localIndex i = 0; i < array.size( 0 ); i++ )
      {
        LvArray::forValuesInSlice( array[ i ],
          [&]( VALUE_TYPE const & value )
          {
            dataFragment[countDataFragment++] = value;
            if( countDataFragment == multiplier )
            {
              stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( dataFragment.data() ), outputString, sizeof( VALUE_TYPE ) * countDataFragment );
              countDataFragment = 0;
            }
          }
        );
      }

      outputString.resize( FindBase64StringLength( sizeof( VALUE_TYPE ) * countDataFragment ) );
      stream << stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( dataFragment.data() ), outputString, sizeof( VALUE_TYPE ) * countDataFragment );
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

    integer FindBase64StringLength( integer dataSize )
    {
      integer base64StringLength = (dataSize * 8) / 6;
      while( base64StringLength % 4 )
      {
        base64StringLength++;
      }
      return base64StringLength;
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
inline void CustomVTUXMLWriter::WriteCellBinaryData( ElementRegionManager::ElementViewAccessor< r1_array > const & dataView, ElementRegionManager const * const elemManager )
{
  std::stringstream stream;
  string outputString;
  outputString.resize(  FindBase64StringLength( sizeof( real64 ) * 3) );
  WriteSize( elemManager->getNumberOfElements< CellElementSubRegion >() * 3, sizeof( real64 ) );
  elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                                 auto const * const elemSubRegion )
    {
      if ( dataView[er][esr].size() > 0 )
        for( localIndex ei = 0; ei  < elemSubRegion->size(); ei++)
        {
          stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( dataView[er][esr][ei].Data() ), outputString, sizeof( real64 ) * 3 );
        }
      else
      {
        real64_array nanArray(3);
        nanArray[0] = nanArray[1] = nanArray[2] = std::nan("0");
        for( localIndex ei = 0; ei  < elemSubRegion->size(); ei++)
        {
          stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( nanArray.data() ), outputString, sizeof( real64 ) * 3 );
        }
      }
    });
  });
  DumpBuffer( stream );
}

template<>
inline void CustomVTUXMLWriter::WriteNodeBinaryData( Wrapper< r1_array > const & dataView )
{
  std::stringstream stream;
  auto & viewRef = dataView.reference();
  string outputString;
  outputString.resize(  FindBase64StringLength( sizeof( real64 ) * 3 ) );
  WriteSize( viewRef.size() * 3, sizeof( real64 ) );
  for( localIndex i = 0; i < viewRef.size(); i++)
  {
    stream <<stringutilities::EncodeBase64( reinterpret_cast< const unsigned char * >( viewRef[i].Data() ), outputString, sizeof( real64 ) * 3 );
  }
  DumpBuffer( stream );
}

VTKFile::VTKFile( string const & name ):
  m_baseName( name ),
  m_binary( false )
{
  int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
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
  int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
  int const mpiSize = MpiWrapper::Comm_size(MPI_COMM_GEOSX);
  ElementRegionManager const * elemManager = domain.getMeshBody(0)->getMeshLevel(0)->getElemManager();
  NodeManager const * nodeManager = domain.getMeshBody(0)->getMeshLevel(0)->getNodeManager();
  string timeStepFolderName = m_baseName + "/" + std::to_string( timeStep );
  if( mpiRank == 0 )    
  {                     
    // Create a directory for this time step
    mode_t mode = 0733;
    mkdir( timeStepFolderName.c_str(), mode );
  }
  MpiWrapper::Barrier();
  string format;
  if( m_binary )
  {
    format = "binary";
  }  
  else
  {
    format = "ascii";
  }

  std::set< std::tuple< string, string, integer, rtTypes::TypeIDs > > cellFields; // First : field name, Second : type, Third : field dimension;
  // Find all cell fields to export
  elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
  {
    elemRegion->forElementSubRegions([&]( auto const * const subRegion )
    {
      for( auto const & wrapperIter : subRegion->wrappers() )
      {
        WrapperBase const * const wrapper = wrapperIter.second;

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
          cellFields.insert(std::make_tuple(fieldName, geosxToVTKTypeMap.at( typeID ), dimension, fieldType) );
        }
      }
    });
  });

  std::set< std::tuple< string, string, integer, rtTypes::TypeIDs > > nodeFields; // First : field name, Second : type, Third : field dimension;
  // Find all node fields to export
  for( auto const & wrapperIter : nodeManager->wrappers() )
  {
    WrapperBase const * const wrapper = wrapperIter.second;
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
       nodeFields.insert( std::make_tuple( fieldName, geosxToVTKTypeMap.at( typeID ), dimension, fieldType ) );
    }
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
    
    // Declare all the point fields
    auto pPointDataNode = pUnstructureGridNode.append_child("PPointData");
    for( auto & nodeField : nodeFields )
    {
      CreatePDataArray(pPointDataNode, std::get<1>(nodeField), std::get<0>(nodeField), std::get<2>(nodeField), format);
    }

    // Declaration of the node PCells
    auto pCellsNode = pUnstructureGridNode.append_child("PCells");
    // .... and its data array defining the connectivities, types, and offsets
    CreatePDataArray( pCellsNode, geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ), "connectivity", 1, format ); //TODO harcoded for the moment
    CreatePDataArray( pCellsNode, geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ), "offsets", 1, format );
    CreatePDataArray( pCellsNode, geosxToVTKTypeMap.at( std::type_index( typeid( integer ) ) ), "types", 1, format );

    // Find all the cell fields to output
    auto pCellDataNode = pUnstructureGridNode.append_child("PCellData");
    for( auto & cellField : cellFields )
    {
      CreatePDataArray(pCellDataNode, std::get<1>(cellField), std::get<0>(cellField), std::get<2>(cellField), format);
    }

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
  localIndex totalNumberOfCells = elemManager->getNumberOfElements< CellElementSubRegion >();
  localIndex totalNumberOfSubRegion = 0;
  elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
  {
    totalNumberOfSubRegion += elemRegion->numSubRegions();
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
  for( auto & nodeField : nodeFields )
  {
    WrapperBase const * const wrapper = nodeManager->getWrapperBase( std::get<0>( nodeField ) );
    vtuWriter.OpenXMLNode( "DataArray", { { "type", std::get<1>( nodeField ) },
                                          { "Name", std::get<0>( nodeField ) },
                                          { "NumberOfComponents", std::to_string( std::get<2>( nodeField ) ) },
                                          { "format", format } } );
    rtTypes::ApplyArrayTypeLambda1( std::get<3>( nodeField ),
                                    [&]( auto type ) -> void
    {
      using cType = decltype(type);
      const Wrapper< cType > & view = Wrapper<cType>::cast( *wrapper );
      vtuWriter.WriteNodeData( view, m_binary );
    });
    vtuWriter.CloseXMLNode( "DataArray" );
  }
  vtuWriter.CloseXMLNode( "PointData" );

  // Definition of the node Cells
  vtuWriter.OpenXMLNode( "Cells", {} );

  // Definition of the node DataArray that will contain the connectivities
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ) },
                                        { "Name", "connectivity" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  
  vtuWriter.WriteCellConnectivities( elemManager, m_binary );

  /*
  elemManager->forElementRegionsComplete< ElementRegion >( [&]( localIndex const er,
                                                                auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellConnectivities( elemSubRegion->GetElementTypeString(), elemSubRegion->nodeList(), m_binary );
    });
  });
  */
  vtuWriter.CloseXMLNode( "DataArray" );


  array1d< std::tuple< integer, localIndex, string > > subRegionsInfo; // First value : cell size, Second value : number of cells, Third value : cell Types
  elemManager->forElementRegionsComplete< CellElementRegion >( [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                                    auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      subRegionsInfo.push_back( std::make_tuple( elemSubRegion->numNodesPerElement(), elemSubRegion->size(), elemSubRegion->GetElementTypeString() ) );
    });
  });

  // Definition of the node DataArray that will contain the offsets
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( localIndex ) ) ) },
                                        { "Name", "offsets" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  vtuWriter.WriteCellOffsets( elemManager, m_binary );
  vtuWriter.CloseXMLNode( "DataArray" );

  // Definition of the node DataArray that will contain the cell types
  vtuWriter.OpenXMLNode( "DataArray", { { "type", geosxToVTKTypeMap.at( std::type_index( typeid( integer ) ) ) },
                                        { "Name", "types" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  vtuWriter.WriteCellTypes( elemManager, m_binary );
  vtuWriter.CloseXMLNode( "DataArray" );

  vtuWriter.CloseXMLNode( "Cells" );
  // Definition of the CellDataArray node that will contains all the data held by the elements
  vtuWriter.OpenXMLNode( "CellData", {} );
  for( auto & cellField : cellFields )
  {
    vtuWriter.OpenXMLNode( "DataArray", { { "type", std::get<1>( cellField )  },
                                           { "Name", std::get<0>( cellField )  },
                                           { "NumberOfComponents", std::to_string( std::get<2>( cellField )  ) },
                                           { "format", format } } );
    rtTypes::ApplyArrayTypeLambda1( std::get<3>( cellField ),
                                    [&]( auto type ) -> void
    {
      using cType = decltype(type);
      auto dataView = elemManager->ConstructViewAccessor< cType >(std::get<0>( cellField ));
      vtuWriter.WriteCellData( dataView, elemManager, m_binary );
    });
    vtuWriter.CloseXMLNode( "DataArray" );
  }
  vtuWriter.CloseXMLNode( "CellData" );
  vtuWriter.CloseXMLNode( "Piece" );
  vtuWriter.CloseXMLNode( "UnstructuredGrid" );
  vtuWriter.CloseXMLNode( "VTKFile" );
}

}
