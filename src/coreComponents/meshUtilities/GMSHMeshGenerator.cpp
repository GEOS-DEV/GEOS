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
 * @file GMSHMeshGenerator.cpp
 */

#include "GMSHMeshGenerator.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "managers/DomainPartition.hpp"

#include <fstream>
#include <math.h>
#include <string>
#include <unordered_map>

#include "mpiCommunications/PartitionBase.hpp"
#include "mpiCommunications/SpatialPartition.hpp"

#include "mesh/MeshBody.hpp"

#include "parmetis.h"

namespace geosx
{
using namespace dataRepository;

GMSHMeshGenerator::GMSHMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_lineNumber(0)
{

  registerWrapper( viewKeyStruct::filePathString, &m_filePath )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "path to the mesh file" );
  registerWrapper( viewKeyStruct::fieldsToImportString, &m_fieldsToImport )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Fields to be imported from the external mesh file" );
  registerWrapper( viewKeyStruct::fieldNamesInGEOSXString, &m_fieldNamesInGEOSX )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the fields within GEOSX" );
  registerWrapper( viewKeyStruct::scaleString, &m_scale )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDefaultValue( 1. )->setDescription( "Scale the coordinates of the vertices" );
  registerWrapper( viewKeyStruct::reverseZString, &m_isZReverse )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDefaultValue( 0 )->setDescription( "0 : Z coordinate is upward, 1 : Z coordinate is downward" );
  registerWrapper( viewKeyStruct::initNbOfProcString, &m_initNbOfProc )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDefaultValue( 1 )->setDescription( "Set the number of processors on which the mesh will be loaded before being balancer on every processors." );
}

GMSHMeshGenerator::~GMSHMeshGenerator()
{}

void GMSHMeshGenerator::GenerateElementRegions( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void GMSHMeshGenerator::PostProcessInput()
{
}

void GMSHMeshGenerator::RemapMesh( dataRepository::Group * const GEOSX_UNUSED_PARAM( domain ) )
{
  return;
}

Group * GMSHMeshGenerator::CreateChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void GMSHMeshGenerator::GetLine(std::ifstream & fileStream, std::string& line) 
{
  do
  {
    if( !std::getline( fileStream, line ) )
    {
      return;
    }
    m_lineNumber++;
  } while (line.empty() );
}

void GMSHMeshGenerator::GenerateMesh( DomainPartition * const GEOSX_UNUSED_PARAM(domain) )
{ 
  int const mpiRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

  std::unordered_map< localIndex, std::string > cellElementRegionsIndexToName;

  // Two following arrays contains temporary informations which 
  // will be stored on the first m_initNbProcs ranks
  array1d< globalIndex > nodeGlobalIndexes;
  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodeCoordinates;

  // Three following array contains temporary information
  // on the mesh whic will be stored by the first m_initNbProcs
  array1d< globalIndex > elementGlobalIndex;  // size : number of element in  the first m_initNbProcs ranks
  array1d< localIndex >  elementPhysicalIds;  // size : number of element in  the first m_initNbProcs ranks
  array1d< idx_t > elementPtr;          // size : number of element + 1 in the first m_initNbProcs ranks 
  array1d< idx_t > elementConnectivity; // size : not known a priori.


  if( mpiRank < m_initNbOfProc )
  {
    // Opening file
    //GEOSX_LOG_RANK_0("Begin import of " << m_filePath );
    std::ifstream inputStream( m_filePath);
    GEOSX_ERROR_IF( !inputStream, "Could not read input file: " << m_filePath );

    std::string lineString;
    while( !inputStream.eof()  )
    { 
     // if( lineString.substr(0,4) == "$End" ) continue;
      GetLine( inputStream, lineString );
      if( mpiRank == 0 ) //Basic information read by rank 0
      {
        if(lineString == "$MeshFormat" ) // Reading mesh format
        {
          GetLine( inputStream, lineString );
          double version;
          int isAscii;
          int dataSize;
          std::istringstream iss(lineString);
          GEOSX_ERROR_IF( !( iss >> version >> isAscii >> dataSize ), " [line " << m_lineNumber << "] 1 float and 2 int are expected (version, ascii info and dataSize)");
        }

        if( lineString == "$PhysicalNames") // Readin physical entities (surface and volume of interests for the simulation)
        {
          GetLine( inputStream, lineString );
          std::istringstream iss(lineString);
          localIndex nbPhysicalEntities;
          int dim;
          int physicalEntityIndex;
          string physicalEntityName;
          GEOSX_ERROR_IF( !( iss >> nbPhysicalEntities ), " [line " << m_lineNumber << "] 1 integer is expected corresponding to the total number of physical entities");
          for( localIndex pE = 0; pE < nbPhysicalEntities; pE++)
          {
            GetLine( inputStream, lineString );
            std::istringstream iss2(lineString);
            GEOSX_ERROR_IF( !( iss2 >> dim >> physicalEntityIndex >> physicalEntityName ), " [line " << m_lineNumber << "] two integers (dimension and index) and one string (name) are expected");
            if( dim == 3 )
            {
              cellElementRegionsIndexToName[physicalEntityIndex] = physicalEntityName;
            }
          }
        }
      }

      // Next part is execute by every m_initNbOfProc
      if( lineString == "$Nodes" )
      {
        GetLine( inputStream, lineString );
        std::istringstream iss(lineString);
        globalIndex nbTotalNodes;
        globalIndex gIndex;
        real64 x;
        real64 y;
        real64 z;
        GEOSX_ERROR_IF( !( iss >> nbTotalNodes ), " [line " << m_lineNumber << "] 1 integer is expected corresponding to the total number of nodes");
        globalIndex nbNodesPerRank = nbTotalNodes / m_initNbOfProc;
        globalIndex remainder = nbTotalNodes % m_initNbOfProc;
        if( mpiRank == m_initNbOfProc  -1)
        {
          nodeCoordinates.resize( nbNodesPerRank + remainder, 3 ) ; // TODO : maybe we can do a better balance?
          nodeGlobalIndexes.resize( nbNodesPerRank + remainder );
        }
        else
        {
          nodeCoordinates.resize( nbNodesPerRank, 3 );
          nodeGlobalIndexes.resize( nbNodesPerRank );
        }
        array1d< globalIndex > offsets( m_initNbOfProc + 1);
        for( int i = 0; i < m_initNbOfProc; i++)
        {
          offsets[i] = i * nbNodesPerRank;
        }
        offsets[m_initNbOfProc] = nbTotalNodes;
        globalIndex tempLocalVertexIndex = 0; //Local vertex index for the m_initNbRanks
        for( globalIndex v = 0; v < nbTotalNodes; v++ )
        {
          GetLine( inputStream, lineString );
          if( v >= offsets[mpiRank] && v < offsets[mpiRank +1] )
          {
            std::istringstream iss2(lineString);
            GEOSX_ERROR_IF( !( iss2 >> gIndex >> x >> y >> z ), "[line " << m_lineNumber << "] expecting node definition with a global index and three float coordinates ");
            nodeCoordinates(tempLocalVertexIndex, 0) = x;
            nodeCoordinates(tempLocalVertexIndex, 1) = y;
            nodeCoordinates(tempLocalVertexIndex, 2) = z;
            nodeGlobalIndexes[tempLocalVertexIndex] = gIndex;
            tempLocalVertexIndex++;
          }
        }
      }

      if( lineString == "$Elements" )
      {
        GetLine( inputStream, lineString );
        std::istringstream iss(lineString);
        globalIndex nbTotalElements;
        globalIndex gIndex;
        globalIndex elementType;
        localIndex nbOfTags;
        localIndex physicalIndex;
        localIndex trash;
        GEOSX_ERROR_IF( !( iss >> nbTotalElements ), " [line " << m_lineNumber << "] 1 integer is expected corresponding to the total number of elements");
        //GEOSX_LOG_RANK( "nb total of elements " << nbTotalElements);

        /// First we need to count the number of volume elements
        localIndex nbTotalVolumeElements = 0;
        std::streampos savePos = inputStream.tellg();
        for( localIndex e = 0; e < nbTotalElements; e++ )
        {
          GetLine( inputStream, lineString );
          std::istringstream iss2(lineString);
          GEOSX_ERROR_IF( !( iss2 >> gIndex >> elementType >> nbOfTags), "[line " << m_lineNumber << "] expecting element index, type and number of tags");
            //GEOSX_LOG_RANK( "element type lol " << elementType );
          if( m_gmshCellTypeToNbNodes.count( elementType ) )
          {
            //GEOSX_LOG_RANK( "ca passe");
            nbTotalVolumeElements++;
          }
        }
        //GEOSX_LOG_RANK( "nb total volume elements " << nbTotalVolumeElements );
        inputStream.seekg( savePos ); // rewind
        globalIndex nbElementsPerRankAPriori = nbTotalElements / m_initNbOfProc;
        globalIndex remainder = nbTotalVolumeElements % m_initNbOfProc;
        if( mpiRank == m_initNbOfProc  -1)
        {
          elementGlobalIndex.reserve( nbElementsPerRankAPriori + remainder ) ; // TODO : maybe we can do a better balance?
          elementPhysicalIds.reserve( nbElementsPerRankAPriori + remainder ) ;
          elementPtr.reserve( nbElementsPerRankAPriori + remainder + 1);
          elementPtr.emplace_back( 0 );
          elementConnectivity.reserve( (nbElementsPerRankAPriori + remainder) * 8 );

        }
        else
        {
          elementGlobalIndex.reserve( nbElementsPerRankAPriori ) ; // TODO : maybe we can do a better balance?
          elementPhysicalIds.reserve( nbElementsPerRankAPriori ) ;
          elementPtr.reserve( nbElementsPerRankAPriori  + 1);
          elementPtr.emplace_back( 0 );
          elementConnectivity.reserve( nbElementsPerRankAPriori );
        }
        array1d< globalIndex > offsets( m_initNbOfProc + 1);
        for( int i = 0; i < m_initNbOfProc; i++)
        {
          offsets[i] = i * nbElementsPerRankAPriori;
        }
        offsets[m_initNbOfProc] = nbTotalElements;
        globalIndex tempLocalElementIndex = 0; //Local vertex index for the m_initNbRanks
        //GEOSX_LOG_RANK( "offsets " << offsets );
        for( globalIndex e = 0; e < nbTotalElements; e++ )
        {
          GetLine( inputStream, lineString );
          if( e >= offsets[mpiRank] && e < offsets[mpiRank +1] )
          {
            std::istringstream iss2(lineString);
            GEOSX_ERROR_IF( !( iss2 >> gIndex >> elementType >> nbOfTags), "[line " << m_lineNumber << "] expecting element index, type and number of tags");
            //GEOSX_LOG_RANK( "element type " << elementType );
            if( m_gmshCellTypeToNbNodes.count( elementType ) )
            {
              if( nbOfTags == 0 )
              {
                physicalIndex = 0;
              }
              else
              {
                GEOSX_ERROR_IF( !( iss2 >> physicalIndex ), "[line " << m_lineNumber << "] expecting element index, type and number of tags" );
                for( localIndex i = 1; i < nbOfTags; i++)
                {
                  GEOSX_ERROR_IF( !( iss2 >> trash ), "[line " << m_lineNumber << "] can't read other element tags" );
                }
              }
              localIndex nbVerticesInElement = m_gmshCellTypeToNbNodes.at(elementType);
              for( localIndex v = 0; v < nbVerticesInElement; v++ )
              {
                globalIndex vGId;
                GEOSX_ERROR_IF( !( iss2 >> vGId ), "[line " << m_lineNumber << "] wrong element connectivity" );
                
                elementConnectivity.emplace_back( vGId );
              }
              //GEOSX_LOG_RANK(1);
              elementGlobalIndex.emplace_back(gIndex);
              //GEOSX_LOG_RANK(2);
              elementPhysicalIds.emplace_back(physicalIndex);
              //GEOSX_LOG_RANK(3);
              elementPtr.emplace_back( elementPtr[elementPtr.size() - 1] + nbVerticesInElement );
              //GEOSX_LOG_RANK(4);
              tempLocalElementIndex++;
            }
          }
        }
        //GEOSX_LOG_RANK( "global indexes " << elementGlobalIndex );
        //GEOSX_LOG_RANK( "physical ids " << elementPhysicalIds );
        //GEOSX_LOG_RANK( "elementPtr " << elementPtr );
        //GEOSX_LOG_RANK( "connectivity " <<elementConnectivity );
      }
    }
  }

  // Call to Parmetis
  
  //Build elmdist
  array1d< idx_t > elmdist( mpiSize + 1);
  idx_t nbElementsPerRank = elementGlobalIndex.size();
  MpiWrapper::Allgather( &nbElementsPerRank, 1, elmdist.data() + 1, 1, MPI_COMM_GEOSX );

  for( int r = 0; r < mpiSize; r ++)
  {
    elmdist[r+1] = elmdist[r+1] + elmdist[r];
  }
  //GEOSX_LOG_RANK( "elmdist " <<elmdist );

  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t ncommonnodes = 4;
  idx_t nparts = mpiSize;
  array1d< real_t > tpwgts( nparts * ncon );
  for(int i = 0; i < tpwgts.size(); i++ )
  {
    tpwgts[i] = 1./nparts;
  }
  real_t ubvec = 1.05; // from doc
  idx_t options = 0;
  idx_t edgecuts;
  array1d< idx_t > part(nbElementsPerRank);
  if( mpiRank >= m_initNbOfProc )
  {
    elementPtr.emplace_back(0);
    elementPtr.emplace_back(0);
    elementConnectivity.emplace_back(-1);// will never be read
    part.emplace_back( -1 ); // will never be read or overwritten
  }
  //GEOSX_LOG_RANK( "part " << part );
  //GEOSX_LOG_RANK( "elmdist " << elmdist );
  //GEOSX_LOG_RANK( "elementPtr " << elementPtr );
  //GEOSX_LOG_RANK( "elementConnectivity " << elementConnectivity );
  MPI_Comm parmetisComm;
  MPI_Group geosxGroup;
  MPI_Group parmetisGroup;
  array1d< int > listOfRankWithElements( m_initNbOfProc);
  for(int r = 0; r < m_initNbOfProc; r++)
  {
    listOfRankWithElements[r] = r;
  }
  MPI_Comm_group( MPI_COMM_GEOSX, &geosxGroup);
  MPI_Group_incl( geosxGroup, m_initNbOfProc, listOfRankWithElements.data(), &parmetisGroup );

  MPI_Comm_create_group(MPI_COMM_GEOSX, parmetisGroup, 0, &parmetisComm );
  if( MPI_COMM_NULL != parmetisComm )
  {
    ParMETIS_V3_PartMeshKway( elmdist.data(), elementPtr.data(), elementConnectivity.data(),
        NULL, &wgtflag, &numflag, &ncon, &ncommonnodes,&nparts,
        tpwgts.data(), &ubvec, &options, &edgecuts, part.data(), &parmetisComm );
    MPI_Comm_free(&parmetisComm);
    MPI_Group_free(&parmetisGroup);
  }
  //GEOSX_LOG_RANK( "part are : " << part );

  // Ordering of elements to be sent to other ranks
  //TODO : it's a little bit overkill. Maybe it's possible to do something better
  array1d< localIndex > nbElementOwnByRank( mpiSize );
  array1d< localIndex > connectivitySizeOwnByRank( mpiSize );
  std::vector< array1d< localIndex > > orderedElementPhysicalIds( mpiSize );
  std::vector< array1d< globalIndex > > orderedElementGlobalIndex( mpiSize );
  std::vector< array1d< globalIndex > > orderedElementConnectivity( mpiSize );
  std::vector< array1d< globalIndex > > orderedElementPtr( mpiSize );
  if(mpiRank < m_initNbOfProc)
  {
    for(int r = 0; r < mpiSize; r++)
    {
      orderedElementPhysicalIds[r].reserve( elementPhysicalIds.size() / mpiSize);
      orderedElementGlobalIndex[r].reserve( elementGlobalIndex.size() / mpiSize);
      orderedElementConnectivity[r].reserve( elementConnectivity.size() / mpiSize );
      orderedElementPtr[r].reserve( elementPtr.size() / mpiSize );
      orderedElementPtr[r].emplace_back( 0 );
    }
    for(int e = 0; e < part.size(); e++)
    {
      orderedElementGlobalIndex[part[e]].emplace_back( elementGlobalIndex[e] );
      orderedElementPhysicalIds[part[e]].emplace_back( elementPhysicalIds[e] );
      for( int v = elementPtr[e]; v < elementPtr[e+1]; v++)
      {
        orderedElementConnectivity[part[e]].emplace_back( elementConnectivity[v] );
      }
      orderedElementPtr[part[e]].emplace_back( orderedElementPtr[part[e]][ orderedElementPtr[part[e]].size()  -1] + elementPtr[e+1] - elementPtr[e] );
    }
    for(int r = 0; r < mpiSize; r++)
    {
      nbElementOwnByRank[r] = orderedElementGlobalIndex[r].size();
      connectivitySizeOwnByRank[r] = orderedElementConnectivity[r].size();
      //GEOSX_LOG_RANK( "to " << r << " phys id " << orderedElementPhysicalIds[r] );
      //GEOSX_LOG_RANK( "to " << r << " g id " << orderedElementGlobalIndex[r] );
      //GEOSX_LOG_RANK( "to " << r << " conn  " << orderedElementConnectivity[r] );
      //GEOSX_LOG_RANK( "to " << r << "  ptr " << orderedElementPtr[r] );
      //localIndex nbElements = orderedElementGlobalIndex[r].size();
      //localIndex connectivitySize = orderedElementConnectivity[r].size();
      // Sending everything !!
      //MpiWrapper::Send( orderedElementPtr[r].data(), orderedElementPtr.size(), r, 0, MPI_COMM_GEOSX );
    }
  }
  //GEOSX_LOG_RANK( "nb elements own by ranl " << nbElementOwnByRank);
  // Receive number of element
  array1d< globalIndex > gatheredNbOfElements( m_initNbOfProc );
  array1d< globalIndex > gatheredConnectivitySize( m_initNbOfProc );
  //GEOSX_LOG_RANK_0("begin gather");
  for(int r = 0; r < mpiSize; r ++)
  {
    MpiWrapper::gather( &nbElementOwnByRank[r], 1, gatheredNbOfElements.data(), 1, r, MPI_COMM_GEOSX );
    MpiWrapper::gather( &connectivitySizeOwnByRank[r], 1, gatheredConnectivitySize.data(), 1, r, MPI_COMM_GEOSX );
  }
  //GEOSX_LOG_RANK("end gather " << gatheredNbOfElements);
  // Compute final number of element for EVERY rank
  localIndex finalNumberOfElements = 0;
  localIndex finalConnectivitySize = 0;
  for( int r = 0; r < m_initNbOfProc; r++ )
  {
    finalNumberOfElements += gatheredNbOfElements[r];
    finalConnectivitySize += gatheredConnectivitySize[r];
  }
  //GEOSX_LOG_RANK( "final number of elements " << finalNumberOfElements);
  //GEOSX_LOG_RANK( "final connectivty size " << finalConnectivitySize);

  // Receiving everything
  array1d< globalIndex > finalElementGlobalIndex( finalNumberOfElements );
  array1d< localIndex >  finalElementPhysicalIds( finalNumberOfElements );
  array1d< globalIndex > finalElementPtr( finalNumberOfElements + 1);         
  array1d< globalIndex > finalElementConnectivity( finalConnectivitySize ); 

  // Send global index
  if ( mpiRank < m_initNbOfProc )
  {
    for( int r = 0; r < mpiSize; r++ )
    {
      MpiWrapper::Send( orderedElementGlobalIndex[r].data(),orderedElementGlobalIndex[r].size(),r,0,MPI_COMM_GEOSX);
    }
  }
  globalIndex offset = 0;
  for( int r = 0; r < m_initNbOfProc; r++ )
  {
    MpiWrapper::Recv(finalElementGlobalIndex.data() + offset,gatheredNbOfElements[r],r,0,MPI_COMM_GEOSX,MPI_STATUS_IGNORE );
    offset += gatheredNbOfElements[r];
  }

  // Send physicalId
  if ( mpiRank < m_initNbOfProc )
  {
    for( int r = 0; r < mpiSize; r++ )
    {
      MpiWrapper::Send( orderedElementPhysicalIds[r].data(),orderedElementPhysicalIds[r].size(),r,0,MPI_COMM_GEOSX);
    }
  }
  offset = 0;
  for( int r = 0; r < m_initNbOfProc; r++ )
  {
    MpiWrapper::Recv(finalElementPhysicalIds.data() + offset,gatheredNbOfElements[r],r,0,MPI_COMM_GEOSX,MPI_STATUS_IGNORE );
    offset += gatheredNbOfElements[r];
  }
  GEOSX_LOG_RANK( "final global Index " << finalElementGlobalIndex );
  GEOSX_LOG_RANK( "final physical id " << finalElementPhysicalIds );
  GEOSX_LOG_RANK( "FINIIIIIIISH");
  MpiWrapper::Barrier();

                    
}

void GMSHMeshGenerator::GetElemToNodesRelationInBox( const std::string & GEOSX_UNUSED_PARAM( elementType ),
                                                       const int GEOSX_UNUSED_PARAM( index )[],
                                                       const int & GEOSX_UNUSED_PARAM( iEle ),
                                                       int GEOSX_UNUSED_PARAM( nodeIDInBox )[],
                                                       const int GEOSX_UNUSED_PARAM( node_size ) )
{}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, GMSHMeshGenerator, std::string const &, Group * const )
}
