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
  VTKFile::VTKFile( string const & name ):
    m_baseName( name ),
    m_binary( false ),
    m_compress( false )
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
      CreatePDataArray( pPointsNode, m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::real64_id), "Position", 3 );

      // Declaration of the node PCells
      auto pCellsNode = pUnstructureGridNode.append_child("PCells");
      // .... and its data array defining the connectivities, types, and offsets
      CreatePDataArray( pCellsNode, m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::localIndex_id), "connectivity", 1 );
      CreatePDataArray( pCellsNode, m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::localIndex_id), "offsets", 1 );
      CreatePDataArray( pCellsNode, m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::integer_id), "types", 1 );

      // Find all the fields to output
      auto pCellDataNode = pUnstructureGridNode.append_child("PCellData");
      elemManager->forElementRegionsComplete< ElementRegion, FaceElementRegion >( [&]( localIndex const er,
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
              if( !m_geosxToVTKTypeMap.count(fieldType) )
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
              CreatePDataArray(pCellDataNode, m_geosxToVTKTypeMap.at(fieldType), fieldName, dimension);
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
  elemManager->forElementRegionsComplete< ElementRegion, FaceElementRegion >( [&]( localIndex const er,
                                                                              auto const * const elemRegion )
  {
    totalNumberOfCells += elemRegion->GetTotalSize();
  });
  vtuWriter.OpenXMLNode( "Piece", { { "NumberOfPoints", std::to_string(nodeManager->size() ) },
                                    { "NumberOfCells", std::to_string( totalNumberOfCells ) } } );

  // Definition of node Points
  vtuWriter.OpenXMLNode( "Points",{} );

  // Definition of the node DataArray that will contain all the node coordinates
  vtuWriter.OpenXMLNode( "DataArray", { { "type", m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::real64_id) },
                                        { "Name", "Position" },
                                        { "NumberOfComponents", "3" },
                                        { "format", format } } );
  vtuWriter.WriteVertices( nodeManager->referencePosition(), m_binary );
  vtuWriter.CloseXMLNode( "DataArray" );
  vtuWriter.CloseXMLNode( "Points" );

  // Definition of the node Cells
  vtuWriter.OpenXMLNode( "Cells", {} );

  // Definition of the node DataArray that will contain the connectivities
  vtuWriter.OpenXMLNode( "DataArray", { { "type", m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::localIndex_id ) },
                                        { "Name", "connectivity" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  elemManager->forElementRegionsComplete< ElementRegion, FaceElementRegion >( [&]( localIndex const er,
                                                                              auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellConnectivities( elemSubRegion->nodeList(), m_binary );
    });
  });
  vtuWriter.CloseXMLNode( "DataArray" );

  // Definition of the node DataArray that will contain the offsets
  vtuWriter.OpenXMLNode( "DataArray", { { "type", m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::localIndex_id ) },
                                        { "Name", "offsets" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  elemManager->forElementRegionsComplete< ElementRegion, FaceElementRegion >( [&]( localIndex const er,
                                                                              auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellOffsets( elemSubRegion->numNodesPerElement(), elemSubRegion->size(), m_binary );
    });
  });
  vtuWriter.CloseXMLNode( "DataArray" );

  // Definition of the node DataArray that will contain the cell types
  vtuWriter.OpenXMLNode( "DataArray", { { "type", m_geosxToVTKTypeMap.at( rtTypes::TypeIDs::localIndex_id ) },
                                        { "Name", "types" },
                                        { "NumberOfComponents", "1" },
                                        { "format", format } } );
  elemManager->forElementRegionsComplete< ElementRegion, FaceElementRegion >( [&]( localIndex const er,
                                                                              auto const * const elemRegion )
  {
    elemRegion->template forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&]( auto const * const elemSubRegion )
    {
      vtuWriter.WriteCellTypes( elemSubRegion->GetElementTypeString(), elemSubRegion->size(), m_binary );
    });
  });

  vtuWriter.CloseXMLNode( "DataArray" );
  vtuWriter.CloseXMLNode( "Cells" );

  // Definition of the CellDataArray node that will contains all the data held by the elements
  vtuWriter.OpenXMLNode( "CellData", {} );
  elemManager->forElementRegionsComplete< ElementRegion, FaceElementRegion >( [&]( localIndex const er,
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
          if( !m_geosxToVTKTypeMap.count(fieldType) )
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
          vtuWriter.OpenXMLNode( "DataArray", { { "type", m_geosxToVTKTypeMap.at( fieldType ) },
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
