/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file PvtuFile.cpp
 */

#include <iostream>
#include <string.h>

#include "PvtuFile.hpp"

#include "common/Logger.hpp"
#include "mesh/MeshBody.hpp"

#if USE_MPI
#include <mpi.h>
#endif

namespace geosx{

// PUBLIC METHODS
void PvtuFile::Load( string const &filename) {
    int mpiSize = 0;
    int mpiRank = 0;
#if USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
    int numFilesPerRank = 0;
    //int nb_total_partitions = 0;
    int remainderFiles = 0;
    std::vector < string > children_files;



    // XML file parsing using pugixml tool
    pugi::xml_document pvtu_doc;
    pvtu_doc.load_file(filename.c_str());

    if( mpiRank == 0) {
        CheckXmlParentFileConsistency(pvtu_doc,filename);
    }
    VtuFilesList(pvtu_doc,children_files);
    // Retrieve the number of partitions
    int const numFiles = children_files.size();

    // Next part of this method is dedicated to the optimization of file loading
    //
    // IF numPartitions > nb_mpi_process : Each process will load 
    // numPartitions / nb_mpi_process. Processes with a rank < numPartitions % nb_mpi_process
    // will load an additional partition
    //
    // IF numPartitions < nb_mpi_process : the first numPartitions process will
    // load ONE partition.
    if(numFiles > mpiSize) {
        if ( mpiRank == 0 ) {
            std::cout << "WARNING : the number of partitions ("
                << numFiles <<") which will be loaded " 
                << "is greater that the number of processes on which GEOSX is running ("
                << mpiSize
                << "). Some processes will hold more than one partition !" << std::endl;
        }
        numFilesPerRank = numFiles / mpiSize;
        remainderFiles = numFiles % mpiSize;
    } else {
        numFilesPerRank = 1;
        remainderFiles = 0;
    }
    if( mpiRank < remainderFiles ) {
        m_vtuFileNames.resize(numFilesPerRank+1);
        m_vtuFiles.resize(numFilesPerRank+1);
        for(int p_index = 0; p_index < numFilesPerRank +1; ++p_index) {
            m_vtuFileNames[p_index] =
                children_files[mpiRank*(numFilesPerRank+1) + p_index];
        }
    } else if( mpiRank < numFiles) {
        m_vtuFileNames.resize(numFilesPerRank);
        m_vtuFiles.resize(numFilesPerRank);
        for(int p_index = 0; p_index < numFilesPerRank; ++p_index) {
            m_vtuFileNames[p_index] =
                children_files[mpiRank*(numFilesPerRank) + p_index+remainderFiles];
        }
    }
    if( mpiRank < numFiles ) {
        for(int p_index = 0; p_index < 
                static_cast< int >(m_vtuFileNames.size()); ++p_index) {
            m_vtuFiles[p_index].Load(m_vtuFileNames[p_index]);
        }
    }
}

void PvtuFile::Save( string const &filename) {
    GEOS_ERROR("pvtu file save is not implemented yet");
}

void VtuFile::Load( string const &filename) {
    pugi::xml_document pvtu_doc;
    pvtu_doc.load_file(filename.c_str());
    CheckXmlChildFileConsistency(pvtu_doc,filename);
    LoadMeshPart(pvtu_doc);
}

void VtuFile::Save( string const &filename) {
    GEOS_ERROR("vtu file save is not implemented yet");
}


//PRIVATE METHODS
void PvtuFile::CheckXmlFileConsistency(pugi::xml_document const & pvtu_doc,
        string const & filename,
        string const & prefix) const {

    // VTKFile is the main node of the pvtufile
    auto const vtk_file =pvtu_doc.child("VTKFile");
    if( vtk_file.empty() ) {
        GEOS_ERROR("Main node VTKFile not found in " + filename);
    }

    string ugrid_name = prefix +"UnstructuredGrid";
    auto const ugrid = vtk_file.child(ugrid_name.c_str());
    if( ugrid.empty() ) {
        GEOS_ERROR("Node " + prefix + "UnstructuredGrid not found or empty in " +
                filename);
    }
    pugi::xml_node main_node;
    if( prefix == "P" ) {
        main_node = ugrid;
    } else {
        main_node = ugrid.child("Piece");
        if (main_node.empty()) {
            GEOS_ERROR("Piece node is missing in " + filename);
        }
    }

    string const point_data_name = prefix + "PointData";
    string const cell_data_name = prefix + "CellData";
    string const points_name = prefix + "Points";
    string const main_node_child_names[3] = {point_data_name,
        cell_data_name,points_name};
    for( auto const & main_node_child_name : main_node_child_names ) {
        auto const main_node_child = main_node.child(main_node_child_name.c_str());
        if(main_node_child.empty()) {
            GEOS_ERROR("Node " +main_node_child_name +
                " not found or empty in " + filename);
        };
    }

    bool points_have_original_index_property = false;
    string const mandatory_attributes[3] =
    {"Name","type","format"};
    for( auto const & pdata_property : main_node.child(point_data_name.c_str()).children() ) {
        if( pdata_property.name() == static_cast< string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    pdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    GEOS_ERROR("Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of " + prefix + "PointData");
                }
            }
            if( pdata_property.attribute("Name").value() == m_strOriginalIndex) {
                points_have_original_index_property = true;
            }
        }
    }
    if (!points_have_original_index_property) {
        GEOS_ERROR("Can't find any DataArray which contains the property " +
            m_strOriginalIndex + " in " + prefix + "PointData in " + filename);
    }

    bool cells_have_original_index_property = false;
    bool cells_have_partition_property = false;
    bool cells_have_region_property = false;
    for( auto const & pdata_property : main_node.child(cell_data_name.c_str()).children() ) {
        if( pdata_property.name() == static_cast< string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    pdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    GEOS_ERROR("Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of " + prefix + "CellData in "
                        + filename);
                }
            }
            if( pdata_property.attribute("Name").value() == m_strOriginalIndex) {
                cells_have_original_index_property = true;
            }
            if( pdata_property.attribute("Name").value() == m_strPartition) {
                cells_have_partition_property = true;
            }
            if( pdata_property.attribute("Name").value() == m_strRegion) {
                cells_have_region_property = true;
            }
        }
    }
    if (!cells_have_original_index_property) {
        GEOS_ERROR("Can't find any DataArray which contains the property "+
                m_strOriginalIndex + " in " + prefix + "CellData in "+filename);
    }
    if (!cells_have_region_property) {
        GEOS_ERROR("Can't find any DataArray which contains the property "+
                m_strRegion + " in "+ prefix + "CellData in " + filename);
    }
    if (!cells_have_partition_property) {
        GEOS_ERROR("Can't find any DataArray which contains the property "+
                m_strPartition + " in " + prefix + "CellData in " + filename);
    }

    bool point_has_a_pdata_array = false;
    bool point_has_a_pdata_array_with_points = false;
    string const data_array_name = prefix + "DataArray";
    for(auto const & point_child : main_node.child(points_name.c_str()).children() ) {
        if( point_child.name() == data_array_name ) {
            point_has_a_pdata_array = true;
            if( point_child.attribute("Name").as_string() ==
                    static_cast< string >("Points")) {
                point_has_a_pdata_array_with_points = true;
                if (point_child.attribute("NumberOfComponents").as_uint() != 3 ) {
                    GEOS_ERROR("GEOSX supports only 3D meshes");
                }
                break;
            }
        }
    }
    if( !point_has_a_pdata_array ) {
        GEOS_ERROR("No " + prefix +"DataArray found in " + prefix + "Points in " + filename);
    }
    if( !point_has_a_pdata_array_with_points ) {
        GEOS_ERROR("No " + prefix +"DataArray named \"Points\" found in " + filename);
    }

}
void PvtuFile::CheckXmlParentFileConsistency(pugi::xml_document const & pvtu_doc,
        string const & filename) const {

        CheckXmlFileConsistency( pvtu_doc, filename, "P");

        bool has_a_piece = false;
        for( auto const & ugrid_child :
                pvtu_doc.child("VTKFile").child("PUnstructuredGrid").children() ) {
            if(ugrid_child.name() == static_cast< string >("Piece")) {
                has_a_piece = true;
                if( ugrid_child.attribute("Source").empty() ) {
                    GEOS_ERROR("Piece nodes has to have an attribute Source not empty.");
                }
            }
        }
        if( !has_a_piece ) {
            GEOS_ERROR("No Piece not found in " + filename );
        }
}
void PvtuFile::VtuFilesList(
        pugi::xml_document const & pvtu_doc,
        std::vector < string > & vtu_files ) const{
    int numPartitions = 0;
    for(auto const & child : pvtu_doc.child("VTKFile").child("PUnstructuredGrid").children()) {
        if( child.name() == static_cast< string > ("Piece") )
        {
            vtu_files.emplace_back(child.attribute("Source").as_string());
        }
    }
}

void VtuFile::CheckXmlChildFileConsistency(pugi::xml_document const & pvtu_doc,
        string const & filename) const {
    CheckXmlFileConsistency( pvtu_doc, filename);

    pugi::xml_node piece_node =
        pvtu_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");
    if( piece_node.attribute("NumberOfPoints").empty() ) {
        GEOS_ERROR("Attribute \"NumberOfPoints\" of Node \"Piece\" is missing or empty in "
                + filename);
    }

    if( piece_node.attribute("NumberOfCells").empty() ) {
        GEOS_ERROR("Attribute \"NumberOfCells\" of Node \"Piece\" is missing or empty in "
                + filename);
    }

    pugi::xml_node cell_node = piece_node.child("Cells");
    string const mandatory_attributes[3] =
    {"Name","type","format"};
    bool cells_have_connectivity = false;
    bool cells_have_type = false;
    bool cells_have_offset = false;
    for( auto cell_att : cell_node.children() ) {
        if( cell_att.name() == static_cast< string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    cell_att.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    GEOS_ERROR("Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of in "
                        + filename);
                }
            }
            if( cell_att.attribute("Name").value() ==
                    static_cast<string>("connectivity")) {
                cells_have_connectivity = true;
            }
            if( cell_att.attribute("Name").value() == static_cast< string >("offsets")) {
                cells_have_offset = true;
            }
            if( cell_att.attribute("Name").value() == static_cast< string > ("types")) {
                cells_have_type = true;
            }
        }
    }
    if (!cells_have_connectivity) {
        GEOS_ERROR("No property \"connectivity\" found in \"Cells\" node of " + filename);
    }
    if (!cells_have_type) {
        GEOS_ERROR("No property \"types\" found in \"Cells\" node of " + filename);
    }
    if (!cells_have_offset) {
        GEOS_ERROR("No property \"offsets\" found in \"Cells\" node of " + filename);
    }
}   

template<typename T, typename Lambda>
void VtuFile::SplitNodeTextString( string const & in,
        std::vector< T >& out,
        Lambda && string_convertor) const {
    std::stringstream stringStream(in);
    string line;
    while(std::getline(stringStream, line))
    {
        std::size_t prev = 0, pos;
        while ((pos = line.find_first_of(" ", prev)) != string::npos)
        {
            if (pos > prev)
                out.push_back(string_convertor(line.substr(prev, pos-prev)));
            prev = pos+1;
        }
        if (prev < line.length())
            out.push_back(string_convertor(line.substr(prev, string::npos)));
    }
}

void VtuFile::LoadMeshPart(pugi::xml_document const & pvtu_doc){
    pugi::xml_node piece_node =
        pvtu_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

    globalIndex numVertices = piece_node.attribute("NumberOfPoints").as_llong();
    m_meshPart.SetNumVertices( numVertices );

    globalIndex numElements = piece_node.attribute("NumberOfCells").as_llong();
    m_meshPart.ReserveNumCellAndPolygons( numElements );

    /// Parse vertices
    pugi::xml_node m_verticesarray =
        piece_node.child("Points").find_child_by_attribute("DataArray","Name","Points");
    assert(!m_verticesarray.empty()) ;
    std::vector< real64 > all_vertices;
    all_vertices.reserve( numVertices*3 );
    SplitNodeTextString(m_verticesarray.text().as_string(), all_vertices, 
            [](string str)-> double {return std::stod(str);});
    assert(static_cast< globalIndex> (all_vertices.size()) / 3 == numVertices);

    /// Parse vertices original index
    pugi::xml_node m_verticesoriginal_indexes_array = 
        piece_node.child("PointData").find_child_by_attribute(
                "DataArray","Name","original_index");
    assert(!m_verticesoriginal_indexes_array.empty());
    std::vector< globalIndex > all_m_verticesoriginal_indexes;
    all_m_verticesoriginal_indexes.reserve( numVertices );
    SplitNodeTextString(m_verticesoriginal_indexes_array.text().as_string(),
            all_m_verticesoriginal_indexes,
            [](string str)-> globalIndex {return std::stoll(str);});
    assert( numVertices == static_cast< globalIndex> (all_m_verticesoriginal_indexes.size() ));

    /// Fill the vertices in the mesh
    for(globalIndex vertexIndex  = 0 ; vertexIndex < numVertices; ++vertexIndex) {
        std::vector< real64 > vertex(3);
        for(localIndex coor = 0 ; coor < 3 ; ++coor) {
            vertex[coor] = all_vertices[3*vertexIndex+coor];
        }
        m_meshPart.SetVertexOriginalIndex(vertexIndex,
                all_m_verticesoriginal_indexes[vertexIndex]);
        m_meshPart.SetVertex(vertexIndex, vertex);
    }

    /// Parse elements types
    pugi::xml_node elements_types_array =
        piece_node.child("Cells").find_child_by_attribute("DataArray","Name","types");
    assert(!elements_types_array.empty());
    std::vector< localIndex > all_elements_types;
    all_elements_types.reserve(numElements);
    SplitNodeTextString( elements_types_array.text().as_string(), all_elements_types,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (all_elements_types.size()) == numElements);

    /// Parse elements regions
    pugi::xml_node elements_regions_array =
        piece_node.child("CellData").find_child_by_attribute("DataArray","Name","region");
    assert(!elements_regions_array.empty());
    std::vector< globalIndex > all_elements_regions;
    all_elements_regions.reserve(numElements);
    SplitNodeTextString( elements_regions_array.text().as_string(), all_elements_regions,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (all_elements_regions.size()) == numElements);
    
    /// Parse elements original_index
    pugi::xml_node elements_original_indexes_array =
        piece_node.child("CellData").find_child_by_attribute("DataArray","Name",
                "original_index");
    assert(!elements_original_indexes_array.empty());
    std::vector< globalIndex > all_elements_original_indexes;
    all_elements_original_indexes.reserve(numElements);
    SplitNodeTextString( elements_original_indexes_array.text().as_string(),
            all_elements_original_indexes,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (all_elements_original_indexes.size()) == numElements);

    /// Parse elements offsets
    std::vector< globalIndex> elements_offsets(numElements);
    pugi::xml_node elements_offsets_array =
        piece_node.child("Cells").find_child_by_attribute("DataArray","Name","offsets");
    assert(!elements_offsets_array.empty());
    std::vector< globalIndex > all_elements_offsets(1,0); // convenient to have 0 as first offset
    all_elements_offsets.reserve(numElements);
    SplitNodeTextString( elements_offsets_array.text().as_string(), all_elements_offsets,
            [](string str)-> globalIndex {return std::stoll(str);});
    assert(static_cast<globalIndex >(all_elements_offsets.size()) == numElements+1);

    /// Parce cells connectivities
    pugi::xml_node elements_connectivity_array =
        piece_node.child("Cells").find_child_by_attribute("DataArray","Name","connectivity");
    assert(!elements_connectivity_array.empty());
    std::vector< globalIndex > all_elements_connectivities;
    all_elements_connectivities.reserve(8 * numElements);
    SplitNodeTextString(elements_connectivity_array.text().as_string(),
            all_elements_connectivities,
            [](string str)-> localIndex {return std::stoi(str);});
    for( globalIndex elementIndex = 0 ; elementIndex < numElements; ++elementIndex) {
        localIndex numVerticesInCell = all_elements_offsets[elementIndex+1] -
            all_elements_offsets[elementIndex];
        std::vector< globalIndex > connectivity(numVerticesInCell);
        for(localIndex co = 0 ; co < numVerticesInCell ; co++) {
            connectivity[co] =
                all_elements_connectivities[all_elements_offsets[elementIndex]+co];
        }
        if( all_elements_types[elementIndex] == 10 ||
                all_elements_types[elementIndex] == 14 ||
                all_elements_types[elementIndex] ==12 ||
                all_elements_types[elementIndex] == 13 ) {
            m_meshPart.AddCell(connectivity);
            m_meshPart.SetCellRegion(elementIndex, all_elements_regions[elementIndex]);
            m_meshPart.SetCellOriginalIndex(elementIndex,
                    all_elements_original_indexes[elementIndex]);
        }
        else if( all_elements_types[elementIndex] == 5 ||
                all_elements_types[elementIndex] == 9) {
            m_meshPart.AddPolygon(connectivity);
            m_meshPart.SetPolygonSurface(elementIndex,
                    all_elements_regions[elementIndex]);
            m_meshPart.SetPolygonOriginalIndex(elementIndex,
                    all_elements_original_indexes[elementIndex]);
        }
    }

    m_meshPart.Finish();
}

void VtuFile::TransferMeshPartToGEOSMesh( MeshLevel * const meshLevel )
{
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  ElementRegionManager * const elemRegMananger = meshLevel->getElemManager();

  arrayView1d<R1Tensor> X = nodeManager->referencePosition();

  nodeManager->resize(m_meshPart.NumVertices());

  real64 const * const vertexData = m_meshPart.m_vertices.data();

  for( localIndex a=0 ; a<m_meshPart.NumVertices() ; ++a )
  {
    real64 * const tensorData = X[a].Data();
    tensorData[0] = vertexData[3*a];
    tensorData[1] = vertexData[3*a+1];
    tensorData[2] = vertexData[3*a+2];
  }

  string regionName; // = mesh_part.regionName();
  CellBlockSubRegion * const cellBlock = elemRegMananger->GetRegion( regionName )->GetSubRegion(0);
  lArray2d & cellToVertex = cellBlock->nodeList();
  cellToVertex.resize( 0, m_meshPart.NumVerticesInCell(0) );
  cellBlock->resize( m_meshPart.NumCells() );

  for( localIndex k=0 ; k<m_meshPart.NumCells() ; ++k )
  {
    for( localIndex a=0 ; a<m_meshPart.NumVerticesInCell(0) ; ++a )
    {
      cellToVertex[k][a] = m_meshPart.CellVertexIndex(k,a);
    }
  }



}


    ////////////////
    /// MESH PART //
    ////////////////

    globalIndex MeshPart::NumVertices() const {
        return m_numVertices;
    }

    globalIndex MeshPart::NumCells() const {
        return m_numCells;
    }

    globalIndex MeshPart::NumPolygons() const {
        return m_numPolygons;
    }

    localIndex MeshPart::NumVerticesInCell( const globalIndex cellIndex ) const {
        assert(cellIndex < m_numCells);
        return m_cellsPtr[cellIndex+1] - m_cellsPtr[cellIndex-1];
    }

    localIndex MeshPart::NumVerticesInPolygon( const globalIndex polygonIndex ) const {
        assert(polygonIndex < m_numPolygons);
        return m_polygonsPtr[polygonIndex+1] - m_polygonsPtr[polygonIndex-1];
    }

    globalIndex MeshPart::CellVertexIndex(globalIndex const cellIndex,
            localIndex const local_corner_index) const {
        assert(local_corner_index < NumVerticesInCell(cellIndex));
        return m_cellsConnectivity[m_cellsPtr[cellIndex] + local_corner_index];
    }

    globalIndex MeshPart::PolygonVertexIndex(globalIndex const polygonIndex,
            localIndex const local_corner_index) const {
        assert(local_corner_index < NumVerticesInPolygon(polygonIndex));
        return m_polygonsConnectivity[m_polygonsPtr[polygonIndex] + local_corner_index];
    }

    std::vector<real64> MeshPart::Vertex(globalIndex const vertexIndex) const {
        assert(vertexIndex < m_numVertices);
        std::vector<real64> vertex(m_vertices.begin()+3*vertexIndex,
                m_vertices.begin() + 3*vertexIndex+3);
        return vertex;
    }

    globalIndex MeshPart::Surface(globalIndex const polygonIndex) const {
        assert(polygonIndex < m_numPolygons);
        return m_surfaceIndexes[polygonIndex];
    }

    globalIndex MeshPart::Region(globalIndex const cellIndex) const {
        assert(cellIndex < m_numCells);
        return m_regionIndexes[cellIndex];
    }

    globalIndex MeshPart::GlobalVertexIndex(globalIndex const vertexIndex) const {
        assert(vertexIndex < m_numVertices);
        return m_originalVertexIndexes[vertexIndex];
    }

    globalIndex MeshPart::GlobalPolygonIndex(globalIndex const polygonIndex) const {
        assert(polygonIndex<m_numPolygons);
        return m_originalPolygonIndexes[polygonIndex];
    }

    globalIndex MeshPart::GlobalCellIndex(globalIndex const cellIndex) const {
        assert(cellIndex < m_numCells);
        return m_originalCellIndexes[cellIndex];
    }

    void MeshPart::SetVertex(globalIndex const vertexIndex,
            std::vector< real64 > const& vertex) {
        assert(vertex.size() == 3); 
        assert(vertexIndex < m_numVertices);
        for( localIndex coor = 0 ; coor < 3 ; coor ++ ) {
            m_vertices[3*vertexIndex+coor] = vertex[coor];
        }
    }

    void MeshPart::SetVertexOriginalIndex( globalIndex const vertexIndex,
            globalIndex const original_index) {
        assert(vertexIndex < m_numVertices);
        m_originalVertexIndexes[vertexIndex] = original_index;
    }

    void MeshPart::SetNumVertices(globalIndex const numVertices) {
        m_numVertices = numVertices;
        m_vertices.resize( numVertices*3);
        m_originalVertexIndexes.resize( numVertices );
    }

    void MeshPart::ReserveNumCellAndPolygons(globalIndex const numElements) {
        m_cellsPtr.reserve( numElements +1);
        m_polygonsPtr.reserve( numElements +1);
        m_cellsConnectivity.reserve( 8 * numElements ); // maximum 8 corners (for an hex)
        m_polygonsConnectivity.reserve( 4 * numElements ); // maximum 4 corners (for a quad)
        m_surfaceIndexes.reserve( numElements);
        m_regionIndexes.reserve( numElements);
        m_originalPolygonIndexes.reserve( numElements);
        m_originalCellIndexes.reserve( numElements);
    }
    
    globalIndex MeshPart::AddCell( std::vector<globalIndex> connectivity ) {
        m_cellsPtr.push_back( connectivity.size() + m_cellsPtr[m_cellsPtr.size()-1]);
        m_numCells++;
        m_regionIndexes.resize(m_numCells);
        m_originalCellIndexes.resize(m_numCells);
        for( localIndex co = 0 ; co < static_cast<localIndex>(connectivity.size()); ++co ) {
            m_cellsConnectivity.push_back(connectivity[co]);
        }
        return m_numCells-1;
    }

    globalIndex MeshPart::AddPolygon( std::vector<globalIndex> connectivity ) {
        m_polygonsPtr.push_back( connectivity.size() + m_polygonsPtr[m_polygonsPtr.size()-1]);
        m_numPolygons++;
        m_surfaceIndexes.resize(m_numPolygons);
        m_originalPolygonIndexes.resize( m_numPolygons);
        for( localIndex co = 0 ; co < static_cast<localIndex>(connectivity.size()); ++co ) {
            m_polygonsConnectivity.push_back(connectivity[co]);
        }
        return m_numPolygons-1;
    }

    void MeshPart::SetCellRegion( globalIndex const cellIndex,
            globalIndex const region_index) {
        assert( cellIndex < static_cast< globalIndex> (m_regionIndexes.size()) );
        m_regionIndexes[cellIndex]  = region_index;
    }

    void MeshPart::SetPolygonSurface( globalIndex const polygonIndex,
            globalIndex const surface_index) {
        assert( polygonIndex < static_cast< globalIndex> (m_surfaceIndexes.size() ) );
        m_surfaceIndexes[polygonIndex]  = surface_index;
    }

    void MeshPart::SetCellOriginalIndex(globalIndex const cellIndexInPartMesh,
                globalIndex const cellIndexInFullMesh) {
        assert( cellIndexInPartMesh < static_cast< globalIndex>
                (m_originalCellIndexes.size() ) );
        m_originalCellIndexes[cellIndexInPartMesh] = cellIndexInFullMesh;
    }

    void MeshPart::SetPolygonOriginalIndex(globalIndex const polygonIndexInPartMesh,
                globalIndex const polygonIndexInFullMesh) {
        assert( polygonIndexInPartMesh <
                static_cast< globalIndex> (m_originalPolygonIndexes.size() ) );
        m_originalPolygonIndexes[polygonIndexInPartMesh] = polygonIndexInFullMesh;
    }   

    void MeshPart::Finish() {
        m_regionIndexes.shrink_to_fit();
        m_surfaceIndexes.shrink_to_fit();
        m_cellsPtr.shrink_to_fit();
        m_polygonsPtr.shrink_to_fit();
        m_cellsConnectivity.shrink_to_fit();
        m_polygonsConnectivity.shrink_to_fit();
        m_originalCellIndexes.shrink_to_fit();
        m_originalPolygonIndexes.shrink_to_fit();
    }
}
