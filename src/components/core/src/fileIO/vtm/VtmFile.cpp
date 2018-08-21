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

#include "VtmFile.hpp"

#include "common/Logger.hpp"
#include "mesh/MeshBody.hpp"

#if USE_MPI
#include <mpi.h>
#endif

namespace geosx{

MeshBlock::MeshBlock( string fileName,
        string blockName ) :
    m_vtuFileName(fileName),
    m_blockName(blockName)
{
}

bool MeshBlock::IsARegionBlock() const  {
    return m_mesh.NumCells() > 0;
}

void MeshBlock::Load() {
    m_vtuFile.Load(m_vtuFileName,m_mesh);    
    m_mesh.SetName(m_blockName);
}

DumbMesh const & MeshBlock::mesh() const {
    return m_mesh;
}

void RankBlock::AddMeshBlock( const MeshBlock& block) {
    m_block.emplace_back(block);
}

void RankBlock::Load() {
    for( auto& block : m_block) {
        block.Load();
    }
}
localIndex RankBlock::NumMeshBlocks() const {
    return static_cast< localIndex > (m_block.size());
}
MeshBlock const & RankBlock::GetMeshBlock(localIndex const meshBlockIndex) const {
    assert( meshBlockIndex < static_cast< localIndex > (m_block.size()));
    return m_block[ meshBlockIndex ];
}
// PUBLIC METHODS

RankBlock const & VtmFile::GetRankBlock(localIndex const rankBlockIndex) const {
    assert( rankBlockIndex < static_cast< localIndex > (m_rankBlocks.size()));
    return m_rankBlocks[rankBlockIndex];
}

localIndex VtmFile::NumRankBlocks() const{
    return static_cast< localIndex >(m_rankBlocks.size());
}
void VtmFile::Load( string const &filename) {
    int mpiSize = 0;
    int mpiRank = 0;
#if USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
    int numFilesPerRank = 0;
    //int nb_total_partitions = 0;
    int remainderFiles = 0;
    // XML file parsing using pugixml tool
    pugi::xml_document vtmDoc;
    vtmDoc.load_file(filename.c_str());

    if( mpiRank == 0) {
        CheckXmlFileConsistency(vtmDoc,filename);
    }
    std::vector< RankBlock > rankBlocks;
    SetRanksAndBlocks(vtmDoc,rankBlocks);
    // Retrieve the number of partitions
    globalIndex const numFiles = rankBlocks.size();

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
        for(int p_index = 0; p_index < numFilesPerRank +1; ++p_index) {
            m_rankBlocks.emplace_back(rankBlocks[mpiRank*(numFilesPerRank+1) + p_index]);
        }
    } else if( mpiRank < numFiles) {
        for(int p_index = 0; p_index < numFilesPerRank; ++p_index) {
            m_rankBlocks.emplace_back(rankBlocks[remainderFiles+mpiRank*(numFilesPerRank) + p_index]);
        }
    }
    if( mpiRank < numFiles ) {
        for(int p_index = 0; p_index < 
                static_cast< int >(m_rankBlocks.size()); ++p_index) {
            m_rankBlocks[p_index].Load();
        }
    }
}

void VtmFile::FromVtmToGEOS(MeshLevel * const meshLevel) {
    for( const auto & rankBlock : m_rankBlocks ) {
        rankBlock.TransferRankBlockToGEOSMesh(meshLevel);
    }
}

void VtuFile::Load( string const &filename, DumbMesh &mesh) {
    pugi::xml_document vtmDoc;
    vtmDoc.load_file(filename.c_str());
    CheckXmlChildFileConsistency(vtmDoc,filename);
    LoadMesh(vtmDoc, mesh);
}

void VtuFile::Save( string const &filename) {
    GEOS_ERROR("vtu file save is not implemented yet");
}


//PRIVATE METHODS
void VtmFile::CheckXmlFileConsistency(pugi::xml_document const & vtmDoc,
        string const & filename) const {

    // VTKFile is the main node of the pvtufile
    auto const & vtkFileNode =vtmDoc.child("VTKFile");
    if( vtkFileNode.empty() ) {
        GEOS_ERROR("Main node VTKFile not found in " + filename);
    }

    if( vtkFileNode.attribute("type").as_string() !=
            static_cast< string >("vtkMultiBlockDataSet")) {
        GEOS_ERROR("VTKFile is not a vtkMultiBlockDataSet" + filename);
    }

    auto const & vtkMultiBlockDataSetNode =  vtkFileNode.child("vtkMultiBlockDataSet");
    if( vtkMultiBlockDataSetNode.empty() ) {
        GEOS_ERROR("vtkMultiBlockDataSet node is not present in " + filename);
    }

    globalIndex countNbRank=0;
    for( auto const & rank : vtkMultiBlockDataSetNode.children()) {
        if( rank.name() == static_cast< string >("Block")) {
            countNbRank++;
            localIndex countNbBlock =0;
            for( auto const & block : rank.children()) {
                if( block.name() == static_cast< string >("DataSet")) {
                    countNbBlock++;
                }
                if( block.attribute("file").empty() ) {
                    GEOS_ERROR(
                            "DataSet "
                            + static_cast< string >(rank.attribute("name").as_string())
                            + " of block "
                            +static_cast< string >(rank.attribute("name").as_string()) +
                            " does not contain a \"file\" attribute");
                }
                if( block.attribute("name").empty() ) {
                    GEOS_ERROR(
                            "DataSet "
                            + static_cast< string >(rank.attribute("name").as_string())
                            + " of block "
                            +static_cast< string >(rank.attribute("name").as_string()) +
                            " does not contain a \"name\" attribute");
                }
            }
            if( countNbBlock == 0 ) {
                GEOS_ERROR(static_cast< string >(rank.attribute("name").as_string()) +
                        " does not contain any DataSet");
            }
        }
    }
    if( countNbRank == 0 ) {
        GEOS_ERROR("There is no block defined in " + filename);
    }

}

void VtmFile::SetRanksAndBlocks(
        pugi::xml_document const & vtmDoc,
                std::vector< RankBlock >& rankBlocks) {
    auto const & vtkFileNode =vtmDoc.child("VTKFile");
    auto const & vtkMultiBlockDataSetNode =  vtkFileNode.child("vtkMultiBlockDataSet");
    
    for(auto const & rank : vtkMultiBlockDataSetNode.children()) {
        if( rank.name() == static_cast< string > ("Block") )
        {
            RankBlock curRankBlock;
            for( auto const & block : rank.children() ) {
                if( block.name() == static_cast< string > ("DataSet")) {
                    MeshBlock curMeshBlock(
                            block.attribute("file").as_string(),
                            block.attribute("name").as_string());
                    curRankBlock.AddMeshBlock(curMeshBlock);
                }
            }
            rankBlocks.emplace_back(curRankBlock);
        }
    }
}

void VtuFile::CheckXmlChildFileConsistency(pugi::xml_document const & vtmDoc,
        string const & filename) const {

    auto const & vtkFileNode =vtmDoc.child("VTKFile");
    auto const & uGridNode = vtkFileNode.child("UnstructuredGrid");
    if( uGridNode.empty() ) {
        GEOS_ERROR("Node UnstructuredGrid not found or empty in " +
                filename);
    }
    pugi::xml_node pieceNode = uGridNode.child("Piece");
    if (pieceNode.empty()) {
        GEOS_ERROR("Piece node is missing in " + filename);
    }

    string const pieceNodeChildNames[4] = {"PointData",
        "CellData","Points","Cells"};
    for( auto const & pieceNodeChildName : pieceNodeChildNames ) {
        auto const & pieceNodeChild = pieceNode.child(pieceNodeChildName.c_str());
        if(pieceNodeChild.empty()) {
            GEOS_ERROR("Node " + pieceNodeChildName +
                " not found or empty in " + filename);
        };
    }

    string const mandatoryDataAttributes[3] =
    {"Name","type","format"};
    for( auto const & dataProperty : pieceNode.child("PointData").children() ) {
        if( dataProperty.name() == static_cast< string >("DataArray")) {
            for( auto & mandatoryAttribute : mandatoryDataAttributes ){
                auto const & attribute =
                    dataProperty.attribute( mandatoryAttribute.c_str());
                if( attribute.empty() ) {
                    GEOS_ERROR("Mandatory attribute " + mandatoryAttribute +
                        " does not exist in a DataArray of PointData");
                }
            }
        }
    }

    bool cellsHasGlobalIndexProperty = false;
    for( auto const & dataProperty : pieceNode.child("CellData").children() ) {
        if( dataProperty.name() == static_cast< string >("DataArray")) {
            for( auto  & mandatoryAttribute : mandatoryDataAttributes ){
                auto const & attribute =
                    dataProperty.attribute( mandatoryAttribute.c_str());
                if( attribute.empty() ) {
                    GEOS_ERROR("Mandatory attribute " + mandatoryAttribute +
                        " does not exist in a DataArray of CellData in "
                        + filename);
                }
            }
            if( dataProperty.attribute("Name").value() ==
                    static_cast<string>("globalIndex")) {
                cellsHasGlobalIndexProperty = true;
            }
        }
    }
    if (!cellsHasGlobalIndexProperty) {
        GEOS_ERROR("Can't find any DataArray which contains the property originalIndex in  CellData in "+filename);
    }

    auto const & pointDataArray = pieceNode.child("Points").find_child_by_attribute("DataArray","Name","Points");
    if(pointDataArray.empty() ) {
        GEOS_ERROR("Can't find DataArray names \"Points\" in node \"Points\" in " + filename);
    }

    if( pieceNode.attribute("NumberOfCells").empty() ) {
        GEOS_ERROR("Attribute \"NumberOfCells\" of Node \"Piece\" is missing or empty in "
                + filename);
    }

    if( pieceNode.attribute("NumberOfPoints").empty() ) {
        GEOS_ERROR("Attribute \"NumberOfPoints\" of Node \"Piece\" is missing or empty in "
                + filename);
    }

    if (pieceNode.child("Cells").find_child_by_attribute(
                "DataArray","Name","connectivity").empty()) {
        GEOS_ERROR("No property \"connectivity\" found in \"Cells\" node of " + filename);
    }
    if (pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","types").empty()) {
        GEOS_ERROR("No property \"types\" found in \"Cells\" node of " + filename);
    }
    if (pieceNode.child("Cells").find_child_by_attribute(
                "DataArray","Name","offsets").empty()) {
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

template < typename Enum >
auto to_underlying_type( Enum e ) ->
typename std::underlying_type< Enum >::type
{
    return static_cast< typename std::underlying_type< Enum >::type >( e );
}

void VtuFile::LoadMesh(pugi::xml_document const & vtmDoc, DumbMesh& mesh){
    pugi::xml_node pieceNode =
        vtmDoc.child("VTKFile").child("UnstructuredGrid").child("Piece");

    globalIndex numVertices = pieceNode.attribute("NumberOfPoints").as_llong();
    mesh.SetNumVertices( numVertices );

    globalIndex numElements = pieceNode.attribute("NumberOfCells").as_llong();

    /// Parse vertices
    pugi::xml_node m_verticesarray =
        pieceNode.child("Points").find_child_by_attribute("DataArray","Name","Points");
    assert(!m_verticesarray.empty()) ;
    std::vector< real64 > allVertices;
    allVertices.reserve( numVertices*3 );
    SplitNodeTextString(m_verticesarray.text().as_string(), allVertices, 
            [](string str)-> double {return std::stod(str);});
    assert(static_cast< globalIndex> (allVertices.size()) / 3 == numVertices);

    /// Fill the vertices in the mesh
    for(globalIndex vertexIndex  = 0 ; vertexIndex < numVertices; ++vertexIndex) {
        std::vector< real64 > vertex(3);
        for(localIndex coor = 0 ; coor < 3 ; ++coor) {
            vertex[coor] = allVertices[3*vertexIndex+coor];
        }
        mesh.SetVertex(vertexIndex, vertex);
    }

    /// Parse elements types
    pugi::xml_node elementsTypesArray =
        pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","types");
    assert(!elementsTypesArray.empty());
    std::vector< localIndex > allElementsTypes;
    allElementsTypes.reserve(numElements);
    SplitNodeTextString( elementsTypesArray.text().as_string(), allElementsTypes,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (allElementsTypes.size()) == numElements);
    globalIndex numTetra = 0;
    globalIndex numHex = 0;
    globalIndex numPrism = 0;
    globalIndex numPyr = 0;
    globalIndex numTri = 0;
    globalIndex numQuad = 0;
    for( auto elementType : allElementsTypes ) {
        if( elementType == 10 ) {
            numTetra++;
        } else if(elementType == 12) {
            numHex++;
        } else if(elementType == 13) {
            numPrism++;
        } else if(elementType == 14) {
            numPyr++;
        } else if(elementType == 5) {
            numTri++;
        } else if(elementType == 9) {
            numQuad++;
        } else {
            GEOS_ERROR("Element type " + std::to_string(elementType) + " not supported");
        }
    }

    mesh.SetNumCellAndPolygons(numTetra,numHex,numPrism,numPyr,numTri,numQuad);

    /// Parse elements globalIndex
    pugi::xml_node elements_original_indexes_array =
        pieceNode.child("CellData").find_child_by_attribute("DataArray","Name",
                "globalIndex");
    assert(!elements_original_indexes_array.empty());
    std::vector< globalIndex > allElementsOriginalIndexes;
    allElementsOriginalIndexes.reserve(numElements);
    SplitNodeTextString( elements_original_indexes_array.text().as_string(),
            allElementsOriginalIndexes,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (allElementsOriginalIndexes.size()) == numElements);

    /// Parse elements offsets
    std::vector< globalIndex> elementsOffsets(numElements);
    pugi::xml_node elementsOffsets_array =
        pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","offsets");
    assert(!elementsOffsets_array.empty());
    std::vector< globalIndex > allElementsOffsets(1,0); // convenient to have 0 as first offset
    allElementsOffsets.reserve(numElements);
    SplitNodeTextString( elementsOffsets_array.text().as_string(), allElementsOffsets,
            [](string str)-> globalIndex {return std::stoll(str);});
    assert(static_cast<globalIndex >(allElementsOffsets.size()) == numElements+1);

    /// Parse cells connectivities
    pugi::xml_node elementsConnectivityArray =
        pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","connectivity");
    assert(!elementsConnectivityArray.empty());
    std::vector< globalIndex > allElementsConnectivities;
    allElementsConnectivities.reserve(8 * numElements);
    SplitNodeTextString(elementsConnectivityArray.text().as_string(),
            allElementsConnectivities,
            [](string str)-> localIndex {return std::stoi(str);});
    for( globalIndex elementIndex = 0 ; elementIndex < numElements; ++elementIndex) {
        localIndex numVerticesInCell = allElementsOffsets[elementIndex+1] -
            allElementsOffsets[elementIndex];
        std::vector< globalIndex > connectivity(numVerticesInCell);
        for(localIndex co = 0 ; co < numVerticesInCell ; co++) {
            connectivity[co] =
                allElementsConnectivities[allElementsOffsets[elementIndex]+co];
        }
        if( allElementsTypes[elementIndex] == 10 ||
                allElementsTypes[elementIndex] == 14 ||
                allElementsTypes[elementIndex] ==12 ||
                allElementsTypes[elementIndex] == 13 ) {
            mesh.AddCell(connectivity);
            mesh.SetCellOriginalIndex(elementIndex,
                    allElementsOriginalIndexes[elementIndex]);
        }
        else if( allElementsTypes[elementIndex] == 5 ||
                allElementsTypes[elementIndex] == 9) {
            mesh.AddPolygon(connectivity);
            mesh.SetPolygonOriginalIndex(elementIndex,
                   allElementsOriginalIndexes[elementIndex]);
        }
        else {
            GEOS_ERROR("Element not recognised");
        }
    }

    mesh.Finish();
}



    ////////////////
    /// MESH PART //
    ////////////////

globalIndex DumbMesh::NumVertices() const {
    return m_numVertices;
}

globalIndex DumbMesh::NumCells() const {
    return m_numCells;
}

globalIndex DumbMesh::NumPolygons() const {
    return m_numPolygons;
}

globalIndex DumbMesh::NumTetra() const {
    return m_numTetra;
}

globalIndex DumbMesh::NumHex() const {
    return m_numHex;
}

globalIndex DumbMesh::NumPrism() const {
    return m_numPrism;
}

globalIndex DumbMesh::NumPyr() const {
    return m_numPyr;
}

globalIndex DumbMesh::NumTri() const {
    return m_numTri;
}

globalIndex DumbMesh::NumQuad() const {
    return m_numQuad;
}

globalIndex DumbMesh::TetIndexToCellIndex(globalIndex const tetIndex) {
    assert(tetIndex < m_numTetra);
    return m_tetIndexToCellIndex[tetIndex];
}

globalIndex DumbMesh::HexIndexToCellIndex(globalIndex const hexIndex) {
    assert(hexIndex < m_numHex);
    return m_hexIndexToCellIndex[hexIndex];
}

globalIndex DumbMesh::PrismIndexToCellIndex(globalIndex const prismIndex) {
    assert(prismIndex < m_numPrism);
    return m_prismIndexToCellIndex[prismIndex];
}

globalIndex DumbMesh::PyrIndexToCellIndex(globalIndex const pyrIndex) {
    assert(pyrIndex < m_numPyr);
    return m_pyrIndexToCellIndex[pyrIndex];
}

globalIndex DumbMesh::TriIndexToPolygonIndex(globalIndex const triIndex) {
    assert(triIndex < m_numTri);
    return m_triIndexToPolygonIndex[triIndex];
}

globalIndex DumbMesh::QuadIndexToPolygonIndex(globalIndex const quadIndex) {
    assert(quadIndex < m_numQuad);
    return m_pyrIndexToCellIndex[quadIndex];
}

localIndex DumbMesh::NumVerticesInCell(globalIndex const cellIndex ) const {
    assert(cellIndex < m_numCells);
    return m_cellsPtr[cellIndex+1] - m_cellsPtr[cellIndex];
}

localIndex DumbMesh::NumVerticesInPolygon(globalIndex const polygonIndex ) const {
    assert(polygonIndex < m_numPolygons);
    return m_polygonsPtr[polygonIndex+1] - m_polygonsPtr[polygonIndex];
}

globalIndex DumbMesh::CellVertexIndex(globalIndex const cellIndex,
        localIndex const local_corner_index) const {
    assert(local_corner_index < NumVerticesInCell(cellIndex));
    return m_cellsConnectivity[m_cellsPtr[cellIndex] + local_corner_index];
}

globalIndex DumbMesh::PolygonVertexIndex(globalIndex const polygonIndex,
        localIndex const local_corner_index) const {
    assert(local_corner_index < NumVerticesInPolygon(polygonIndex));
    return m_polygonsConnectivity[m_polygonsPtr[polygonIndex] + local_corner_index];
}

real64 const * DumbMesh::Vertex(globalIndex const vertexIndex) const {
    assert(vertexIndex < m_numVertices);
    return &(m_vertices[3*vertexIndex]);
}

globalIndex DumbMesh::GlobalPolygonIndex(globalIndex const polygonIndex) const {
    assert(polygonIndex<m_numPolygons);
    return m_globalPolygonIndexes[polygonIndex];
}

globalIndex DumbMesh::GlobalCellIndex(globalIndex const cellIndex) const {
    assert(cellIndex < m_numCells);
    return m_globalCellIndexes[cellIndex];
}

void DumbMesh::SetVertex(globalIndex const vertexIndex,
        std::vector< real64 > const& vertex) {
    assert(vertex.size() == 3); 
    assert(vertexIndex < m_numVertices);
    for( localIndex coor = 0 ; coor < 3 ; coor ++ ) {
        m_vertices[3*vertexIndex+coor] = vertex[coor];
    }
}

void DumbMesh::SetNumVertices(globalIndex const numVertices) {
    m_numVertices = numVertices;
    m_vertices.resize( numVertices*3);
}

void DumbMesh::SetNumCellAndPolygons(globalIndex numTetra,
                                       globalIndex numHex,
                                       globalIndex numPrism,
                                       globalIndex numPyr,
                                       globalIndex numTri,
                                       globalIndex numQuad) {
    m_numTetra = numTetra;
    m_numHex = numHex;
    m_numPrism = numPrism;
    m_numPyr = numPyr;
    m_numTri = numTri;
    m_numQuad = numQuad;
    m_numCells = numTetra + numHex + numPrism + numPyr;
    m_numPolygons = numTri + numQuad;
    m_cellsPtr.reserve(m_numCells+1);
    m_polygonsPtr.reserve( m_numPolygons +1);
    m_cellsConnectivity.reserve( 4*numTetra + 8* numHex + 6*numPrism + 5*numPyr );
    m_polygonsConnectivity.reserve( 3*numTri + 4*numQuad );
    m_globalPolygonIndexes.resize( m_numPolygons);
    m_globalCellIndexes.resize( m_numCells);
    m_tetIndexToCellIndex.reserve(m_numTetra);
    m_hexIndexToCellIndex.reserve(m_numHex);
    m_prismIndexToCellIndex.reserve(m_numPrism);
    m_pyrIndexToCellIndex.reserve(m_numPyr);
    m_triIndexToPolygonIndex.reserve(m_numTri);
    m_quadIndexToPolygonIndex.reserve(m_numQuad);
}

void DumbMesh::AddCell( std::vector<globalIndex> const & connectivity ) {
    m_cellsPtr.push_back( connectivity.size() + m_cellsPtr[m_cellsPtr.size()-1]);
    globalIndex numVerticesInCell = static_cast<localIndex>(connectivity.size());
    for( localIndex co = 0 ; co < numVerticesInCell; ++co ) {
        m_cellsConnectivity.push_back(connectivity[co]);
    }

    if( connectivity.size() == 4 ) {
        m_tetIndexToCellIndex.push_back(m_cellsPtr.size()-2);
    } else if (connectivity.size() == 8) {
        m_hexIndexToCellIndex.emplace_back(m_cellsPtr.size()-2);
    } else if (connectivity.size() == 6) {
        m_prismIndexToCellIndex.emplace_back(m_cellsPtr.size()-2);
    } else if (connectivity.size() == 5 ) {
        m_pyrIndexToCellIndex.emplace_back(m_cellsPtr.size()-2);
    }
}

void DumbMesh::AddPolygon( std::vector<globalIndex> const & connectivity ) {
    m_polygonsPtr.push_back( connectivity.size() + m_polygonsPtr[m_polygonsPtr.size()-1]);
    for( localIndex co = 0 ; co < static_cast<localIndex>(connectivity.size()); ++co ) {
        m_polygonsConnectivity.push_back(connectivity[co]);
    }

    if( connectivity.size() == 3 ) {
        m_triIndexToPolygonIndex.push_back(m_polygonsPtr.size()-2);
    } else if (connectivity.size() == 4) {
        m_quadIndexToPolygonIndex.emplace_back(m_polygonsPtr.size()-2);
    } 
}

void DumbMesh::SetCellOriginalIndex(globalIndex const cellIndexInPartMesh,
        globalIndex const cellIndexInFullMesh) {
    assert( cellIndexInPartMesh < static_cast< globalIndex>
            (m_globalCellIndexes.size() ) );
    m_globalCellIndexes[cellIndexInPartMesh] = cellIndexInFullMesh;
}

void DumbMesh::SetPolygonOriginalIndex(globalIndex const polygonIndexInPartMesh,
        globalIndex const polygonIndexInFullMesh) {
    assert( polygonIndexInPartMesh <
            static_cast< globalIndex> (m_globalPolygonIndexes.size() ) );
    m_globalPolygonIndexes[polygonIndexInPartMesh] = polygonIndexInFullMesh;
}   

void DumbMesh::SetName(string const & name) {
    m_name = name;
}

void DumbMesh::Finish() {
    assert(static_cast<globalIndex>(m_cellsPtr.size()) == m_numCells +1);
    assert(static_cast<globalIndex>(m_polygonsPtr.size()) == m_numPolygons +1);
    assert(static_cast<globalIndex>(m_tetIndexToCellIndex.size()) == m_numTetra);
    std::cout << "m_numHex " << m_numHex << std::endl;
    std::cout << "m_hexIndexToCellIndex.size()) " << m_hexIndexToCellIndex.size() <<std::endl;
    assert(static_cast<globalIndex>(m_hexIndexToCellIndex.size()) == m_numHex);
    assert(static_cast<globalIndex>(m_prismIndexToCellIndex.size()) == m_numPrism);
    assert(static_cast<globalIndex>(m_pyrIndexToCellIndex.size()) == m_numPyr);
    assert(static_cast<globalIndex>(m_triIndexToPolygonIndex.size()) == m_numTri);
    assert(static_cast<globalIndex>(m_quadIndexToPolygonIndex.size()) == m_numQuad);
}

string const & DumbMesh::Name() const {
    return m_name;
}

void RankBlock::TransferRankBlockToGEOSMesh( MeshLevel * const meshLevel ) const
{
  for (auto& meshBlock : m_block) {
      const auto & mesh = meshBlock.mesh();
      if( mesh.NumPolygons() == 0 ) {
          NodeManager * const nodeManager = meshLevel->getNodeManager();
          ElementRegionManager * const elemRegMananger = meshLevel->getElemManager();


          nodeManager->resize(mesh.NumVertices());
          arrayView1d<R1Tensor> X = nodeManager->referencePosition();
          std::cout << "node Manager pointer : " << nodeManager << std::endl;
          std::cout << "X pointer : " << X[0].Data() << std::endl;

          for( globalIndex a=0 ; a< mesh.NumVertices() ; ++a )
          {
              real64 * const tensorData = X[a].Data();
              tensorData[0] = mesh.Vertex(a)[0];
              tensorData[1] = mesh.Vertex(a)[1];
              tensorData[2] = mesh.Vertex(a)[2];
          }


          for( int i = 0 ; i < elemRegMananger->GetRegions().size() ; i++) {
              std::cout<< "REGION DBUT : " << i << std::endl;
              std::cout << elemRegMananger->GetRegion(i)->getName() << std::endl;;
              std::cout<< "REGION FIN : " << i << std::endl;
          }
          CellBlockSubRegion * const subRegion = elemRegMananger->GetRegion( "Region2" )->RegisterGroup<CellBlockSubRegion>("cb1");
          CellBlockSubRegion * const cellBlock =
              elemRegMananger->GetRegion( "Region2" )->GetSubRegion("HEX");
          std::cout << "HA" << std::endl;
          for(int i = 0 ; i < elemRegMananger->GetRegion( "Region2" )->GetSubRegions().size();i++){
              std::cout << "SUBREGION : " << i << std::endl;
          }
          std::cout << "BE" << std::endl;

          std::cout << "cellBlock pointer : " << cellBlock << std::endl;

          subRegion->resize( mesh.NumCells() );
          auto & cellToVertex = subRegion->nodeList();
          cellToVertex.resize(mesh.NumCells(),  mesh.NumVerticesInCell(0) );

          for( localIndex k=0 ; k<mesh.NumCells() ; ++k )
          {
              for( localIndex a=0 ; a< mesh.NumVerticesInCell(k) ; ++a )
              {
                  cellToVertex[k][a] = mesh.CellVertexIndex(k,a);
              }
          }
      }
  }
}

}
