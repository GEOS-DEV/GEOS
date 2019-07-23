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
 * @file PvtuFile.cpp
 */

#include <iostream>
#include <string.h>

#include "VtmFile.hpp"

#include "Logger.hpp"
#include "mesh/MeshBody.hpp"

#ifdef GEOSX_USE_MPI
#include <mpi.h>
#endif

namespace geosx{

    /*
void DumbPropertyManager::SetName(string const & name ) {
    m_name = name;
}

void DumbPropertyManager::SetSize(globalIndex const size) {
    m_values.resize(size);
}

void DumbPropertyManager::SetValue(globalIndex const elementIndex, double const value) {
    assert(elementIndex < static_cast<globalIndex>(m_values.size()));
    m_values[elementIndex] = value;
}

string DumbPropertyManager::GetName() const {
    return m_name;
}

double DumbPropertyManager::GetValue(globalIndex const elementIndex) const {
    assert(elementIndex < static_cast<globalIndex>(m_values.size()));
    return m_values[elementIndex];
}
*/

MeshBlock::MeshBlock( string fileName,
        string blockName ) :
    m_vtuFileName(fileName),
    m_blockName(blockName)
{
}

bool MeshBlock::IsARegionBlock() const  {
    return m_mesh.NumCells() > 0;
}

void MeshBlock::Load(bool loadMesh, bool loadProperties) {
    if( loadMesh ) {
    m_vtuFile.Load(m_vtuFileName,m_mesh);    
    m_mesh.SetName(m_blockName);
    }
    if( loadProperties ) {
    m_vtuFile.Load(m_vtuFileName,m_properties);    
    }
}

std::map< string, real64_array > const & MeshBlock::PropertyMap() const {
    return m_properties;
}

DumbMesh const & MeshBlock::mesh() const {
    return m_mesh;
}

void RankBlock::AddMeshBlock( const MeshBlock& block) {
    m_block.push_back(block);
}

void RankBlock::Load(bool loadMesh, bool loadProperties) {
    for( auto& block : m_block) {
        block.Load(loadMesh, loadProperties);
    }
}
localIndex RankBlock::NumMeshBlocks() const {
    return m_block.size();
}
MeshBlock const & RankBlock::GetMeshBlock(localIndex const meshBlockIndex) const {
    return m_block[ meshBlockIndex ];
}
// PUBLIC METHODS

RankBlock const & VtmFile::GetRankBlock(localIndex const rankBlockIndex) const {
    return m_rankBlocks[rankBlockIndex];
}

localIndex VtmFile::NumRankBlocks() const{
    return m_rankBlocks.size();
}
void VtmFile::Load( string const &fileName, bool loadMesh, bool loadProperties) {
    int mpiSize = 0;
    int mpiRank = 0;
    m_fileName = fileName;
#ifdef GEOSX_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
    int numFilesPerRank = 0;
    //int nb_total_partitions = 0;
    int remainderFiles = 0;
    // XML file parsing using pugixml tool
    pugi::xml_document vtmDoc;
    vtmDoc.load_file(fileName.c_str());

    if( mpiRank == 0) {
        CheckXmlFileConsistency(vtmDoc,fileName);
    }
    array1d< RankBlock > rankBlocks;
    SetRanksAndBlocks(vtmDoc,rankBlocks);
    // Retrieve the number of partitions
    int const numFiles = integer_conversion<int>(rankBlocks.size());

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
            m_rankBlocks.push_back(rankBlocks[mpiRank*(numFilesPerRank+1) + p_index]);
        }
    } else if( mpiRank < numFiles) {
        for(int p_index = 0; p_index < numFilesPerRank; ++p_index) {
            m_rankBlocks.push_back(rankBlocks[remainderFiles+mpiRank*(numFilesPerRank) + p_index]);
        }
    }
    if( mpiRank < numFiles ) {
        for(localIndex p_index = 0; p_index < 
                m_rankBlocks.size(); ++p_index) {
            m_rankBlocks[p_index].Load(loadMesh, loadProperties);
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

void VtuFile::Load( string const &filename,
        std::map< string, real64_array > & properties) {
    pugi::xml_document vtmDoc;
    vtmDoc.load_file(filename.c_str());
    CheckXmlChildFileConsistency(vtmDoc,filename);
    LoadProperties(vtmDoc, properties);
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

    localIndex countNbRank=0;
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
                array1d< RankBlock >& rankBlocks) {
    auto const & vtkFileNode =vtmDoc.child("VTKFile");
    auto const & vtkMultiBlockDataSetNode =  vtkFileNode.child("vtkMultiBlockDataSet");
    string dirPath;
    if( m_fileName.rfind("/") != std::string::npos) {
        dirPath = m_fileName.substr(0, m_fileName.rfind("/"));
    }
    std::cout << "dirpath " << dirPath << std::endl;
    
    for(auto const & rank : vtkMultiBlockDataSetNode.children()) {
        if( rank.name() == static_cast< string > ("Block") )
        {
            RankBlock curRankBlock;
            for( auto const & block : rank.children() ) {
                if( block.name() == static_cast< string > ("DataSet")) {
                    string filename;
                    if( dirPath.empty() ) {
                        filename = block.attribute("file").as_string();
                    }
                    else {
                        filename = dirPath+"/"+block.attribute("file").as_string();
                    }
                    std::cout << filename << std::endl;
                    MeshBlock curMeshBlock(
                            filename,
                            block.attribute("name").as_string());
                    curRankBlock.AddMeshBlock(curMeshBlock);
                }
            }
            rankBlocks.push_back(curRankBlock);
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
        }
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
        //GEOS_ERROR("Can't find any DataArray which contains the property originalIndex in  CellData in "+filename);
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
        array1d< T >& out,
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

    localIndex numVertices = pieceNode.attribute("NumberOfPoints").as_int();
    mesh.SetNumVertices( numVertices );

    localIndex numElements = pieceNode.attribute("NumberOfCells").as_int();

    /// Parse vertices
    pugi::xml_node m_verticesarray =
        pieceNode.child("Points").find_child_by_attribute("DataArray","Name","Points");
    assert(!m_verticesarray.empty()) ;
    real64_array allVertices;
    allVertices.reserve( numVertices*3 );
    SplitNodeTextString(m_verticesarray.text().as_string(), allVertices, 
            [](string str)-> double {return std::stod(str);});
    assert(allVertices.size() / 3 == numVertices);

    /// Fill the vertices in the mesh
    for(localIndex vertexIndex  = 0 ; vertexIndex < numVertices; ++vertexIndex) {
        real64_array vertex(3);
        for(localIndex coor = 0 ; coor < 3 ; ++coor) {
            vertex[coor] = allVertices[3*vertexIndex+coor];
        }
        mesh.SetVertex(vertexIndex, vertex);
    }

    /// Parse elements types
    pugi::xml_node elementsTypesArray =
        pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","types");
    assert(!elementsTypesArray.empty());
    localIndex_array allElementsTypes;
    allElementsTypes.reserve(numElements);
    SplitNodeTextString( elementsTypesArray.text().as_string(), allElementsTypes,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(allElementsTypes.size() == numElements);
    localIndex numTetra = 0;
    localIndex numHex = 0;
    localIndex numPrism = 0;
    localIndex numPyr = 0;
    localIndex numTri = 0;
    localIndex numQuad = 0;
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
    /*
    pugi::xml_node elements_original_indexes_array =
        pieceNode.child("CellData").find_child_by_attribute("DataArray","Name",
                "globalIndex");
    assert(!elements_original_indexes_array.empty());
    globalIndex_array allElementsOriginalIndexes;
    allElementsOriginalIndexes.reserve(numElements);
    SplitNodeTextString( elements_original_indexes_array.text().as_string(),
            allElementsOriginalIndexes,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (allElementsOriginalIndexes.size()) == numElements);
    */

    /// Parse elements offsets
    localIndex_array elementsOffsets(numElements);
    pugi::xml_node elementsOffsets_array =
        pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","offsets");
    assert(!elementsOffsets_array.empty());
    localIndex_array allElementsOffsets(1);
    allElementsOffsets[0]=1; // convenient to have 0 as first offset
    allElementsOffsets.reserve(numElements);
    SplitNodeTextString( elementsOffsets_array.text().as_string(), allElementsOffsets,
            [](string str)-> localIndex {return std::stoi(str);});
    assert(allElementsOffsets.size() == numElements+1);

    /// Parse cells connectivities
    pugi::xml_node elementsConnectivityArray =
        pieceNode.child("Cells").find_child_by_attribute("DataArray","Name","connectivity");
    assert(!elementsConnectivityArray.empty());
    localIndex_array allElementsConnectivities;
    allElementsConnectivities.reserve(8 * numElements);
    SplitNodeTextString(elementsConnectivityArray.text().as_string(),
            allElementsConnectivities,
            [](string str)-> localIndex {return std::stoi(str);});
    for( localIndex elementIndex = 0 ; elementIndex < numElements; ++elementIndex) {
        localIndex numVerticesInCell = allElementsOffsets[elementIndex+1] -
            allElementsOffsets[elementIndex];
        localIndex_array connectivity(numVerticesInCell);
        for(localIndex co = 0 ; co < numVerticesInCell ; co++) {
            connectivity[co] =
                allElementsConnectivities[allElementsOffsets[elementIndex]+co];
        }
        if( allElementsTypes[elementIndex] == 10 ||
                allElementsTypes[elementIndex] == 14 ||
                allElementsTypes[elementIndex] ==12 ||
                allElementsTypes[elementIndex] == 13 ) {
            mesh.AddCell(connectivity);
            /*
            mesh.SetCellOriginalIndex(elementIndex,
                    allElementsOriginalIndexes[elementIndex]);
                    */
        }
        else if( allElementsTypes[elementIndex] == 5 ||
                allElementsTypes[elementIndex] == 9) {
            mesh.AddPolygon(connectivity);
            /*
            mesh.SetPolygonOriginalIndex(elementIndex,
                   allElementsOriginalIndexes[elementIndex]);
                   */
        }
        else {
            GEOS_ERROR("Element not recognised");
        }
    }

    mesh.Finish();
}


void VtuFile::LoadProperties(pugi::xml_document const & vtmDoc,
        std::map< string, real64_array >& properties){
    pugi::xml_node cellDataNode =
        vtmDoc.child("VTKFile").child("UnstructuredGrid").child("Piece").child("CellData");
    for( auto& property : cellDataNode.children()) {
        std::string propertyName = property.attribute("Name").as_string();
        std::string propertyType = property.attribute("type").as_string();
        if(propertyName != "globalIndex") {
            if(propertyType == "Float64" || propertyType == "Float32") {
                real64_array propertyToBeInserted;
                properties.insert(std::pair< string, real64_array >(propertyName,propertyToBeInserted));
                auto& curProperty = properties.find(propertyName)->second;
                real64_array propertyValues; // convenient to have 0 as first offset
                SplitNodeTextString( property.text().as_string(), propertyValues,
                        [](string str)-> double {return std::stod(str);});
                curProperty.resize(propertyValues.size());
                for( localIndex elementIndex = 0 ;
                        elementIndex < propertyValues.size(); elementIndex++) {
                    curProperty[elementIndex] = propertyValues[elementIndex];
                }
            }
        }
    }
}

    ////////////////
    /// MESH PART //
    ////////////////

localIndex DumbMesh::NumVertices() const {
    return m_numVertices;
}

localIndex DumbMesh::NumCells() const {
    return m_numCells;
}

localIndex DumbMesh::NumPolygons() const {
    return m_numPolygons;
}

localIndex DumbMesh::NumTetra() const {
    return m_numTetra;
}

localIndex DumbMesh::NumHex() const {
    return m_numHex;
}

localIndex DumbMesh::NumPrism() const {
    return m_numPrism;
}

localIndex DumbMesh::NumPyr() const {
    return m_numPyr;
}

localIndex DumbMesh::NumTri() const {
    return m_numTri;
}

localIndex DumbMesh::NumQuad() const {
    return m_numQuad;
}

localIndex DumbMesh::TetIndexToCellIndex(localIndex const tetIndex) {
    assert(tetIndex < m_numTetra);
    return m_tetIndexToCellIndex[tetIndex];
}

localIndex DumbMesh::HexIndexToCellIndex(localIndex const hexIndex) {
    assert(hexIndex < m_numHex);
    return m_hexIndexToCellIndex[hexIndex];
}

localIndex DumbMesh::PrismIndexToCellIndex(localIndex const prismIndex) {
    assert(prismIndex < m_numPrism);
    return m_prismIndexToCellIndex[prismIndex];
}

localIndex DumbMesh::PyrIndexToCellIndex(localIndex const pyrIndex) {
    assert(pyrIndex < m_numPyr);
    return m_pyrIndexToCellIndex[pyrIndex];
}

localIndex DumbMesh::TriIndexToPolygonIndex(localIndex const triIndex) {
    assert(triIndex < m_numTri);
    return m_triIndexToPolygonIndex[triIndex];
}

localIndex DumbMesh::QuadIndexToPolygonIndex(localIndex const quadIndex) {
    assert(quadIndex < m_numQuad);
    return m_pyrIndexToCellIndex[quadIndex];
}

localIndex DumbMesh::NumVerticesInCell(localIndex const cellIndex ) const {
    assert(cellIndex < m_numCells);
    return m_cellsPtr[cellIndex+1] - m_cellsPtr[cellIndex];
}

localIndex DumbMesh::NumVerticesInPolygon(localIndex const polygonIndex ) const {
    assert(polygonIndex < m_numPolygons);
    return m_polygonsPtr[polygonIndex+1] - m_polygonsPtr[polygonIndex];
}

localIndex DumbMesh::CellVertexIndex(localIndex const cellIndex,
        localIndex const local_corner_index) const {
    assert(local_corner_index < NumVerticesInCell(cellIndex));
    return m_cellsConnectivity[m_cellsPtr[cellIndex] + local_corner_index];
}

localIndex DumbMesh::PolygonVertexIndex(localIndex const polygonIndex,
        localIndex const local_corner_index) const {
    assert(local_corner_index < NumVerticesInPolygon(polygonIndex));
    return m_polygonsConnectivity[m_polygonsPtr[polygonIndex] + local_corner_index];
}

real64 const * DumbMesh::Vertex(localIndex const vertexIndex) const {
    assert(vertexIndex < m_numVertices);
    return &(m_vertices[3*vertexIndex]);
}

globalIndex DumbMesh::GlobalPolygonIndex(localIndex const polygonIndex) const {
    assert(polygonIndex<m_numPolygons);
    return m_globalPolygonIndexes[polygonIndex];
}

globalIndex DumbMesh::GlobalCellIndex(localIndex const cellIndex) const {
    assert(cellIndex < m_numCells);
    return m_globalCellIndexes[cellIndex];
}

void DumbMesh::SetVertex(localIndex const vertexIndex,
        real64_array const & vertex) {
    assert(vertex.size() == 3); 
    assert(vertexIndex < m_numVertices);
    for( localIndex coor = 0 ; coor < 3 ; coor ++ ) {
        m_vertices[3*vertexIndex+coor] = vertex[coor];
    }
}

void DumbMesh::SetNumVertices(localIndex const numVertices) {
    m_numVertices = numVertices;
    m_vertices.resize( numVertices*3);
}

void DumbMesh::SetNumCellAndPolygons(localIndex numTetra,
                                       localIndex numHex,
                                       localIndex numPrism,
                                       localIndex numPyr,
                                       localIndex numTri,
                                       localIndex numQuad) {
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

void DumbMesh::AddCell( localIndex_array const & connectivity ) {
    m_cellsPtr.push_back( connectivity.size() + m_cellsPtr[m_cellsPtr.size()-1]);
    localIndex numVerticesInCell = connectivity.size();
    for( localIndex co = 0 ; co < numVerticesInCell; ++co ) {
        m_cellsConnectivity.push_back(connectivity[co]);
    }

    if( connectivity.size() == 4 ) {
        m_tetIndexToCellIndex.push_back(m_cellsPtr.size()-2);
    } else if (connectivity.size() == 8) {
        m_hexIndexToCellIndex.push_back(m_cellsPtr.size()-2);
    } else if (connectivity.size() == 6) {
        m_prismIndexToCellIndex.push_back(m_cellsPtr.size()-2);
    } else if (connectivity.size() == 5 ) {
        m_pyrIndexToCellIndex.push_back(m_cellsPtr.size()-2);
    }
}

void DumbMesh::AddPolygon( localIndex_array const & connectivity ) {
    m_polygonsPtr.push_back( connectivity.size() + m_polygonsPtr[m_polygonsPtr.size()-1]);
    for( localIndex co = 0 ; co < connectivity.size(); ++co ) {
        m_polygonsConnectivity.push_back(connectivity[co]);
    }

    if( connectivity.size() == 3 ) {
        m_triIndexToPolygonIndex.push_back(m_polygonsPtr.size()-2);
    } else if (connectivity.size() == 4) {
        m_quadIndexToPolygonIndex.push_back(m_polygonsPtr.size()-2);
    } 
}

void DumbMesh::SetCellOriginalIndex(localIndex const cellIndexInPartMesh,
        globalIndex const cellIndexInFullMesh) {
    assert( cellIndexInPartMesh < 
            m_globalCellIndexes.size() );
    m_globalCellIndexes[cellIndexInPartMesh] = cellIndexInFullMesh;
}

void DumbMesh::SetPolygonOriginalIndex(localIndex const polygonIndexInPartMesh,
        globalIndex const polygonIndexInFullMesh) {
    assert( polygonIndexInPartMesh <
            m_globalPolygonIndexes.size() );
    m_globalPolygonIndexes[polygonIndexInPartMesh] = polygonIndexInFullMesh;
}   

void DumbMesh::SetName(string const & name) {
    m_name = name;
}

void DumbMesh::Finish() {
    assert(m_cellsPtr.size() == m_numCells +1);
    assert(m_polygonsPtr.size() == m_numPolygons +1);
    assert(m_tetIndexToCellIndex.size() == m_numTetra);
    assert(m_hexIndexToCellIndex.size() == m_numHex);
    assert(m_prismIndexToCellIndex.size() == m_numPrism);
    assert(m_pyrIndexToCellIndex.size() == m_numPyr);
    assert(m_triIndexToPolygonIndex.size() == m_numTri);
    assert(m_quadIndexToPolygonIndex.size() == m_numQuad);
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
          ElementRegionManager * const elemRegManager = meshLevel->getElemManager();


          nodeManager->resize(mesh.NumVertices());
          arrayView1d<R1Tensor> & X = nodeManager->referencePosition();
          std::cout << "node Manager pointer : " << nodeManager << std::endl;
          std::cout << "X pointer : " << X[0].Data() << std::endl;

          for( localIndex a=0 ; a< mesh.NumVertices() ; ++a )
          {
              real64 * const tensorData = X[a].Data();
              tensorData[0] = mesh.Vertex(a)[0];
              tensorData[1] = mesh.Vertex(a)[1];
              tensorData[2] = mesh.Vertex(a)[2];
          }


          for( int i = 0 ; i < elemRegManager->GetRegions().size() ; i++) {
              std::cout<< "REGION DBUT : " << i << std::endl;
              std::cout << elemRegManager->GetRegion(i)->getName() << std::endl;
              std::cout<< "REGION FIN : " << i << std::endl;
          }
          CellElementSubRegion * const subRegion = elemRegManager->GetRegion( "Region2" )->RegisterGroup<CellElementSubRegion>("cb1");
          CellElementSubRegion * const cellBlock =
              elemRegManager->GetRegion( "Region2" )->GetSubRegion<CellElementSubRegion>("HEX");
          std::cout << "HA" << std::endl;
          for(int i = 0 ; i < elemRegManager->GetRegion( "Region2" )->GetSubRegions().size();i++){
              std::cout << "SUBREGION : " << i << std::endl;
          }
          std::cout << "BE" << std::endl;

          std::cout << "cellBlock pointer : " << cellBlock << std::endl;

          subRegion->resize( mesh.NumCells());
          auto & cellToVertex = subRegion->nodeList();
          cellToVertex.resize(mesh.NumCells(),mesh.NumVerticesInCell(0));

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
