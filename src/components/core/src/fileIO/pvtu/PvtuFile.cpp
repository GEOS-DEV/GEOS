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
void PvtuFile::load( std::string const &filename) {
    int mpiSize = 0;
    int mpiRank = 0;
#if USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
    int numFilesPerRank = 0;
    //int nb_total_partitions = 0;
    int remainderFiles = 0;
    std::vector < std::string > children_files;



    // XML file parsing using pugixml tool
    pugi::xml_document pvtu_doc;
    pvtu_doc.load_file(filename.c_str());

    if( mpiRank == 0) {
        check_xml_parent_file_consistency(pvtu_doc,filename);
    }
    vtu_files_list(pvtu_doc,children_files);
    // Retrieve the number of partitions
    int const numFiles = children_files.size();

    // Next part of this method is dedicated to the optimization of file loading
    //
    // IF nb_partitions > nb_mpi_process : Each process will load 
    // nb_partitions / nb_mpi_process. Processes with a rank < nb_partitions % nb_mpi_process
    // will load an additional partition
    //
    // IF nb_partitions < nb_mpi_process : the first nb_partitions process will
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
        vtu_file_names_.resize(numFilesPerRank+1);
        vtu_files_.resize(numFilesPerRank+1);
        for(int p_index = 0; p_index < numFilesPerRank +1; ++p_index) {
            vtu_file_names_[p_index] =
                children_files[mpiRank*(numFilesPerRank+1) + p_index];
        }
    } else if( mpiRank < numFiles) {
        vtu_file_names_.resize(numFilesPerRank);
        vtu_files_.resize(numFilesPerRank);
        for(int p_index = 0; p_index < numFilesPerRank; ++p_index) {
            vtu_file_names_[p_index] =
                children_files[mpiRank*(numFilesPerRank) + p_index+remainderFiles];
        }
    }
    if( mpiRank < numFiles ) {
        for(int p_index = 0; p_index < 
                static_cast< int >(vtu_file_names_.size()); ++p_index) {
            vtu_files_[p_index].load(vtu_file_names_[p_index]);
        }
    }
}

void PvtuFile::save( std::string const &filename) {
    geos_abort("pvtu file save is not implemented yet");
}

void VtuFile::load( std::string const &filename) {
    pugi::xml_document pvtu_doc;
    pvtu_doc.load_file(filename.c_str());
    check_xml_child_file_consistency(pvtu_doc,filename);
    load_mesh_part(pvtu_doc);
}

void VtuFile::save( std::string const &filename) {
    geos_abort("vtu file save is not implemented yet");
}


//PRIVATE METHODS
void PvtuFile::check_xml_file_consistency(pugi::xml_document const & pvtu_doc,
        std::string const & filename,
        std::string const & prefix) const {

    // VTKFile is the main node of the pvtufile
    auto const vtk_file =pvtu_doc.child("VTKFile");
    if( vtk_file.empty() ) {
        geos_abort("Main node VTKFile not found in " + filename);
    }

    std::string ugrid_name = prefix +"UnstructuredGrid";
    auto const ugrid = vtk_file.child(ugrid_name.c_str());
    if( ugrid.empty() ) {
        geos_abort("Node " + prefix + "UnstructuredGrid not found or empty in " +
                filename);
    }
    pugi::xml_node main_node;
    if( prefix == "P" ) {
        main_node = ugrid;
    } else {
        main_node = ugrid.child("Piece");
        if (main_node.empty()) {
            geos_abort("Piece node is missing in " + filename);
        }
    }

    std::string const point_data_name = prefix + "PointData";
    std::string const cell_data_name = prefix + "CellData";
    std::string const points_name = prefix + "Points";
    std::string const main_node_child_names[3] = {point_data_name,
        cell_data_name,points_name};
    for( auto main_node_child_name : main_node_child_names ) {
        auto const main_node_child = main_node.child(main_node_child_name.c_str());
        if(main_node_child.empty()) {
            geos_abort("Node " +main_node_child_name +
                " not found or empty in " + filename);
        };
    }

    bool points_have_original_index_property = false;
    std::string const mandatory_attributes[3] =
    {"Name","type","format"};
    for( auto pdata_property : main_node.child(point_data_name.c_str()).children() ) {
        if( pdata_property.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    pdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    geos_abort("Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of " + prefix + "PointData");
                }
            }
            if( pdata_property.attribute("Name").value() == str_original_index_) {
                points_have_original_index_property = true;
            }
        }
    }
    if (!points_have_original_index_property) {
        geos_abort("Can't find any DataArray which contains the property " +
            str_original_index_ + " in " + prefix + "PointData in " + filename);
    }

    bool cells_have_original_index_property = false;
    bool cells_have_partition_property = false;
    bool cells_have_region_property = false;
    for( auto pdata_property : main_node.child(cell_data_name.c_str()).children() ) {
        if( pdata_property.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    pdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    geos_abort("Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of " + prefix + "CellData in "
                        + filename);
                }
            }
            if( pdata_property.attribute("Name").value() == str_original_index_) {
                cells_have_original_index_property = true;
            }
            if( pdata_property.attribute("Name").value() == str_partition_) {
                cells_have_partition_property = true;
            }
            if( pdata_property.attribute("Name").value() == str_region_) {
                cells_have_region_property = true;
            }
        }
    }
    if (!cells_have_original_index_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_original_index_ + " in " + prefix + "CellData in "+filename);
    }
    if (!cells_have_region_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_region_ + " in "+ prefix + "CellData in " + filename);
    }
    if (!cells_have_partition_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_partition_ + " in " + prefix + "CellData in " + filename);
    }

    bool point_has_a_pdata_array = false;
    bool point_has_a_pdata_array_with_points = false;
    std::string const data_array_name = prefix + "DataArray";
    for(auto point_child : main_node.child(points_name.c_str()).children() ) {
        if( point_child.name() == data_array_name ) {
            point_has_a_pdata_array = true;
            if( point_child.attribute("Name").as_string() ==
                    static_cast< std::string >("Points")) {
                point_has_a_pdata_array_with_points = true;
                if (point_child.attribute("NumberOfComponents").as_uint() != 3 ) {
                    geos_abort("GEOSX supports only 3D meshes");
                }
                break;
            }
        }
    }
    if( !point_has_a_pdata_array ) {
        geos_abort("No " + prefix +"DataArray found in " + prefix + "Points in " + filename);
    }
    if( !point_has_a_pdata_array_with_points ) {
        geos_abort("No " + prefix +"DataArray named \"Points\" found in " + filename);
    }

}
void PvtuFile::check_xml_parent_file_consistency(pugi::xml_document const & pvtu_doc,
        std::string const & filename) const {

        check_xml_file_consistency( pvtu_doc, filename, "P");

        bool has_a_piece = false;
        for( auto ugrid_child :
                pvtu_doc.child("VTKFile").child("PUnstructuredGrid").children() ) {
            if(ugrid_child.name() == static_cast< std::string >("Piece")) {
                has_a_piece = true;
                if( ugrid_child.attribute("Source").empty() ) {
                    geos_abort("Piece nodes has to have an attribute Source not empty.");
                }
            }
        }
        if( !has_a_piece ) {
            geos_abort("No Piece not found in " + filename );
        }
}
void PvtuFile::vtu_files_list(
        pugi::xml_document const & pvtu_doc,
        std::vector < std::string > & vtu_files ) const{
    int nb_partitions = 0;
    for(auto child : pvtu_doc.child("VTKFile").child("PUnstructuredGrid").children()) {
        if( child.name() == static_cast< std::string > ("Piece") )
        {
            vtu_files.emplace_back(child.attribute("Source").as_string());
        }
    }
}

void VtuFile::check_xml_child_file_consistency(pugi::xml_document const & pvtu_doc,
        std::string const & filename) const {
    check_xml_file_consistency( pvtu_doc, filename);

    pugi::xml_node piece_node =
        pvtu_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");
    if( piece_node.attribute("NumberOfPoints").empty() ) {
        geos_abort("Attribute \"NumberOfPoints\" of Node \"Piece\" is missing or empty in "
                + filename);
    }

    if( piece_node.attribute("NumberOfCells").empty() ) {
        geos_abort("Attribute \"NumberOfCells\" of Node \"Piece\" is missing or empty in "
                + filename);
    }

    pugi::xml_node cell_node = piece_node.child("Cells");
    std::string const mandatory_attributes[3] =
    {"Name","type","format"};
    bool cells_have_connectivity = false;
    bool cells_have_type = false;
    bool cells_have_offset = false;
    for( auto cell_att : cell_node.children() ) {
        if( cell_att.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    cell_att.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    geos_abort("Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of in "
                        + filename);
                }
            }
            if( cell_att.attribute("Name").value() ==
                    static_cast<std::string>("connectivity")) {
                cells_have_connectivity = true;
            }
            if( cell_att.attribute("Name").value() == static_cast< std::string >("offsets")) {
                cells_have_offset = true;
            }
            if( cell_att.attribute("Name").value() == static_cast< std::string > ("types")) {
                cells_have_type = true;
            }
        }
    }
    if (!cells_have_connectivity) {
        geos_abort("No property \"connectivity\" found in \"Cells\" node of " + filename);
    }
    if (!cells_have_type) {
        geos_abort("No property \"types\" found in \"Cells\" node of " + filename);
    }
    if (!cells_have_offset) {
        geos_abort("No property \"offsets\" found in \"Cells\" node of " + filename);
    }
}   

template<typename T, typename Lambda>
void VtuFile::split_node_text_string( std::string const & in,
        std::vector< T >& out,
        Lambda && string_convertor) const {
    std::stringstream stringStream(in);
    std::string line;
    while(std::getline(stringStream, line))
    {
        std::size_t prev = 0, pos;
        while ((pos = line.find_first_of(" ", prev)) != std::string::npos)
        {
            if (pos > prev)
                out.push_back(string_convertor(line.substr(prev, pos-prev)));
            prev = pos+1;
        }
        if (prev < line.length())
            out.push_back(string_convertor(line.substr(prev, std::string::npos)));
    }
}

void VtuFile::load_mesh_part(pugi::xml_document const & pvtu_doc){
    pugi::xml_node piece_node =
        pvtu_doc.child("VTKFile").child("UnstructuredGrid").child("Piece");

    globalIndex nb_vertices = piece_node.attribute("NumberOfPoints").as_llong();
    mesh_part_.set_nb_vertices( nb_vertices );

    globalIndex nb_elements = piece_node.attribute("NumberOfCells").as_llong();
    mesh_part_.reserve_nb_cells_and_polygons( nb_elements );

    /// Parse vertices
    pugi::xml_node vertices_array =
        piece_node.child("Points").find_child_by_attribute("DataArray","Name","Points");
    assert(!vertices_array.empty()) ;
    std::vector< real64 > all_vertices;
    all_vertices.reserve( nb_vertices*3 );
    split_node_text_string(vertices_array.text().as_string(), all_vertices, 
            [](std::string str)-> double {return std::stod(str);});
    assert(static_cast< globalIndex> (all_vertices.size()) / 3 == nb_vertices);

    /// Parse vertices original index
    pugi::xml_node vertices_original_indexes_array = 
        piece_node.child("PointData").find_child_by_attribute(
                "DataArray","Name","original_index");
    assert(!vertices_original_indexes_array.empty());
    std::vector< globalIndex > all_vertices_original_indexes;
    all_vertices_original_indexes.reserve( nb_vertices );
    split_node_text_string(vertices_original_indexes_array.text().as_string(),
            all_vertices_original_indexes,
            [](std::string str)-> globalIndex {return std::stoll(str);});
    assert( nb_vertices == static_cast< globalIndex> (all_vertices_original_indexes.size() ));

    /// Fill the vertices in the mesh
    for(globalIndex vertex_index  = 0 ; vertex_index < nb_vertices; ++vertex_index) {
        std::vector< real64 > vertex(3);
        for(localIndex coor = 0 ; coor < 3 ; ++coor) {
            vertex[coor] = all_vertices[3*vertex_index+coor];
        }
        mesh_part_.set_vertex_original_index(vertex_index,
                all_vertices_original_indexes[vertex_index]);
        mesh_part_.set_vertex(vertex_index, vertex);
    }

    /// Parse elements types
    pugi::xml_node elements_types_array =
        piece_node.child("Cells").find_child_by_attribute("DataArray","Name","types");
    assert(!elements_types_array.empty());
    std::vector< localIndex > all_elements_types;
    all_elements_types.reserve(nb_elements);
    split_node_text_string( elements_types_array.text().as_string(), all_elements_types,
            [](std::string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (all_elements_types.size()) == nb_elements);

    /// Parse elements regions
    pugi::xml_node elements_regions_array =
        piece_node.child("CellData").find_child_by_attribute("DataArray","Name","region");
    assert(!elements_regions_array.empty());
    std::vector< globalIndex > all_elements_regions;
    all_elements_regions.reserve(nb_elements);
    split_node_text_string( elements_regions_array.text().as_string(), all_elements_regions,
            [](std::string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (all_elements_regions.size()) == nb_elements);
    
    /// Parse elements original_index
    pugi::xml_node elements_original_indexes_array =
        piece_node.child("CellData").find_child_by_attribute("DataArray","Name",
                "original_index");
    assert(!elements_original_indexes_array.empty());
    std::vector< globalIndex > all_elements_original_indexes;
    all_elements_original_indexes.reserve(nb_elements);
    split_node_text_string( elements_original_indexes_array.text().as_string(),
            all_elements_original_indexes,
            [](std::string str)-> localIndex {return std::stoi(str);});
    assert(static_cast< globalIndex> (all_elements_original_indexes.size()) == nb_elements);

    /// Parse elements offsets
    std::vector< globalIndex> elements_offsets(nb_elements);
    pugi::xml_node elements_offsets_array =
        piece_node.child("Cells").find_child_by_attribute("DataArray","Name","offsets");
    assert(!elements_offsets_array.empty());
    std::vector< globalIndex > all_elements_offsets(1,0); // convenient to have 0 as first offset
    all_elements_offsets.reserve(nb_elements);
    split_node_text_string( elements_offsets_array.text().as_string(), all_elements_offsets,
            [](std::string str)-> globalIndex {return std::stoll(str);});
    assert(static_cast<globalIndex >(all_elements_offsets.size()) == nb_elements+1);

    /// Parce cells connectivities
    pugi::xml_node elements_connectivity_array =
        piece_node.child("Cells").find_child_by_attribute("DataArray","Name","connectivity");
    assert(!elements_connectivity_array.empty());
    std::vector< globalIndex > all_elements_connectivities;
    all_elements_connectivities.reserve(8 * nb_elements);
    split_node_text_string(elements_connectivity_array.text().as_string(),
            all_elements_connectivities,
            [](std::string str)-> localIndex {return std::stoi(str);});
    for( globalIndex element_index = 0 ; element_index < nb_elements; ++element_index) {
        localIndex nb_vertices_in_cell = all_elements_offsets[element_index+1] -
            all_elements_offsets[element_index];
        std::vector< globalIndex > connectivity(nb_vertices_in_cell);
        for(localIndex co = 0 ; co < nb_vertices_in_cell ; co++) {
            connectivity[co] =
                all_elements_connectivities[all_elements_offsets[element_index]+co];
        }
        if( all_elements_types[element_index] == 10 ||
                all_elements_types[element_index] == 14 ||
                all_elements_types[element_index] ==12 ||
                all_elements_types[element_index] == 13 ) {
            mesh_part_.add_cell(connectivity);
            mesh_part_.set_cell_region(element_index, all_elements_regions[element_index]);
            mesh_part_.set_cell_original_index(element_index,
                    all_elements_original_indexes[element_index]);
        }
        else if( all_elements_types[element_index] == 5 ||
                all_elements_types[element_index] == 9) {
            mesh_part_.add_polygon(connectivity);
            mesh_part_.set_polygon_surface(element_index,
                    all_elements_regions[element_index]);
            mesh_part_.set_polygon_original_index(element_index,
                    all_elements_original_indexes[element_index]);
        }
    }

    mesh_part_.finish();
}

void VtuFile::TransferMeshPartToGEOSMesh( MeshLevel * const meshLevel )
{
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  ElementRegionManager * const elemRegMananger = meshLevel->getElemManager();

  arrayView1d<R1Tensor> X = nodeManager->referencePosition();

  nodeManager->resize(mesh_part_.nb_vertices());

  real64 const * const vertexData = mesh_part_.vertices_.data();

  for( localIndex a=0 ; a<mesh_part_.nb_vertices() ; ++a )
  {
    real64 * const tensorData = X[a].Data();
    tensorData[0] = vertexData[3*a];
    tensorData[1] = vertexData[3*a+1];
    tensorData[2] = vertexData[3*a+2];
  }

  string regionName; // = mesh_part.regionName();
  CellBlockSubRegion * const cellBlock = elemRegMananger->GetRegion( regionName )->GetSubRegion(0);
  lArray2d & cellToVertex = cellBlock->nodeList();
  cellToVertex.resize( 0, mesh_part_.nb_vertices_in_cell(0) );
  cellBlock->resize( mesh_part_.nb_cells() );

  for( localIndex k=0 ; k<mesh_part_.nb_cells() ; ++k )
  {
    for( localIndex a=0 ; a<mesh_part_.nb_vertices_in_cell(0) ; ++a )
    {
      cellToVertex[k][a] = mesh_part_.cell_vertex_index(k,a);
    }
  }



}


    ////////////////
    /// MESH PART //
    ////////////////

    globalIndex MeshPart::nb_vertices() const {
        return nb_vertices_;
    }

    globalIndex MeshPart::nb_cells() const {
        return nb_cells_;
    }

    globalIndex MeshPart::nb_polygons() const {
        return nb_polygons_;
    }

    localIndex MeshPart::nb_vertices_in_cell( const globalIndex cell_index ) const {
        assert(cell_index < nb_cells_);
        return cells_ptr_[cell_index+1] - cells_ptr_[cell_index-1];
    }

    localIndex MeshPart::nb_vertices_in_polygon( const globalIndex polygon_index ) const {
        assert(polygon_index < nb_polygons_);
        return polygons_ptr_[polygon_index+1] - polygons_ptr_[polygon_index-1];
    }

    globalIndex MeshPart::cell_vertex_index(globalIndex const cell_index,
            localIndex const local_corner_index) const {
        assert(local_corner_index < nb_vertices_in_cell(cell_index));
        return cells_connectivity_[cells_ptr_[cell_index] + local_corner_index];
    }

    globalIndex MeshPart::polygon_vertex_index(globalIndex const polygon_index,
            localIndex const local_corner_index) const {
        assert(local_corner_index < nb_vertices_in_polygon(polygon_index));
        return polygons_connectivity_[polygons_ptr_[polygon_index] + local_corner_index];
    }

    std::vector<real64> MeshPart::vertex(globalIndex const vertex_index) const {
        assert(vertex_index < nb_vertices_);
        std::vector<real64> vertex(vertices_.begin()+3*vertex_index,
                vertices_.begin() + 3*vertex_index+3);
        return vertex;
    }

    globalIndex MeshPart::surface(globalIndex const polygon_index) const {
        assert(polygon_index < nb_polygons_);
        return surfaces_indexes_[polygon_index];
    }

    globalIndex MeshPart::region(globalIndex const cell_index) const {
        assert(cell_index < nb_cells_);
        return regions_indexes_[cell_index];
    }

    globalIndex MeshPart::global_vertex_index(globalIndex const vertex_index) const {
        assert(vertex_index < nb_vertices_);
        return original_vertices_indexes_[vertex_index];
    }

    globalIndex MeshPart::global_polygon_index(globalIndex const polygon_index) const {
        assert(polygon_index<nb_polygons_);
        return original_polygons_indexes_[polygon_index];
    }

    globalIndex MeshPart::global_cell_index(globalIndex const cell_index) const {
        assert(cell_index < nb_cells_);
        return original_cells_indexes_[cell_index];
    }

    void MeshPart::set_vertex(globalIndex const vertex_index,
            std::vector< real64 > const& vertex) {
        assert(vertex.size() == 3); 
        assert(vertex_index < nb_vertices_);
        for( localIndex coor = 0 ; coor < 3 ; coor ++ ) {
            vertices_[3*vertex_index+coor] = vertex[coor];
        }
    }

    void MeshPart::set_vertex_original_index( globalIndex const vertex_index,
            globalIndex const original_index) {
        assert(vertex_index < nb_vertices_);
        original_vertices_indexes_[vertex_index] = original_index;
    }

    void MeshPart::set_nb_vertices(globalIndex const nb_vertices) {
        nb_vertices_ = nb_vertices;
        vertices_.resize( nb_vertices*3);
        original_vertices_indexes_.resize( nb_vertices );
    }

    void MeshPart::reserve_nb_cells_and_polygons(globalIndex const nb_elements) {
        cells_ptr_.reserve( nb_elements +1);
        polygons_ptr_.reserve( nb_elements +1);
        cells_connectivity_.reserve( 8 * nb_elements ); // maximum 8 corners (for an hex)
        polygons_connectivity_.reserve( 4 * nb_elements ); // maximum 4 corners (for a quad)
        surfaces_indexes_.reserve( nb_elements);
        regions_indexes_.reserve( nb_elements);
        original_polygons_indexes_.reserve( nb_elements);
        original_cells_indexes_.reserve( nb_elements);
    }
    
    globalIndex MeshPart::add_cell( std::vector<globalIndex> connectivity ) {
        cells_ptr_.push_back( connectivity.size() + cells_ptr_[cells_ptr_.size()-1]);
        nb_cells_++;
        regions_indexes_.resize(nb_cells_);
        original_cells_indexes_.resize(nb_cells_);
        for( localIndex co = 0 ; co < static_cast<localIndex>(connectivity.size()); ++co ) {
            cells_connectivity_.push_back(connectivity[co]);
        }
        return nb_cells_-1;
    }

    globalIndex MeshPart::add_polygon( std::vector<globalIndex> connectivity ) {
        polygons_ptr_.push_back( connectivity.size() + polygons_ptr_[polygons_ptr_.size()-1]);
        nb_polygons_++;
        surfaces_indexes_.resize(nb_polygons_);
        original_polygons_indexes_.resize( nb_polygons_);
        for( localIndex co = 0 ; co < static_cast<localIndex>(connectivity.size()); ++co ) {
            polygons_connectivity_.push_back(connectivity[co]);
        }
        return nb_polygons_-1;
    }

    void MeshPart::set_cell_region( globalIndex const cell_index,
            globalIndex const region_index) {
        assert( cell_index < static_cast< globalIndex> (regions_indexes_.size()) );
        regions_indexes_[cell_index]  = region_index;
    }

    void MeshPart::set_polygon_surface( globalIndex const polygon_index,
            globalIndex const surface_index) {
        assert( polygon_index < static_cast< globalIndex> (surfaces_indexes_.size() ) );
        surfaces_indexes_[polygon_index]  = surface_index;
    }

    void MeshPart::set_cell_original_index(globalIndex const cell_index_in_part_mesh,
                globalIndex const cell_index_in_full_mesh) {
        assert( cell_index_in_part_mesh < static_cast< globalIndex>
                (original_cells_indexes_.size() ) );
        original_cells_indexes_[cell_index_in_part_mesh] = cell_index_in_full_mesh;
    }

    void MeshPart::set_polygon_original_index(globalIndex const polygon_index_in_part_mesh,
                globalIndex const polygon_index_in_full_mesh) {
        assert( polygon_index_in_part_mesh <
                static_cast< globalIndex> (original_polygons_indexes_.size() ) );
        original_polygons_indexes_[polygon_index_in_part_mesh] = polygon_index_in_full_mesh;
    }   

    void MeshPart::finish() {
        regions_indexes_.shrink_to_fit();
        surfaces_indexes_.shrink_to_fit();
        cells_ptr_.shrink_to_fit();
        polygons_ptr_.shrink_to_fit();
        cells_connectivity_.shrink_to_fit();
        polygons_connectivity_.shrink_to_fit();
        original_cells_indexes_.shrink_to_fit();
        original_polygons_indexes_.shrink_to_fit();
    }
}
