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

#include "PvtuFile.hpp"
#include "common/Logger.hpp"
#include <iostream>
#include <string.h>

#if USE_MPI
#include <mpi.h>
#endif

namespace geosx{

// PUBLIC METHODS
void PvtuFile::load( std::string const &filename) {
    int size = 0;
    int rank = 0;
#if USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    int nb_partition_per_core = 0;
    //int nb_total_partitions = 0;
    int remainder = 0;
    std::vector < std::string > children_files;



    // XML file parsing using pugixml tool
    pugi::xml_document pvtu_doc;
    pvtu_doc.load_file(filename.c_str());

    if( rank == 0) {
        check_xml_parent_file_consistency(pvtu_doc,filename);
    }
    vtu_files_list(pvtu_doc,children_files);
    // Retrieve the number of partitions
    int const nb_total_partitions = children_files.size();

    // Next part of this method is dedicated to the optimization of file loading
    //
    // IF nb_partitions > nb_mpi_process : Each process will load 
    // nb_partitions / nb_mpi_process. Processes with a rank < nb_partitions % nb_mpi_process
    // will load an additional partition
    //
    // IF nb_partitions < nb_mpi_process : the first nb_partitions process will
    // load ONE partition.
    if(nb_total_partitions > size) {
        if ( rank == 0 ) {
            std::cout << "WARNING : the number of partitions ("
                << nb_total_partitions <<") which will be loaded " 
                << "is greater that the number of processes on which GEOSX is running ("
                << size
                << "). Some processes will hold more than one partition !" << std::endl;
        }
        nb_partition_per_core = nb_total_partitions / size;
        remainder = nb_total_partitions % size;
    } else {
        nb_partition_per_core = 1;
        remainder = 0;
    }
    if( rank < remainder ) {
        vtu_file_names_.resize(nb_partition_per_core+1);
        vtu_files_.resize(nb_partition_per_core+1);
        for(int p_index = 0; p_index < nb_partition_per_core +1; ++p_index) {
            vtu_file_names_[p_index] =
                children_files[rank*(nb_partition_per_core+1) + p_index];
        }
    } else if( rank < nb_total_partitions) {
        vtu_file_names_.resize(nb_partition_per_core);
        vtu_files_.resize(nb_partition_per_core);
        for(int p_index = 0; p_index < nb_partition_per_core; ++p_index) {
            vtu_file_names_[p_index] =
                children_files[rank*(nb_partition_per_core) + p_index+remainder];
        }
    }
    if( rank < nb_total_partitions ) {
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
}
