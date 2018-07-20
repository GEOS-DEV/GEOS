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
namespace geosx{

// PUBLIC METHODS
void PvtuFile::load( std::string const &filename) {

    // XML file parsing using pugixml tool
    pvtu_doc_.load_file(filename.c_str());

    // Retrieve the number of partitions
    int const nb_total_partitions = nb_partitions();
    std::cout << "There are " << nb_total_partitions << " partitions" << std::endl;

    check_parent_xml_file_consistency();
}

void PvtuFile::save( std::string const &fileName) {
    geos_abort(" pvtu file save is not implemented yet");
}


//PRIVATE METHODS
void PvtuFile::check_parent_xml_file_consistency() const {

    // VTKFile is the main node of the pvtufile
    auto const vtk_file =pvtu_doc_.child("VTKFile");
    if( vtk_file.empty() ) {
        geos_abort("Main node VTKFile not found or empty in the parent file");
    }

    auto const pugrid = vtk_file.child("PUnstructuredGrid");
    if( pugrid.empty() ) {
        geos_abort("Node PUnstructuredGrid not found or empty in the parent file");
    }

    std::string const pugrid_child_names[4] = {"PPointData","PCellData","PPoints","Piece"};
    for( auto pugrid_child_name : pugrid_child_names ) {
        auto const pugrid_child = pugrid.child(pugrid_child_name.c_str());
        if(pugrid_child.empty()) {
            std::string const message = "Node " + pugrid_child_name +
                " not found or empty in the parent file";
            geos_abort(message);
        };
    }

    bool points_have_original_index_property = false;
    std::string const mandatory_attributes[4] =
    {"Name","type","format","NumberOfComponents"};
    for( auto ppdata_property : pugrid.child("PPointData").children() ) {
        if( ppdata_property.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    ppdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    std::string const message ="Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of PPointData";
                    geos_abort(message);
                }
            }
            if( ppdata_property.attribute("Name").value() == str_original_index_) {
                points_have_original_index_property = true;
            }
        }
    }
    if (!points_have_original_index_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_original_index_ + " in PPointData");
    }

    bool cells_have_original_index_property = false;
    bool cells_have_partition_property = false;
    bool cells_have_region_property = false;
    for( auto ppdata_property : pugrid.child("PCellData").children() ) {
        if( ppdata_property.name() == static_cast< std::string >("DataArray")) {
            for( auto mandatory_attribute : mandatory_attributes ){
                auto const attribute =
                    ppdata_property.attribute( mandatory_attribute.c_str());
                if( attribute.empty() ) {
                    std::string const message ="Mandatory attribute " + mandatory_attribute +
                        " does not exist in a DataArray of PCellData";
                    geos_abort(message);
                }
            }
            if( ppdata_property.attribute("Name").value() == str_original_index_) {
                cells_have_original_index_property = true;
            }
            if( ppdata_property.attribute("Name").value() == str_partition_) {
                cells_have_partition_property = true;
            }
            if( ppdata_property.attribute("Name").value() == str_region_) {
                cells_have_region_property = true;
            }
        }
    }
    if (!cells_have_original_index_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_original_index_ + " in PCellData");
    }
    if (!cells_have_region_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_region_ + " in PCellData");
    }
    if (!cells_have_partition_property) {
        geos_abort("Can't find any DataArray which contains the property "+
                str_partition_ + " in PCellData");
    }

    bool ppoint_has_a_pdata_array = false;
    bool ppoint_has_a_pdata_array_with_points = false;
    for(auto ppoint_child : pugrid.child("PPoints").children() ) {
        if( ppoint_child.name() == static_cast< std::string >("PDataArray") ) {
            ppoint_has_a_pdata_array = true;
            if( ppoint_child.attribute("Name").as_string() ==
                    static_cast< std::string >("Points")) {
                ppoint_has_a_pdata_array_with_points = true;
                if (ppoint_child.attribute("NumberOfComponents").as_uint() != 3 ) {
                    geos_abort("GEOSX supports only 3D meshes");
                }
                break;
            }
        }
    }
    if( !ppoint_has_a_pdata_array ) {
        geos_abort("No PDataArray found in PPoints.");
    }
    if( !ppoint_has_a_pdata_array_with_points ) {
        geos_abort("No PDataArray named \"Points\" found");
    }

    for( auto pugrid_child : pugrid.children() ) {
        if(pugrid_child.name() == static_cast< std::string >("Piece")) {
            if( pugrid_child.attribute("Source").empty() ) {
                geos_abort("Piece nodes has to have an attribute Source not empty.");
            }
        }
    }
}
    int PvtuFile::nb_partitions() const {
        int nb_partitions = 0;
        for(auto child : pvtu_doc_.child("VTKFile").child("PUnstructuredGrid").children()) {
            if( child.name() == static_cast< std::string > ("Piece") )
            {
                nb_partitions++;
            }
        }
        return nb_partitions;
    }
}
