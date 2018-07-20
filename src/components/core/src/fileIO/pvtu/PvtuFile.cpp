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
    vtup_doc_.load_file(filename.c_str());

    check_parent_xml_file_consistency();
}

void PvtuFile::save( std::string const &fileName) {
    geos_abort(" vtup file save is not implemented yet");
}


//PRIVATE METHODS
void PvtuFile::check_parent_xml_file_consistency() const {

    // VTKFile is the main node of the vtupfile
    auto const vtk_file =vtup_doc_.child("VTKFile");
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
}
}
