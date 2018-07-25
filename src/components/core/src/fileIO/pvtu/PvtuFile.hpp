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
 * @file PvtuFile.hpp
 */

#ifndef VTUPFILE_HPP_
#define VTUPFILE_HPP_

#include "pugixml.hpp"
#include "common/DataTypes.hpp"
#include <vector>

/*!
 * @brief this class stands for the I/O of pvtu file
 * @details vtu(p) files is an extension fully supported by the VTK/Paraview
 * framework. Details can be found here : www.vtk.org/VTK/img/file-formats.pdf
 * A pvtu file is named here as the "parent" file. This file make reference
 * to several "child" files (.vtu). Each one contains a part of the mesh.
 * @todo the export.
 */
namespace geosx{
class VtuFile;

class PvtuFile {
    public:
        PvtuFile() {
        }

        virtual ~PvtuFile(){
        }

        /*!
         * @brief load a .pvtu file
         * @param[in] filename the name of the XML pvtu file to be loaded
         */
        virtual void load( std::string const & filename);

        /*!
         * @brief save a .pvtu file
         * @param[in] filename the name of the XML pvtu file to be saved
         */
        virtual void save( std::string const & filename);

    protected:
        /*!
         * @brief check if the XML file contains the right nodes
         * @param[in] pvtu_doc the XML document
         * @param[in] filename name of the file being loaded
         * @param[in] prefix will be "P" for the parent file
         */
        virtual void check_xml_file_consistency(pugi::xml_document const & pvtu_doc,
                std::string const & filename,
                std::string const & prefix = "") const;

    private:
        /*!
         * @brief Specific check for parent file
         * @param[in] pvtu_doc the XML document
         * @param[in] filename name of the file being loaded
         */
         void check_xml_parent_file_consistency(pugi::xml_document const & pvtu_doc,
                std::string const & filename) const;
        /*!
         * @brief retrieve the list of the vtu files
         * @param[in] pvtu_doc the XML document
         * @param[in,out] vtu_files vector containing the name of the vtu files (one
         * for each partition)
         */
        void vtu_files_list(
                pugi::xml_document const & pvtu_doc,
                std::vector < std::string > & vtu_files ) const;
    protected:
        /*
        /// This is the parent XML document
        pugi::xml_document pvtu_doc_;
        */

        /// Name of the Vertices/Elements attribute on which
        /// the original index is stored
        std::string const str_original_index_ { "original_index" };

        /// Name of the Elements attribute on which
        /// the partition index is stored
        std::string const str_partition_ { "partition" };

        /// Name of the Elements attribute on which
        /// the partition index is stored
        std::string const str_region_ { "region" };

        std::vector< VtuFile > vtu_files_;

        std::vector< std::string > vtu_file_names_;
};

/*!
 * @brief This class may be temporary and stock a mesh with vertices, polygons (triangles
 * and quads) and cells (pyramids,  prisms, tetraedra, hexaedra)
 */
class MeshPart {
    public:
        MeshPart() {
        }

        //////////////////////
        /// MESH ACCESSORS ///
        //////////////////////

        globalIndex nb_vertices() const;

        globalIndex nb_cells() const;

        globalIndex nb_polygons() const;

        /*!
         * @brief return the number of vertices in a polygon
         * @param[in] polygon_index the polygon index within the mesh part
         * @return the number of vertices within the polygon
         */
        localIndex nb_vertices_in_polygon(globalIndex const polygon_index);

        /*!
         * @brief return the number of vertices in a cell
         * @param[in] cell_index the cell index within the mesh part
         * @return the number of vertices within the cell
         */
        localIndex nb_vertices_in_cell(globalIndex const cell_index);

        /*!
         * @brief return the index of a cell corner
         * @param[in] cell_index the cell index on this mesh part
         * @param[in] local_corner_index the index of the corner within the cell
         * @return the index of the vertex
         */
        globalIndex cell_vertex_index(globalIndex const cell_index,
                localIndex const local_corner_index) const;

        /*!
         * @brief return the index of a polygon corner
         * @param[in] polygon_index the polygon index on this mesh part
         * @param[in] local_corner_index the index of the corner within the polygon
         * @return the index of the veertex
         */
        globalIndex polygon_vertex_index(globalIndex const polygon_index,
                localIndex const local_corner_index) const;

        /*!
         * @brief returns the vertex coordinate
         * @param[in] vertex_index vertex index on this mesh part
         * @return a pointer to the first coordinate (X)
         */
        std::vector< real64 > vertex(globalIndex const vertex_index) const;

        /*!
         * @brief return the surface of a polygon
         * @param[in] polygon_index polygon index on this mesh part
         * @return the surface on which the polygon belongs
         */
        globalIndex surface(globalIndex const polygon_index) const;

        /*!
         * @brief return the region of a cell
         * @param[in] cell_index cell index
         * @return the cell on which the region belongs
         */
        globalIndex region(globalIndex const cell_index) const;

        /*!
         * @brief return the global index of a vertex in the full mesh
         * @param[in] vertex_index the vertex index on this mesh part
         * @return the global index in the full mesh
         */
        globalIndex global_vertex_index(globalIndex const vertex_index) const;

        /*!
         * @brief return the global index of a cell in the full mesh
         * @param[in] cell_index the cell index on this mesh part
         * @return the cell index in the full mesh
         */
        globalIndex global_cell_index(globalIndex const cell_index) const;

        /*!
         * @brief return the global index of a polygon in the full mesh
         * @param[in] polygon_index the polygon index on this mesh part
         * @return the polygon index in the full mesh
         */
        globalIndex global_polygon_index(globalIndex const polygon_index) const;

        /////////////////////
        /// MESH BUILDERS ///
        /////////////////////
        
        /*!
         * @brief Set the number of vertices to resize the vector before filling it
         * @details Because we work in 3D, the size of the vector containing the
         * vertices will be 3* nb_vertices
         * @param[in] nb_vertices the number of vertices
         */
        void set_nb_vertices(globalIndex const nb_vertices); 

        /*!
         * @brief Set vertex coordinates
         * @param[in] vertex_index the index of the vertex
         * @param[in] vertex a vector containing the three coordinates
         */
        void set_vertex(globalIndex const vertex_index,std::vector< real64 > const & vertex);

        /*!
         * @brief set the original index of a vertex
         * @param[in] vertex_index the index of the vertex
         * @param[in] original_index the index of the vertex in the full mesh
         */
        void set_vertex_original_index(globalIndex const vertex_index,
                globalIndex const original_index);

        /*!
         * @brief Reserve the number of cells and polygons
         * @details In pvtu/vtu file format, polygons and cells are stored as "Elements"
         * As a consequence, we do not know a priori the number of polygons and the
         * number of cells. To optimize the std::vector storage, we reserve a number
         * of cells and polygons corresponding to the total number of elements in the
         * mesh part. The connectivity vectors are also reserved.
         */
        void reserve_nb_cells_and_polygons(globalIndex const nb_elements);

        /*!
         * @brief Set the number of cells to resize the vector before filling it
         * @param[in] nb_cells the number of cells
         */
        //void set_nb_cells(globalIndex const nb_cells );

        /*!
         * @brief Set the number of polygons to resize the vector before filling it
         * @param[in] nb_polygons the number of polygons
         */
        //void set_nb_polygons(globalIndex const nb_polygons );

        /*!
         * @brief set the offset in the cells_ptr_vector
         * @details basically, the offset is the number of corners of the
         * cell_index cell index within the part mesh.
         */
        /*
        void set_cell_ptr_offset(globalIndex const cell_index, globalIndex const offset) ;
        */

        /*!
         * @brief set the offset in the polygons_ptr_vector
         * @details basically, the offset is the number of corners of the
         * polygon_index polygon index within the part mesh
         */
        /*
        void set_polygon_ptr_offset(globalIndex const polygon_index, globalIndex const offset) ;
        */

        /*!
         * @brief Set the connectivity of a cell
         * @param[in] cell_index the index of the cell
         * @param[in] local_corner_index the index of the vertex within the cell
         * @param[in] vertex_index the index of the vertex within the mesh part
         */
        /*
        void set_cell_corner( globalIndex const cell_index,
                localIndex const local_corner_index,
                globalIndex const vertex_index);
                */

        /*!
         * @brief Set the connectivity of a polygon
         * @param[in] polygon_index the index of the polygon on this mesh part
         * @param[in] local_corner_index the index of the vertex within the polygon
         * @param[in] vertex_index the index of the vertex within the mesh part
         */
        /*
        void set_polygon_corner(globalIndex const polygon_index,
                localIndex const local_corner_index,
                globalIndex const vertex_index);
                */

        /*!
         * @brief add a cell to the mesh
         * @param[in] connectivity vector containing the indexes of the vertices
         * which compose the cell
         * @return the index of the cell
         */
        globalIndex add_cell( std::vector<globalIndex> connectivity );

        /*!
         * @brief add a polygon to the mesh
         * @param[in] connectivity vector containing the indexes of the vertices
         * which compose the polygon
         * @return the index of the polygon
         */
        globalIndex add_polygon( std::vector<globalIndex> connectivity );

        /*! 
         * @brief Set the surface index of a polygon
         * @param[in] polygon_index the index of the polygon within the mesh part
         * @param[in] surface_index the surface index (common to the full mesh)
         */
        void set_polygon_surface( globalIndex const polygon_index,
                globalIndex const surface_index);

        /*! 
         * @brief Set the region index of a cell
         * @param[in] cell_index the index of the polygon within the mesh part
         * @param[in] region_index the region index (common to the full mesh)
         */
        void set_cell_region( globalIndex const cell_index,
                globalIndex const region_index);

        /*!
         * @brief Set the original index, i.e. the index of a cell in the full mesh
         * @param[in] cell_index_in_part_mesh the cell index within the part mesh
         * @param[in] cell_index_in_full_mesh the cell index within the full mesh
         */
        void set_cell_original_index(globalIndex const cell_index_in_part_mesh,
                globalIndex const cell_index_in_full_mesh);

        /*!
         * @brief Set the original index, i.e. the index of a polygon in the full mesh
         * @param[in] polygon_index_in_part_mesh the polygon index within the part mesh
         * @param[in] polygon_index_in_full_mesh the polygon index within the full mesh
         */
        void set_polygon_original_index(globalIndex const polygon_index_in_part_mesh,
                globalIndex const polygon_index_in_full_mesh);

        /*!
         * @brief Shrink the vector to size
         */
        void finish();
    private:
        /// Number of vertices
        globalIndex nb_vertices_{0};

        /// Number of cells
        globalIndex  nb_cells_{0};

        /// Number of polygons
        globalIndex nb_polygons_{0};

        /// Contains the 3D coordinates of the vertices
        std::vector< real64> vertices_;
        
        /// Contains the cells connectivity
        std::vector< globalIndex > cells_connectivity_;
        
        /// Contains the polygons connectivity
        std::vector< globalIndex > polygons_connectivity_;

        /// Size : nb_cells +1. Contains the first indexes to look at in cells_connectivity_
        std::vector< globalIndex > cells_ptr_{1,0};
        
        /// Size : nb_polygons +1. Contains the first indexes to look at in
        /// polygons_connectivity
        std::vector< globalIndex > polygons_ptr_{1,0};

        /// Contains the original indexes of the vertices in the full mesh
        std::vector< globalIndex > original_vertices_indexes_;
        
        /// Contains the original indexes of the polygons in the full mesh
        std::vector< globalIndex > original_polygons_indexes_;
        
        /// Contains the original indexes of the cells in the full mesh
        std::vector< globalIndex > original_cells_indexes_;

        /// Contains the surface indexes on which the polygons belong (size : nb_polygons)
        std::vector< globalIndex > surfaces_indexes_;

        /// Contains the region indexes on which the cells belong (size : nb_cells)
        std::vector< globalIndex > regions_indexes_;
};
/*!
 * @brief this class stands for the I/O of vtu file
 * @details vtu(p) files is an extension fully supported by the VTK/Paraview
 * framework. Details can be found here : www.vtk.org/VTK/img/file-formats.pdf
 * @todo the export.
 */
class VtuFile : public PvtuFile{
    public:
        VtuFile() {
        }

        /*!
         * @brief load a .vtu file
         * @param[in] filename the name of the XML pvtu file to be loaded
         */
        void load( std::string const & filename) final;

        /*!
         * @brief save a .vtu file
         * @param[in] filename the name of the XML vtu file to be saved
         */
        void save( std::string const & filename) final;

    private:
        /*!
         * @brief check if the XML file contains the right nodes
         void check_xml_file_consistency(pugi::xml_document const & pvtu_doc,
         std::string const & filename,
         std::string const & prefix = "") const final;
         */

        /*!
         * @brief Specific check for child files
         * @param[in] pvtu_doc the XML document
         * @param[in] filename name of the file being loaded
         */
        void check_xml_child_file_consistency(pugi::xml_document const & pvtu_doc,
                std::string const & filename) const;

        /*!
         * @brief load the mesh part form the XML document
         * @param[in] pvtu_doc the XML document
         */
        void load_mesh_part(pugi::xml_document const & pvtu_doc);

        /*!
         * @brief split a "big" string contained in a DataArray node of a vtu file
         * @details memory should be reserved for the vector of string
         * @param[in] in the string to be splitter
         * @param[in,out] out the vector of string
         */
        template<typename T, typename Lambda>
            void split_node_text_string( std::string const & in,
                    std::vector< T >& out,
                    Lambda && string_convertor) const;
    private:
         MeshPart mesh_part_;
};


}
#endif /*PvtuFile.hpp*/
