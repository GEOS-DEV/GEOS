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
 * @brief this class stands for the I/O of vtm file
 * @details vtm files is an extension fully supported by the VTK/Paraview
 * framework. Details can be found here : www.vtk.org/VTK/img/file-formats.pdf
 * A vtm file is named here as the "parent" file. This file make reference
 * to several "child" files (.vtu). Each one contains a part of the mesh.
 */
namespace geosx{
class MeshLevel;

/*!
 * @brief This class may be temporary and stock a mesh with vertices, polygons (triangles
 * and quads) and cells (pyramids,  prisms, tetraedra, hexaedra)
 */
class DumbMesh {
    public:
        DumbMesh() {
        }

        //////////////////////
        /// MESH ACCESSORS ///
        //////////////////////

        globalIndex NumVertices() const;

        globalIndex NumCells() const;

        globalIndex NumPolygons() const;

        /*!
         * @brief return the number of vertices in a polygon
         * @param[in] polygonIndex the polygon index within the mesh part
         * @return the number of vertices within the polygon
         */
        localIndex NumVerticesInPolygon(globalIndex const polygonIndex) const;

        /*!
         * @brief return the number of vertices in a cell
         * @param[in] cellIndex the cell index within the mesh part
         * @return the number of vertices within the cell
         */
        localIndex NumVerticesInCell(globalIndex const cellIndex) const;

        /*!
         * @brief return the index of a cell corner
         * @param[in] cellIndex the cell index on this mesh part
         * @param[in] localCornerIndex the index of the corner within the cell
         * @return the index of the vertex
         */
        globalIndex CellVertexIndex(globalIndex const cellIndex,
                localIndex const localCornerIndex) const;

        /*!
         * @brief return the index of a polygon corner
         * @param[in] polygonIndex the polygon index on this mesh part
         * @param[in] localCornerIndex the index of the corner within the polygon
         * @return the index of the veertex
         */
        globalIndex PolygonVertexIndex(globalIndex const polygonIndex,
                localIndex const localCornerIndex) const;

        /*!
         * @brief returns the vertex coordinate
         * @param[in] vertexIndex vertex index on this mesh part
         * @return a pointer to the first coordinate (X)
         */
        std::vector< real64 > Vertex(globalIndex const vertexIndex) const;

        /*!
         * @brief return the global index of a cell in the full mesh
         * @param[in] cellIndex the cell index on this mesh part
         * @return the cell index in the full mesh
         */
        globalIndex GlobalCellIndex(globalIndex const cellIndex) const;

        /*!
         * @brief return the global index of a polygon in the full mesh
         * @param[in] polygonIndex the polygon index on this mesh part
         * @return the polygon index in the full mesh
         */
        globalIndex GlobalPolygonIndex(globalIndex const polygonIndex) const;

        /////////////////////
        /// MESH BUILDERS ///
        /////////////////////
        
        /*!
         * @brief Set the number of vertices to resize the vector before filling it
         * @details Because we work in 3D, the size of the vector containing the
         * vertices will be 3* numVertices
         * @param[in] numVertices the number of vertices
         */
        void SetNumVertices(globalIndex const numVertices); 

        /*!
         * @brief Set vertex coordinates
         * @param[in] vertexIndex the index of the vertex
         * @param[in] vertex a vector containing the three coordinates
         */
        void SetVertex(globalIndex const vertexIndex,std::vector< real64 > const & vertex);

        /*!
         * @brief Reserve the number of cells and polygons
         * @details In pvtu/vtu file format, polygons and cells are stored as "Elements"
         * As a consequence, we do not know a priori the number of polygons and the
         * number of cells. To optimize the std::vector storage, we reserve a number
         * of cells and polygons corresponding to the total number of elements in the
         * mesh part. The connectivity vectors are also reserved.
         */
        void ReserveNumCellAndPolygons(globalIndex const numElements);

        /*!
         * @brief add a cell to the mesh
         * @param[in] connectivity vector containing the indexes of the vertices
         * which compose the cell
         * @return the index of the cell
         */
        globalIndex AddCell( std::vector<globalIndex> connectivity );

        /*!
         * @brief add a polygon to the mesh
         * @param[in] connectivity vector containing the indexes of the vertices
         * which compose the polygon
         * @return the index of the polygon
         */
        globalIndex AddPolygon( std::vector<globalIndex> connectivity );

        /*!
         * @brief Set the original index, i.e. the index of a cell in the full mesh
         * @param[in] cellIndexInPartMesh the cell index within the part mesh
         * @param[in] cellIndexInFullMesh the cell index within the full mesh
         */
        void SetCellOriginalIndex(globalIndex const cellIndexInPartMesh,
                globalIndex const cellIndexInFullMesh);

        /*!
         * @brief Set the original index, i.e. the index of a polygon in the full mesh
         * @param[in] polygonIndexInPartMesh the polygon index within the part mesh
         * @param[in] polygonIndexInFullMesh the polygon index within the full mesh
         */
        void SetPolygonOriginalIndex(globalIndex const polygonIndexInPartMesh,
                globalIndex const polygonIndexInFullMesh);

        /*!
         * @brief Shrink the vector to size
         */
        void Finish();
//    private:
        /// Number of vertices
        globalIndex m_numVertices{0};

        /// Number of cells
        globalIndex m_numCells{0};

        /// Number of polygons
        globalIndex m_numPolygons{0};

        /// Contains the 3D coordinates of the vertices
        array< real64> m_vertices;
        
        /// Contains the cells connectivity
        std::vector< globalIndex > m_cellsConnectivity;
        
        /// Contains the polygons connectivity
        std::vector< globalIndex > m_polygonsConnectivity;

        /// Size : numCells +1. Contains the first indexes to look at in cells_connectivity_
        std::vector< globalIndex > m_cellsPtr{1,0};
        
        /// Size : numPolygons +1. Contains the first indexes to look at in
        /// polygons_connectivity
        std::vector< globalIndex > m_polygonsPtr{1,0};
        
        /// Contains the original indexes of the polygons in the full mesh
        std::vector< globalIndex > m_globalPolygonIndexes;
        
        /// Contains the original indexes of the cells in the full mesh
        std::vector< globalIndex > m_globalCellIndexes;
};

/*!
 * @brief this class stands for the I/O of vtu file
 * @details vtu(p) files is an extension fully supported by the VTK/Paraview
 * framework. Details can be found here : www.vtk.org/VTK/img/file-formats.pdf
 * @todo the export.
 */
class VtuFile {
    public:
        VtuFile() {
        }

        /*!
         * @brief load a .vtu file
         * @param[in] fileName the name of the XML pvtu file to be loaded
         */
        void Load( string const & fileName, DumbMesh & mesh);

        /*!
         * @brief save a .vtu file
         * @param[in] fileName the name of the XML vtu file to be saved
         */
        void Save( string const & fileName);

        void TransferDumbMeshToGEOSMesh( MeshLevel * const meshLevel );

    private:
        /*!
         * @brief check if the XML file contains the right nodes
         void check_xml_file_consistency(pugi::xml_document const & pvtuDoc,
         string const & fileName,
         string const & prefix = "") const final;
         */

        /*!
         * @brief Specific check for child files
         * @param[in] pvtuDoc the XML document
         * @param[in] fileName name of the file being loaded
         */
        void CheckXmlChildFileConsistency(pugi::xml_document const & pvtuDoc,
                string const & fileName) const;

        /*!
         * @brief load the mesh part form the XML document
         * @param[in] pvtuDoc the XML document
         */
        void LoadMesh(pugi::xml_document const & pvtuDoc, DumbMesh& mesh);

        /*!
         * @brief split a "big" string contained in a DataArray node of a vtu file
         * @details memory should be reserved for the vector of string
         * @param[in] in the string to be splitter
         * @param[in,out] out the vector of string
         */
        template<typename T, typename Lambda>
            void SplitNodeTextString( string const & in,
                    std::vector< T >& out,
                    Lambda && stringConvertor) const;
};

class MeshBlock {
    public:
        MeshBlock( string fileName,
                string blockName);
        void Load();
    private:
        string m_vtuFileName;
        string m_blockName;
        VtuFile m_vtuFile;
        DumbMesh m_mesh;
};

class RankBlock {
    public:
        void AddMeshBlock( const MeshBlock& block);
        void Load();
    private:
        std::vector< MeshBlock > m_block;
};
class VtmFile {
    public:
        VtmFile() {
        }

        ~VtmFile(){
        }

        /*!
         * @brief load a .vtm file
         * @param[in] fileName the name of the XML vtm file to be loaded
         */
         void Load( string const & fileName);

         void FromVtmToGEOS();

    protected:
        /*!
         * @brief check if the XML file contains the right nodes
         * @param[in] pvtuDoc the XML document
         * @param[in] fileName name of the file being loaded
         * @param[in] prefix will be "P" for the parent file
         */
         void CheckXmlFileConsistency(pugi::xml_document const & vtmDoc,
                string const & fileName ) const;

    private:
        /*!
         * @brief retrieve the list of the vtu files, organized per ranks and per blocks
         * @param[in] vtmDoc the XML document
         */
        void SetRanksAndBlocks(
                pugi::xml_document const & vtmDoc,
                std::vector< RankBlock >& rankBlocks);
    private:
        std::vector< RankBlock > m_rankBlocks;
};

}
#endif /*PvtuFile.hpp*/
