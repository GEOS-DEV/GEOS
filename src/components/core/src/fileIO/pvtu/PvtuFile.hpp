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
class MeshLevel;

class PvtuFile {
    public:
        PvtuFile() {
        }

        virtual ~PvtuFile(){
        }

        /*!
         * @brief load a .pvtu file
         * @param[in] fileName the name of the XML pvtu file to be loaded
         */
        virtual void Load( string const & fileName);

        /*!
         * @brief save a .pvtu file
         * @param[in] fileName the name of the XML pvtu file to be saved
         */
        virtual void Save( string const & fileName);

    protected:
        /*!
         * @brief check if the XML file contains the right nodes
         * @param[in] pvtuDoc the XML document
         * @param[in] fileName name of the file being loaded
         * @param[in] prefix will be "P" for the parent file
         */
        virtual void CheckXmlFileConsistency(pugi::xml_document const & pvtuDoc,
                string const & fileName,
                string const & prefix = "") const;

    private:
        /*!
         * @brief Specific check for parent file
         * @param[in] pvtuDoc the XML document
         * @param[in] fileName name of the file being loaded
         */
         void CheckXmlParentFileConsistency(pugi::xml_document const & pvtuDoc,
                string const & fileName) const;
        /*!
         * @brief retrieve the list of the vtu files
         * @param[in] pvtuDoc the XML document
         * @param[in,out] vtuFiles vector containing the name of the vtu files (one
         * for each partition)
         */
        void VtuFilesList(
                pugi::xml_document const & pvtuDoc,
                std::vector < string > & vtuFiles ) const;
    protected:
        /*
        /// This is the parent XML document
        pugi::xml_document pvtuDoc_;
        */

        /// Name of the Vertices/Elements attribute on which
        /// the original index is stored
        string const m_strOriginalIndex { "original_index" };

        /// Name of the Elements attribute on which
        /// the partition index is stored
        string const m_strPartition { "partition" };

        /// Name of the Elements attribute on which
        /// the partition index is stored
        string const m_strRegion { "region" };

        std::vector< VtuFile > m_vtuFiles;

        std::vector< string > m_vtuFileNames;
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
         * @brief return the surface of a polygon
         * @param[in] polygonIndex polygon index on this mesh part
         * @return the surface on which the polygon belongs
         */
        globalIndex Surface(globalIndex const polygonIndex) const;

        /*!
         * @brief return the region of a cell
         * @param[in] cellIndex cell index
         * @return the cell on which the region belongs
         */
        globalIndex Region(globalIndex const cellIndex) const;

        /*!
         * @brief return the global index of a vertex in the full mesh
         * @param[in] vertexIndex the vertex index on this mesh part
         * @return the global index in the full mesh
         */
        globalIndex GlobalVertexIndex(globalIndex const vertexIndex) const;

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
         * @brief set the original index of a vertex
         * @param[in] vertexIndex the index of the vertex
         * @param[in] originalIndex the index of the vertex in the full mesh
         */
        void SetVertexOriginalIndex(globalIndex const vertexIndex,
                globalIndex const originalIndex);

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
         * @brief Set the surface index of a polygon
         * @param[in] polygonIndex the index of the polygon within the mesh part
         * @param[in] surfaceIndex the surface index (common to the full mesh)
         */
        void SetPolygonSurface( globalIndex const polygonIndex,
                globalIndex const surfaceIndex);

        /*! 
         * @brief Set the region index of a cell
         * @param[in] cellIndex the index of the polygon within the mesh part
         * @param[in] regionIndex the region index (common to the full mesh)
         */
        void SetCellRegion( globalIndex const cellIndex,
                globalIndex const regionIndex);

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

        /// Contains the original indexes of the vertices in the full mesh
        std::vector< globalIndex > m_originalVertexIndexes;
        
        /// Contains the original indexes of the polygons in the full mesh
        std::vector< globalIndex > m_originalPolygonIndexes;
        
        /// Contains the original indexes of the cells in the full mesh
        std::vector< globalIndex > m_originalCellIndexes;

        /// Contains the surface indexes on which the polygons belong (size : numPolygons)
        std::vector< globalIndex > m_surfaceIndexes;

        /// Contains the region indexes on which the cells belong (size : numCells)
        std::vector< globalIndex > m_regionIndexes;
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
         * @param[in] fileName the name of the XML pvtu file to be loaded
         */
        void Load( string const & fileName) final;

        /*!
         * @brief save a .vtu file
         * @param[in] fileName the name of the XML vtu file to be saved
         */
        void Save( string const & fileName) final;

        void TransferMeshPartToGEOSMesh( MeshLevel * const meshLevel );

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
        void LoadMeshPart(pugi::xml_document const & pvtuDoc);

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
    private:
         MeshPart m_meshPart;
};


}
#endif /*PvtuFile.hpp*/
