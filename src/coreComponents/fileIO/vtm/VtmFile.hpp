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

/*!
 * @brief this class stands for the I/O of vtm file
 * @details vtm files is an extension fully supported by the VTK/Paraview
 * framework. Details can be found here : www.vtk.org/VTK/img/file-formats.pdf
 * A vtm file is named here as the "parent" file. This file make reference
 * to several "child" files (.vtu). Each one contains a part of the mesh.
 */
namespace geosx{
class MeshLevel;

/*
class DumbPropertyManager {
    public:
        void SetName(string const & name);
        void SetValue(globalIndex const elementIndex, double const value);
        void SetSize(globalIndex const size);

        string GetName() const;
        double GetValue(globalIndex const elementIndex) const;
    private:
    string m_name;
    std::vector< double > m_values;
};
*/

/*!
 * @brief This class may be temporary and stock a mesh with vertices, polygons (triangles
 * and quads) and cells (pyramids,  prisms, tetraedrons, hexahedrons)
 */
class DumbMesh {
    public:
        DumbMesh() {
        }

        //////////////////////
        /// MESH ACCESSORS ///
        //////////////////////

        /*!
         * @brief Returns the number of vertices
         */
        localIndex NumVertices() const;

        /*!
         * @brief Returns the number of cells
         */
        localIndex NumCells() const;

        /*!
         * @brief Returns the number of tetrahedrons
         */
        localIndex NumTetra() const;

        /*!
         * @brief Returns the number of hexahedrons
         */
        localIndex NumHex() const;

        /*!
         * @brief Returns the number of prisms
         */
        localIndex NumPrism() const;

        /*!
         * @brief Returns the number of pyramids
         */
        localIndex NumPyr() const;

        /*!
         * @brief Returns the number of triangles
         */
        localIndex NumTri() const;

        /*!
         * @brief Returns the number of quads
         */
        localIndex NumQuad() const;

        /*!
         * @brief Returns the number of polygons
         */
        localIndex NumPolygons() const;

        /*!
         * @brief Returns the cell index within this mesh knowing the tetrahedron index
         * @param[in] tetIndex the index of the tetrahedron within this mesh
         */
        localIndex TetIndexToCellIndex(localIndex const localIndex);

        /*!
         * @brief Returns the cell index within this mesh knowing the hexahedron index
         * @param[in] hexIndex the index of the hexahedron within this mesh
         */
        localIndex HexIndexToCellIndex(localIndex const hexIndex);

        /*!
         * @brief Returns the cell index within this mesh knowing the prism index
         * @param[in] prismIndex the index of the prism within this mesh
         */
        localIndex PrismIndexToCellIndex(localIndex const prismIndex);

        /*!
         * @brief Returns the pyramid index within this mesh knowing the pyramid index
         * @param[in] pyrIndex the index of the tetrahedron within this mesh
         */
        localIndex PyrIndexToCellIndex(localIndex const pyrIndex);

        /*!
         * @brief Returns the polygon index within this mesh knowing the triangle index
         * @param[in] triIndex the index of the triangle within this mesh
         */
        localIndex TriIndexToPolygonIndex(localIndex const triIndex);

        /*!
         * @brief Returns the polygon index within this mesh knowing the quad index
         * @param[in] quadIndex the index of the quad within this mesh
         */
        localIndex QuadIndexToPolygonIndex(localIndex const quadIndex);

        /*!
         * @brief return the number of vertices in a polygon
         * @param[in] polygonIndex the polygon index within the mesh part
         * @return the number of vertices within the polygon
         */
        localIndex NumVerticesInPolygon(localIndex const polygonIndex) const;

        /*!
         * @brief return the number of vertices in a cell
         * @param[in] cellIndex the cell index within the mesh part
         * @return the number of vertices within the cell
         */
        localIndex NumVerticesInCell(localIndex const cellIndex) const;

        /*!
         * @brief return the index of a cell corner
         * @param[in] cellIndex the cell index on this mesh part
         * @param[in] localCornerIndex the index of the corner within the cell
         * @return the index of the vertex
         */
        localIndex CellVertexIndex(localIndex const cellIndex,
                localIndex const localCornerIndex) const;

        /*!
         * @brief return the index of a polygon corner
         * @param[in] polygonIndex the polygon index on this mesh part
         * @param[in] localCornerIndex the index of the corner within the polygon
         * @return the index of the veertex
         */
        localIndex PolygonVertexIndex(localIndex const polygonIndex,
                localIndex const localCornerIndex) const;

        /*!
         * @brief returns the vertex coordinate
         * @param[in] vertexIndex vertex index on this mesh part
         * @return a pointer to the first coordinate (X)
         */
        real64 const * Vertex(localIndex const vertexIndex) const;

        /*!
         * @brief return the global index of a cell in the full mesh
         * @param[in] cellIndex the cell index on this mesh part
         * @return the cell index in the full mesh
         */
        globalIndex GlobalCellIndex(localIndex const cellIndex) const;

        /*!
         * @brief return the global index of a polygon in the full mesh
         * @param[in] polygonIndex the polygon index on this mesh part
         * @return the polygon index in the full mesh
         */
        globalIndex GlobalPolygonIndex(localIndex const polygonIndex) const;

        string const & Name() const;

        /////////////////////
        /// MESH BUILDERS ///
        /////////////////////
        
        /*!
         * @brief Set the number of vertices to resize the vector before filling it
         * @details Because we work in 3D, the size of the vector containing the
         * vertices will be 3* numVertices
         * @param[in] numVertices the number of vertices
         */
        void SetNumVertices(localIndex const numVertices); 

        /*!
         * @brief Set vertex coordinates
         * @param[in] vertexIndex the index of the vertex
         * @param[in] vertex a vector containing the three coordinates
         */
        void SetVertex(localIndex const vertexIndex,real64_array const & vertex);

        /*!
         * @brief Set the number of cells and polygons
         */
        void SetNumCellAndPolygons(localIndex numTetra,
                                   localIndex numHex,
                                   localIndex numPrism,
                                   localIndex numPyr,
                                   localIndex numTri,
                                   localIndex numQuad);

        /*!
         * @brief add a cell to the mesh
         * @param[in] connectivity vector containing the indexes of the vertices
         * which compose the cell
         * @return the index of the cell
         */
        void AddCell(localIndex_array const & connectivity );

        /*!
         * @brief add a polygon to the mesh
         * @param[in] connectivity vector containing the indexes of the vertices
         * which compose the polygon
         * @return the index of the polygon
         */
        void AddPolygon(localIndex_array const & connectivity );

        /*!
         * @brief Set the original index, i.e. the index of a cell in the full mesh
         * @param[in] cellIndexInPartMesh the cell index within the part mesh
         * @param[in] cellIndexInFullMesh the cell index within the full mesh
         */
        void SetCellOriginalIndex(localIndex const cellIndexInPartMesh,
                globalIndex const cellIndexInFullMesh);

        /*!
         * @brief Set the original index, i.e. the index of a polygon in the full mesh
         * @param[in] polygonIndexInPartMesh the polygon index within the part mesh
         * @param[in] polygonIndexInFullMesh the polygon index within the full mesh
         */
        void SetPolygonOriginalIndex(localIndex const polygonIndexInPartMesh,
                globalIndex const polygonIndexInFullMesh);

        void SetName( string const & name);

        void Finish();

    private:
        /// Number of vertices
        localIndex m_numVertices{0};

        /// Number of cells
        localIndex m_numCells{0};

        /// Number of Tetra
        localIndex m_numTetra{0};

        /// Number of Hex
        localIndex m_numHex{0};

        /// Number of Prisms
        localIndex m_numPrism{0};

        /// Number of Pyramids
        localIndex m_numPyr{0};

        /// Number of Triangles
        localIndex m_numTri{0};

        /// Number of Quads
        localIndex m_numQuad{0};

        /// Number of polygons
        localIndex m_numPolygons{0};

        /// Contains the 3D coordinates of the vertices
        real64_array m_vertices;
        
        /// Contains the cells connectivity
        localIndex_array m_cellsConnectivity;
        
        /// Contains the polygons connectivity
        localIndex_array m_polygonsConnectivity;

        /// Size : numCells +1. Contains the first indexes to look at in m_cellsConnectivity
        localIndex_array m_cellsPtr{0};
        
        /// Size : numPolygons +1. Contains the first indexes to look at in m_cellsConnectivity
        localIndex_array m_polygonsPtr{0};

        /// Contains the global indexes of the polygons in the full mesh
        globalIndex_array m_globalPolygonIndexes;
        
        /// Contains the global indexes of the cells in the full mesh
        globalIndex_array m_globalCellIndexes;
        
        /// Maps the index of a tetrahedron to its corresponding cell index in this mesih
        localIndex_array m_tetIndexToCellIndex;
        
        /// Maps the index of a hexahedron to its corresponding cell index in this mesh
        localIndex_array m_hexIndexToCellIndex;

        /// Maps the index of a prism to its corresponding cell index in this mesh
        localIndex_array m_prismIndexToCellIndex;

        /// Maps the index of a pyramid to its corresponding cell index in this mesh
        localIndex_array m_pyrIndexToCellIndex;

        /// Maps the index of a triangle to its corresponding cell index in this mesh
        localIndex_array m_triIndexToPolygonIndex;
        
        /// Maps the index of a quad to its corresponding cell index in this mesh
        localIndex_array m_quadIndexToPolygonIndex;

        /// The name of the mesh
        string m_name;

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
        void Load( string const & fileName,
        std::map<string, real64_array > & properties);

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
         * @brief load the mesh part from the XML document
         * @param[in] pvtuDoc the XML document
         */
        void LoadMesh(pugi::xml_document const & pvtuDoc, DumbMesh& mesh);

        /*!
         * @brief load the properties from the XML document
         * @param[in] pvtuDoc the XML document
         */
        void LoadProperties(pugi::xml_document const & pvtuDoc,
                std::map<string, real64_array >& properties);


        /*!
         * @brief split a "big" string contained in a DataArray node of a vtu file
         * @details memory should be reserved for the vector of string
         * @param[in] in the string to be splitter
         * @param[in,out] out the vector of string
         */
        template<typename T, typename Lambda>
            void SplitNodeTextString( string const & in,
                    array1d< T >& out,
                    Lambda && stringConvertor) const;
};


class MeshBlock {
    public:
        MeshBlock( string fileName,
                string blockName);
        MeshBlock() {}
        void Load(bool loadMesh, bool loadProperties);
        DumbMesh const & mesh() const;
        bool IsARegionBlock() const;
        std::map< string, real64_array > const & PropertyMap() const;
    private:
        string m_vtuFileName;
        string m_blockName;
        VtuFile m_vtuFile;
        DumbMesh m_mesh;
        std::map< string, real64_array > m_properties;
};

class RankBlock {
    public:
        void AddMeshBlock( const MeshBlock& block);
        void Load(bool loadMesh, bool loadProperties);
        void TransferRankBlockToGEOSMesh( MeshLevel * const meshLevel ) const;
        localIndex NumMeshBlocks() const;
        MeshBlock const & GetMeshBlock(localIndex const meshBlockIndex) const;
    private:
        array1d< MeshBlock > m_block;
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
         void Load( string const & fileName, bool loadMesh, bool loadProperties);

         void FromVtmToGEOS(MeshLevel * const meshLevel);

         RankBlock const & GetRankBlock(localIndex const rankBlockIndex) const;
         localIndex NumRankBlocks() const;

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
                array1d< RankBlock >& rankBlocks);
    private:
        array1d< RankBlock > m_rankBlocks;
        string m_fileName {""};
};

}
#endif /*PvtuFile.hpp*/
