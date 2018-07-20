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
 * @file VtupFile.hpp
 */

#ifndef VTUPFILE_HPP_
#define VTUPFILE_HPP_

#include "pugixml.hpp"

/*!
 * @brief this class stands for the I/O of vtup file
 * @details vtu(p) files is an extension fully supported by the VTK/Paraview
 * framework. Details can be found here : www.vtk.org/VTK/img/file-formats.pdf
 * A vtup file is named here as the "parent" file. This file make reference
 * to several "child" files (.vtu). Each one contains a part of the mesh.
 * @todo the export.
 */
namespace geosx{
class VtupFile {
    public:
        VtupFile() {
        }

        /*!
         * @brief load a .vtup file
         * @param[in] filename the name of the XML vtup file to be loaded
         */
        void load( std::string const & filename);

        /*!
         * @brief save a .vtup file
         * @param[in] filename the name of the XML vtup file to be saved
         */
        void save( std::string const & filename);

    private:
        /*!
         * @brief check if the XML file contains the right nodes
         * @return the number of children
         */
        int check_parent_xml_file_consistency() const;
    private:
        /// This is the parent XML document
        pugi::xml_document vtup_doc_;

        /// Name of the Vertices/Elements attribute on which
        /// the original index is stored
        std::string const str_original_index_ { "original_index" };

        /// Name of the Elements attribute on which
        /// the partition index is stored
        std::string const str_partition_ { "partition" };

        /// Name of the Elements attribute on which
        /// the partition index is stored
        std::string const str_region_ { "region" };
};
}
#endif /*VtupFile.hpp*/
