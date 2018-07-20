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
         */
        virtual void check_xml_file_consistency(pugi::xml_document const & pvtu_doc) const;

    private:
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
         */
        void check_xml_file_consistency(pugi::xml_document const & pvtu_doc) const final;
};
}
#endif /*PvtuFile.hpp*/
