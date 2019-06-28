/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file VTKFile.hpp
 */

#ifndef VTKFILE_HPP_
#define VTKFILE_HPP_

#include "common/DataTypes.hpp"
#include "pugixml.hpp"
#include "dataRepository/RestartFlags.hpp" 

#ifdef GEOSX_USE_MPI
#include <mpi.h>
#endif

namespace geosx
{

class VTKFile
{
  public:
  VTKFile() = delete;

  /*!
   * @brief Initialize the VTK file writer
   * @details This constructor will construct the root file for
   * the output (.pvd)
   * @param[in] name the name of the pvd file
   */
  VTKFile( string const & name );

  /*!
   * @brief Set the plot level
   * @param[in] plotLevel the plot level. All fields flagged with
   * a plot level inferior or equal to this prescribed value will
   * be output
   */
  void SetPlotLevel( const int plotLevel )
  {
    m_plotLevel = dataRepository::IntToPlotLevel(plotLevel);
  }

  /*!
   * @brief Output a file for one time step
   * @param[in] cycle the cycle number
   */
  void Write( double const timeStep );

  private:
    /*!
     * @brief Create a XML Node for PDataArray
     * @param[in,out] parent the parent XML node
     * @param[in] type a string containing the type of the field
     * @param[in] name the name of the field
     * @param[in] nbComponents dimension of the field
     */
  void CreatePDataArray( pugi::xml_node & parent,
                         string const & type,
                         string const & name,
                         int const & nbComponents )
  {
    auto pDataArrayNode = parent.append_child("PDataArray");
    pDataArrayNode.append_attribute("type") = type.c_str();
    pDataArrayNode.append_attribute("Name") = name.c_str();
    pDataArrayNode.append_attribute("NumberOfComponents") = std::to_string(nbComponents).c_str();
  }

  private:
    /// Root file ( .pvd )
    pugi::xml_document m_rootFile;

    /// Unstructured file gathering all vtu files for a time step ( .pvtu )
    pugi::xml_document m_pvtuFile;

    /// Plot level
    dataRepository::PlotLevel m_plotLevel;

};
}
#endif /* VTKFILE_H_ */
