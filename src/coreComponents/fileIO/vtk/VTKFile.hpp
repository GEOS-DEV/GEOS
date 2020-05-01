/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKFile.hpp
 */

#ifndef GEOSX_FILEIO_VTK_VTKFILE_HPP_
#define GEOSX_FILEIO_VTK_VTKFILE_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "mesh/InterObjectRelation.hpp"
#include "codingUtilities/StringUtilities.hpp"

#include "dataRepository/xmlWrapper.hpp"


namespace geosx
{

class DomainPartition;

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
    m_plotLevel = dataRepository::toPlotLevel( plotLevel );
  }

  /*!
   * @brief Set the binary mode
   * @param[in] binary the binary mode
   */
  void SetBinaryMode( const bool binary )
  {
    m_binary = binary;
  }

  /*!
   * @brief Output a file for one time step
   * @param[in] cycle the cycle number
   */
  void Write( double const timeStep,
              DomainPartition const & domain );

private:
  /*!
   * @brief Create a XML Node for DataArray
   * @param[in,out] parent the parent XML node
   * @param[in] type a string containing the type of the field
   * @param[in] name the name of the field
   * @param[in] nbComponents dimension of the field
   * @param[in] p is a parallel data array
   * @return the corresponding xml node
   */
  xmlWrapper::xmlNode CreateDataArray( pugi::xml_node & parent,
                                       string const & type,
                                       string const & name,
                                       int const & nbComponents,
                                       string const & format = "ascii",
                                       bool p = false )
  {
    xmlWrapper::xmlNode dataArrayNode;
    if( p )
    {
      dataArrayNode = parent.append_child( "PDataArray" );
    }
    else
    {
      dataArrayNode = parent.append_child( "DataArray" );
    }
    dataArrayNode.append_attribute( "type" ) = type.c_str();
    dataArrayNode.append_attribute( "Name" ) = name.c_str();
    dataArrayNode.append_attribute( "NumberOfComponents" ) = std::to_string( nbComponents ).c_str();
    dataArrayNode.append_attribute( "format" ) = format.c_str();
    return dataArrayNode;
  }

  /*!
   * @brief Create a XML Node for PDataArray
   * @param[in,out] parent the parent XML node
   * @param[in] type a string containing the type of the field
   * @param[in] name the name of the field
   * @param[in] nbComponents dimension of the field
   * @return the corresponding xml node
   */
  xmlWrapper::xmlNode CreatePDataArray( xmlWrapper::xmlNode & parent,
                                        string const & type,
                                        string const & name,
                                        int const & nbComponents,
                                        string const & format = "ascii" )
  {
    return CreateDataArray( parent, type, name, nbComponents, format, true );
  }

private:
  /// Root file ( .pvd )
  xmlWrapper::xmlDocument m_rootFile;

  /// Unstructured file gathering all vtu files for a time step ( .pvtu )
  pugi::xml_document m_pvtuFile;

  /// Plot level
  dataRepository::PlotLevel m_plotLevel;

  /// Base name of the output
  string m_baseName;

  /// Tells wether or not the output is binary
  bool m_binary;
};
}
#endif /* GEOSX_FILEIO_VTK_VTKFILE_HPP_ */
