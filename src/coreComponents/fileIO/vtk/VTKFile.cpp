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
 * @file VTKFile.cpp
 */

#include "VTKFile.hpp"
#include <sys/stat.h>

namespace geosx
{
  VTKFile::VTKFile( string const & name )
  {
    int mpiRank;
    MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
    if( mpiRank == 0 )
    {
      // Declaration of XML version
      auto declarationNode = m_rootFile.append_child(pugi::node_declaration);
      declarationNode.append_attribute("version") = "1.0";

      // Declaration of the node VTKFile
      auto vtkFileNode = m_rootFile.append_child("VTKFile");
      vtkFileNode.append_attribute("type") = "Collection";
      vtkFileNode.append_attribute("version") = "0.1";
      vtkFileNode.append_attribute("byteOrder") = "LittleEndian";
      vtkFileNode.append_attribute("compressor") = "vtkZLibDataCompressor";

      // Declaration of the node Collection
      vtkFileNode.append_child("Collection");
    }

    //temp
    m_rootFile.save_file("test.pvd");

  }

  void VTKFile::Write( double const timeStep )
  {
    int mpiRank;
    MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
    if( mpiRank == 0 )
    {
      /// Add the new entry to the pvd root file
      auto collectionNode = m_rootFile.child("VTKFile").child("Collection");
      auto dataSetNode = collectionNode.append_child("DataSet");
      dataSetNode.append_attribute("timestep") = std::to_string( timeStep ).c_str();
      dataSetNode.append_attribute("group") = "";
      dataSetNode.append_attribute("part") = "0";
      mode_t mode = 0733;

      /// Create the pvtu file for this time step
      {
        // Create a directory for this time step
        mkdir( std::to_string( timeStep ).c_str(), mode );
        pugi::xml_document pvtuFile;

        // Declaration of XML version
        auto declarationNode = pvtuFile.append_child(pugi::node_declaration);
        declarationNode.append_attribute("version") = "1.0";

        // Declaration of the node VTKFile
        auto vtkFileNode = pvtuFile.append_child("VTKFile");
        vtkFileNode.append_attribute("type") = "PUnstructuredGrid";
        vtkFileNode.append_attribute("version") = "0.1";
        vtkFileNode.append_attribute("byteOrder") = "LittleEndian";

        // Declaration of the node PUnstructuredGrid
        auto pUnstructureGridNode = vtkFileNode.append_child("PUnstructuredGrid");
        pUnstructureGridNode.append_attribute("GhostLevel") = "0";

        // Declaration the node PPoints
        auto pPointsNode = pUnstructureGridNode.append_child("PPoints");
        // .... and the data array containg the positions
        CreatePDataArray( pPointsNode, "Float64", "Position", 3 );

        // Declaration of the node PCells
        auto pCellsNode = pUnstructureGridNode.append_child("PCells");
        // .... and its data array defining the connectivities, types, and offsets
        CreatePDataArray( pCellsNode, "Int32", "Position", 3 );
        CreatePDataArray( pCellsNode, "Int32", "Position", 3 );
        CreatePDataArray( pCellsNode, "Uint8", "Position", 3 );

      }

    }


  }


}
