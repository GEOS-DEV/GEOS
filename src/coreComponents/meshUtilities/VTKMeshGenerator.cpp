/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKMeshGenerator.cpp
 */

#include "VTKMeshGenerator.hpp"

#include "Elements/Element.hpp"
#include "MeshDataWriters/Variable.hpp"
#include "managers/DomainPartition.hpp"

#include <math.h>

#include "mpiCommunications/PartitionBase.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "mesh/MeshBody.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>

namespace geosx
{
using namespace dataRepository;



VTKMeshGenerator::VTKMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{

  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the mesh file" );
}

VTKMeshGenerator::~VTKMeshGenerator()
{}

void VTKMeshGenerator::postProcessInput()
{
  GEOSX_LOG_RANK_0("debut");
  string_array filePathTokenized = stringutilities::tokenize( m_filePath, "."); //TODO maybe code a method in Path to get the file extension?
  string extension = filePathTokenized[filePathTokenized.size() - 1];
  if( extension == "vtk")
  {
    vtkSmartPointer<vtkUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	vtkUgReader->SetFileName( m_filePath.c_str() );
    vtkUgReader->Update();
	m_vtkMesh = vtkUgReader->GetOutput();
  }
  else if( extension == "vtu")
  {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    GEOSX_LOG_RANK_0(m_filePath);
	vtkUgReader->SetFileName( m_filePath.c_str() );
    vtkUgReader->Update();
	m_vtkMesh = vtkUgReader->GetOutput();
  }
  else
  {
    GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
  }
}

Group * VTKMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  GEOSX_LOG_RANK_0("Writing " << m_vtkMesh->GetNumberOfPoints() << " vertices");
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  nodeManager.resize( m_vtkMesh->GetNumberOfPoints() );
  for( vtkIdType v = 0; v < m_vtkMesh->GetNumberOfPoints(); v++)
  {
    double * point = m_vtkMesh->GetPoint( v );
    for( integer i = 0; i < 3; i++)
    {
      X(v,i) = point[i];
    }
  }


}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )
}
