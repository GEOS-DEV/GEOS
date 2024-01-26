/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKHierarchicalDataRepository.cpp
 */

#include "mesh/generators/VTKHierarchicalDataRepository.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include <vtkXMLPartitionedDataSetCollectionReader.h>

namespace geos
{
using namespace dataRepository;

VTKHierarchicalDataRepository::VTKHierarchicalDataRepository( string const & name,
                                                              Group * const parent )
  : ExternalDataRepositoryBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setApplyDefaultValue( "attribute" ).
    setDescription( "Path to the mesh file" );
}

void VTKHierarchicalDataRepository::open()
{
  string const extension = m_filePath.extension();
  GEOS_ERROR_IF( extension != "vtpc", "Unsupported vtk extension. File must be a vtpc file" );

  vtkNew< vtkXMLPartitionedDataSetCollectionReader > reader;
  reader->SetFileName( m_filePath.c_str());
  reader->Update();

  m_collection = vtkSmartPointer< vtkPartitionedDataSetCollection >( vtkPartitionedDataSetCollection::SafeDownCast( reader->GetOutput()));
  m_dataAssembly = vtkSmartPointer< vtkDataAssembly >( m_collection->GetDataAssembly() );
}

vtkSmartPointer< vtkPartitionedDataSet >
VTKHierarchicalDataRepository::search( string const & path )
{
  int node = m_dataAssembly->GetFirstNodeByPath( path.c_str());
  GEOS_ERROR_IF( node == -1, "Node doesn't exist" );
  GEOS_ERROR_IF( m_dataAssembly->GetNumberOfChildren( node ) > 0, "only leaf nodes can be queried." );

  return m_collection->GetPartitionedDataSet( node );

  // for(auto & sub_node : m_dataAssembly->GetChildNodes(node, false))
  // {
  //   std::vector<unsigned int> datasets = m_dataAssembly->GetDataSetIndices(sub_node, false);
  //   dataset = pdsc.GetPartitionedDataSet(node)
  // }
  // {
  //   d
  //   for d in datasets
  //   {
  //     dataset = pdsc.GetPartitionedDataSet(d)
  //     grid = pv.wrap(dataset.GetPartition(0))
  //     # grid.scale([1.0, 1.0, args.Zamplification], inplace=True)
  //     region_engine.add_mesh(grid)
  //   }
  // }
  //TODO
}

REGISTER_CATALOG_ENTRY( ExternalDataRepositoryBase, VTKHierarchicalDataRepository, string const &, Group * const )


}
