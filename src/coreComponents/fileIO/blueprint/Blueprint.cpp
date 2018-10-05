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

#include "common/GeosxConfig.hpp"

#include "Blueprint.hpp"
#include "common/Logger.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

#ifdef GEOSX_USE_ATK
#include "sidre/sidre.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"
#include "sidre/IOManager.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"
#endif

#include <cstring>
#include <unordered_map>
#include <utility>
#include <mpi.h>
#include <iostream>

using namespace axom::sidre;

namespace geosx
{

const std::unordered_map< localIndex, const std::string > Blueprint::numNodesToElemName =
{
  { 1, "point" },
  { 2, "line" },
  { 3, "triangle" },
  { 4, "tet" },
  { 8, "hex" }
};



Blueprint::Blueprint( const NodeManager& node_manager, const ElementRegionManager& elem_reg_manager,
                      const std::string& output_path, MPI_Comm comm, const std::string& coord_name,
                      const std::string& topo_name):
#ifdef GEOSX_USE_ATK
  m_node_manager( node_manager ),
  m_elem_reg_manager( elem_reg_manager ),
//  m_comm( comm ),
#endif
  m_output_path( output_path ),
  m_coord_name( coord_name ),
  m_topo_name( topo_name )
{}



void Blueprint::write(int cycle,
                      real64 const & eventProgress) const
{
#ifdef GEOSX_USE_ATK
  const string mesh_name = "bp_mesh";

  DataStore ds;
  Group* root = ds.getRoot()->createGroup( mesh_name );
  Group* coords = root->createGroup( "coordsets/" + m_coord_name );
  Group* topo = root->createGroup( "topologies/" + m_topo_name );
  Group* fields = root->createGroup( "fields");

  addNodes( coords, fields );
  addCells( topo, fields );

  conduit::Node mesh_node;
  conduit::Node info;
  ds.getRoot()->createNativeLayout( mesh_node );
  if( !conduit::blueprint::verify( "mesh", mesh_node[ mesh_name ], info ) )
  {
    std::cout << "does not conform to the blueprint:";
    info.print();
    std::cout << std::endl;

    GEOS_ERROR( "Does not conform to the blueprint. See above errors" );
  }

  conduit::Node root_node;
  Node & index = root_node[ "blueprint_index" ];

  conduit::blueprint::mesh::generate_index( mesh_node[ mesh_name ], mesh_name, 1, index[ mesh_name ] );

  info.reset();
  if ( !conduit::blueprint::mesh::index::verify( index[ mesh_name ], info ) )
  {
    std::cout << "index does not conform to the blueprint:";
    info.print();
    std::cout << std::endl;

    GEOS_ERROR( "Does not conform to the blueprint. See above errors" );
  }


  // Build the file-name
  char baseFileName[200] = { 0 };

  integer eventProgressPercent = static_cast<integer>(eventProgress * 100.0);
  sprintf(baseFileName, "_%03d_%06d", eventProgressPercent, cycle);

  const std::string root_output_path = m_output_path + baseFileName + ".root";
  const std::string output_path = m_output_path + baseFileName + ".hdf5";

  root_node[ "protocol/name" ] = "conduit_hdf5";
  root_node[ "protocol/version" ] = "0.1";

  root_node[ "number_of_files" ] = 1;
  root_node[ "number_of_trees" ] = 1;
  root_node[ "file_pattern" ] = output_path;
  root_node[ "tree_pattern" ] = "/";

  conduit::relay::io::save( root_node, root_output_path, "hdf5" );
  conduit::relay::io::save( mesh_node, output_path );
#endif /* GEOSX_USE_ATK */
}


void Blueprint::addNodes( Group* coords, Group* fields ) const
{
#ifdef GEOSX_USE_ATK
  coords->createView( "type" )->setString( "explicit" );

  const r1_array& position = m_node_manager.referencePosition();
  const localIndex n_nodes = position.size();

  View* x_view = coords->createView( "values/x" );
  x_view->setExternalDataPtr( const_cast< R1Tensor* >( position.data() ) );
  x_view->apply( axom::sidre::TypeID::DOUBLE_ID, n_nodes, 0, 3 );

  View* y_view = coords->createView( "values/y" );
  y_view->setExternalDataPtr( const_cast< R1Tensor* >( position.data() ) );
  y_view->apply( axom::sidre::TypeID::DOUBLE_ID, n_nodes, 1, 3 );

  View* z_view = coords->createView( "values/z" );
  z_view->setExternalDataPtr( const_cast< R1Tensor* >( position.data() ) );
  z_view->apply( axom::sidre::TypeID::DOUBLE_ID, n_nodes, 2, 3 );

  for ( const std::pair< const std::string, const dataRepository::ViewWrapperBase* >& pair : m_node_manager.wrappers() )
  {
    const dataRepository::ViewWrapperBase* view = pair.second;
    // if ( view->sizedFromParent() == 1 &&
    //      view->size() > 0 &&
    //      view->shouldRegisterDataPtr() &&
    //      view->getName() != m_node_manager.viewKeys.referencePosition.Key() )
    if ( view->getName() == "Acceleration" || view->getName() == "Mass" )
    {
      Group* cur_field_group = fields->createGroup( "nodes_" + view->getName() );
      cur_field_group->createView( "association" )->setString( "vertex" );
      cur_field_group->createView( "topology" )->setString( m_topo_name );
      cur_field_group->createView( "volume_dependent" )->setString( "false" );
      View* data_view = cur_field_group->createGroup( "values" )->createView( "data" );
      view->registerDataPtr( data_view );
    }
  }
#endif /* GEOSX_USE_ATK */
}


void Blueprint::addCells( Group* topo, Group* fields ) const
{
#ifdef GEOSX_USE_ATK
  if ( m_elem_reg_manager.numCellBlocks() != 1 )
  {
    GEOS_ERROR( "Blueprint IO currently only works in problems with one cell block." );
  }

  const ElementRegion* elem_region = m_elem_reg_manager.GetRegion(0);
  const CellBlockSubRegion* cell_block = elem_region->GetSubRegion(0);
  const array2d<localIndex>& connectivity = cell_block->nodeList().Base();
  const localIndex n_cells = connectivity.size(0);
  const localIndex n_nodes_per_cell = connectivity.size(1);
  const std::string& elem_name = this->numNodesToElemName.at( n_nodes_per_cell );

  topo->createView( "coordset" )->setString( m_coord_name );
  topo->createView( "type" )->setString( "unstructured" );
  Group* elements_group = topo->createGroup( "elements" );
  elements_group->createView( "shape" )->setString( elem_name );
  
  View* connec_view = elements_group->createView( "connectivity" );
  connec_view->setExternalDataPtr( const_cast< localIndex* >( connectivity.data() ) );
  connec_view->apply( detail::SidreTT< localIndex >::id, n_cells * n_nodes_per_cell );

  // for ( const std::pair< const std::string, const ViewWrapperBase* >& pair : cell_block->wrappers() )
  // {
  //   const ViewWrapperBase* view = pair.second;
  //   if ( view->sizedFromParent() == 1 &&
  //        view->size() > 0 &&
  //        view->shouldRegisterDataPtr() &&
  //        view->getName() != CellBlockSubRegion::viewKeyStruct::nodeListString &&
  //        view->getName() != CellBlockSubRegion::viewKeyStruct::edgeListString &&
  //        view->getName() != CellBlockSubRegion::viewKeyStruct::faceListString &&
  //        view->getName() != CellBlockSubRegion::viewKeyStruct::constitutiveGroupingString &&
  //        view->getName() != CellBlockSubRegion::viewKeyStruct::constitutiveMapString &&
  //        view->getName() != CellBlockSubRegion::viewKeyStruct::dNdXString ) 
  //   {
  //     Group* cur_field_group = fields->createGroup( "cell_" + view->getName() );
  //     cur_field_group->createView( "association" )->setString( "element" );
  //     cur_field_group->createView( "topology" )->setString( m_topo_name );
  //     cur_field_group->createView( "volume_dependent" )->setString( "false" );
  //     View* data_view = cur_field_group->createGroup( "values" )->createView( "data" );
  //     view->registerDataPtr( data_view );
  //   }
  // }
#endif /* GEOSX_USE_ATK */
}


} /* end namespace geosx */
