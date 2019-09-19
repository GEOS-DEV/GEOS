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

#include "common/GeosxConfig.hpp"

#include "Blueprint.hpp"

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/NodeManager.hpp"

#ifdef GEOSX_USE_ATK
#include "axom/sidre/core/sidre.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"
#endif

#include <cstring>
#include <unordered_map>
#include <utility>
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



Blueprint::Blueprint( const NodeManager& node_manager,
                      const ElementRegionManager& elem_reg_manager,
                      const std::string& output_path,
                      MPI_Comm GEOSX_UNUSED_ARG( comm ),
                      const std::string& coord_name,
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



void Blueprint::write( int cycle,
                       integer const eventCounter ) const
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
    std::ostringstream msg;
    msg << "does not conform to the blueprint:";
    info.to_json_stream(msg);
    GEOS_LOG_RANK(msg.str());

    GEOS_ERROR(msg.str());
  }

  conduit::Node root_node;
  Node & index = root_node[ "blueprint_index" ];

  conduit::blueprint::mesh::generate_index( mesh_node[ mesh_name ], mesh_name, 1, index[ mesh_name ] );

  info.reset();
  if ( !conduit::blueprint::mesh::index::verify( index[ mesh_name ], info ) )
  {
    std::ostringstream msg;
    msg << "index does not conform to the blueprint:";
    info.to_json_stream(msg);
    GEOS_LOG_RANK(msg.str());

    GEOS_ERROR(msg.str());
  }


  // Build the file-name
  char baseFileName[200] = { 0 };

  sprintf(baseFileName, "_%06d%03d", cycle, eventCounter);

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

  for ( const std::pair< const std::string, const dataRepository::WrapperBase* >& pair : m_node_manager.wrappers() )
  {
    const dataRepository::WrapperBase* view = pair.second;
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


void Blueprint::addCells( Group* GEOSX_UNUSED_ARG( topo ), Group* GEOSX_UNUSED_ARG( fields ) ) const
{
  GEOSX_DEBUG_VAR( m_elem_reg_manager );
#ifdef GEOSX_USE_ATK
  // if ( m_elem_reg_manager.numCellBlocks() != 1 )
  // {
  //   GEOS_ERROR( "Blueprint IO currently only works in problems with one cell block." );
  // }

  // const ElementRegionBase* elem_region = m_elem_reg_manager.GetRegion(0);
  // const CellElementSubRegion* cell_block = elem_region->GetSubRegion<CellElementSubRegion>(0);
  // const array2d<localIndex>& connectivity = cell_block->nodeList().Base();
  // const localIndex n_cells = connectivity.size(0);
  // const localIndex n_nodes_per_cell = connectivity.size(1);
  // const std::string& elem_name = this->numNodesToElemName.at( n_nodes_per_cell );

  // topo->createView( "coordset" )->setString( m_coord_name );
  // topo->createView( "type" )->setString( "unstructured" );
  // Group* elements_group = topo->createGroup( "elements" );
  // elements_group->createView( "shape" )->setString( elem_name );
  
  // View* connec_view = elements_group->createView( "connectivity" );
  // connec_view->setExternalDataPtr( const_cast< localIndex* >( connectivity.data() ) );
  // connec_view->apply( detail::SidreTT< localIndex >::id, n_cells * n_nodes_per_cell );

  // for ( const std::pair< const std::string, const WrapperBase* >& pair : cell_block->wrappers() )
  // {
  //   const WrapperBase* view = pair.second;
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
