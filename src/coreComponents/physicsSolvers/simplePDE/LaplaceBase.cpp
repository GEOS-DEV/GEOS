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
 * @file LaplaceBase.cpp
 */

#include "LaplaceBase.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
}

using namespace dataRepository;

/* 
 * TODO copy and update guide from LaplaceFEM.cpp
*/
//START_SPHINX_INCLUDE_01
LaplaceBase::LaplaceBase( const string & name,
                          Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_timeIntegrationOption( TimeIntegrationOption::ImplicitTransient )
{
  registerWrapper( laplaceBaseViewKeys.timeIntegrationOption.key(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( laplaceBaseViewKeys.fieldVarName.key(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of field variable" );

}
//END_SPHINX_INCLUDE_01


// Destructor
LaplaceBase::~LaplaceBase()
{
  // TODO Auto-generated destructor stub
}

/* REGISTER THE PDE SOLUTION DATA ON THE MESH
   In the LaplaceFEM solver, we compute the solution of the partial differential equation "numerically".
   This means that we are not solving this PDE everywhere,
   we are computing the solution at specific locations in space.
   To do that, we have to register the Laplace solver so that it works on nodes of a mesh.
   This registration process is done here, in three steps:
   1 - for each mesh body (if the mesh is split), we collect he "nodes" (nodes carry their location information),
   2 - On nodes, we register a new property called m_fieldName and give it a type (here, the type is an array of real64)
   3 - We set some information for this property on the nodes: what is their "PlotLevel"? how can they be described?
     The PlotLevel is a flag that instructs GEOSX to export values of this property for instance so that they can be plotted.
     All properties mounted on nodes carry a certain PlotLevel value. This way, every time GEOSX triggers an
     output event (a request to "print out" data), all properties at or above a certain PlotLevel are automatically exported.
     The description here is simply an additional metadata for the newly mounted property.
 */
//START_SPHINX_INCLUDE_02
void LaplaceBase::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    NodeManager & nodes = meshBody.getMeshLevel( 0 ).getNodeManager();

    nodes.registerWrapper< real64_array >( m_fieldName ).
      setApplyDefaultValue( 0.0 ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setDescription( "Primary field variable" );
  } );
}
//END_SPHINX_INCLUDE_02

} // namespace geosx
