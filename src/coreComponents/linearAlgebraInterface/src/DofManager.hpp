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
 * @file DofManager.hpp
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"

namespace geosx
{

/**
 * Defines a simple GEOSX sparsity pattern.  
 * It is intended to be a lightweight way to communicate between the
 * DoFManager and the LAI implementation, without relying on a specific
 * LAI choice.
 *
 * TODO: replace with proper CrsArray 
 */
struct SparsityPattern
{
  array1d<localIndex> rowLengths;  //<! row lengths, size numLocalRows
  array1d<globalIndex> colIndices; //<! packed column indices, size numLocalNonZeros
};


/**
 * The DoFManager is responsible for allocating global dofs, constructing 
 * sparsity patterns, and generally simplifying the interaction between
 * PhysicsSolvers and linear algebra operations.
 */
class DofManager
{
public:

  /**
   * Constructor
   */
  DofManager(MeshLevel const & meshLevel);

  /**
   * Destructor
   */
  ~DofManager() = default;

  /**
   * Enumeration of geometric objects for support location.  Note that this enum
   * is nearly identical to Connectivity, but we keep both for code readability
   * in function calls.
   */
  enum class Location {Elem,Face,Node};    

  /**
   * Enumeration of geometric objects for connectivity type.  Note that this enum
   * is nearly identical to Location, but we keep both for code readability
   * in function calls.
   */
  enum class Connectivity {Elem,Face,Node,None};

  /**
   * Add fields.
   * The user can add a field with a support location, connectivity type, string key, and number of scalar components.
   *
   * The connectivity type is used to infer the sparsity pattern that connects degrees of freedom.
   * If LC denotes a boolean connectivity graph between support locations L and connectors C, the desired sparsity
   * pattern will be computed as LC*CL.  For example, for a TPFA discretization we have dofs located at cell centers,
   * and connected through adjacent faces.  In this example, LC is the cell-to-face connectivity, and LC*CL is the 
   * desired TPFA sparsity pattern.  More generally,
   *
   * Example 1 = add("displacement",NODE,ELEM,3) for a Q1 finite-element interpolation for elasticity
   * Example 2 = add("pressure",ELEM,FACE,1) for a scalar TPFA-type approximation
   * Example 3 = add("pressure",ELEM,NODE,1) for a scalar MPFA-type approximation 
   * Example 4 = add("mass",ELEM,NONE,1) for a diagonal-only sparsity pattern (no connectivitys)
   *
   * When the number of components is greater than one, we always assume they are tightly coupled to one another
   * and form a dense block.  The sparsity pattern LC*CL is then interpreted as the super-node pattern, containing
   * dense sub-blocks.
   */
  void add(string const & field, 
           Location location, 
           Connectivity connectivity, 
           integer const components = 1 );

  /**
   * Close manager.  
   * After all desired fields have been added, this function closes the manager and actually determines
   * global indices for the registered fields.  
   */
  void close();

  /**
   * Return global number of dofs across all processors. If field argument is "all", return monolithic size.
   */
  globalIndex numGlobalDofs( string const & field = "all" ) const;

  /**
   * Return local number of dofs on this processor.  If field argument is "all", return monolithic size
   */
  localIndex numLocalDofs( string const & field = "all" ) const;

  /**  
   * Get a sparsity pattern.  Without additional arguments, this function routines the sparsity pattern for the
   * monolithic matrix.  Sub-patterns can be extracted, however, using row and column field keys.
   */
  void getSparsityPattern( SparsityPattern & pattern, 
                           string const & rowField = "all", 
                           string const & colField = "all") const;

  /**
   * Get global indices for dofs connected by the connector type.  We have two versions, since cells need
   * three indices while faces and nodes only need two.  This keeps the interface the same, but we will only 
   * implement appropriate combinations.
   *
   * Example 1 = getIndices(indices,ELEM,er,esr,ei,"pressure") = get pressure indices connected to this cell 
   * Example 2 = getIndices(indices,FACE,fi,"pressure") = get pressure indices connected to this face
   * Example 3 = getIndices(indices,NODE,ni,"pressure") = get pressure indices connected to this node
   */
  void getIndices( array1d<globalIndex> & indices,
                   Connectivity connectivity,
                   localIndex const region,
                   localIndex const subregion,
                   localIndex const index,
                   string const & field = "all") const;

  /**
   * Get global indices for dofs connected by the connector type.  We have two versions, since cells need
   * three indices while faces and nodes only need two.  This keeps the interface the same, but we will only 
   * implement appropriate combinations.
   *
   * Example 1 = getIndices(indices,ELEM,er,esr,ei,"pressure") = get pressure indices connected to this cell 
   * Example 2 = getIndices(indices,FACE,fi,"pressure") = get pressure indices connected to this face
   * Example 3 = getIndices(indices,NODE,ni,"pressure") = get pressure indices connected to this node
   */
  void getIndices( array1d<globalIndex> & indices,
                   Connectivity connectivity,
                   localIndex const index,
                   string const & field = "all") const;

private:
  /**
   * Reference to corresponding MeshLevel
   */
  MeshLevel const & m_meshLevel;

  /**
   * Field description
   */
  struct FieldDescription
  {
    string name;
    Location location;
    Connectivity connectivity;
    integer components;
  };

  /**
   * Array of field descriptions
   */
  array1d<FieldDescription> m_field;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_ */

