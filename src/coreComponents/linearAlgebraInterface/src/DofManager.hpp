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
#include "managers/DomainPartition.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "mesh/MeshLevel.hpp"
#include "dataRepository/ManagedGroup.hpp"

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
  DofManager();

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
   * Assign a mesh
   */
  void setMesh(DomainPartition * const domain, 
               localIndex const meshLevelIndex = 0, 
               localIndex const meshBodyIndex = 0);

  /**
   * Add fields.
   * The user can add a field with a support location, connectivity type, string key, number of scalar components,
   * and a list of element regions over which the field is active.  If the region list is empty, it is assumed
   * the field exists on all regions.
   *
   * The connectivity type is used to infer the sparsity pattern that connects degrees of freedom.
   * If LC denotes a boolean connectivity graph between support locations L and connectors C, the desired sparsity
   * pattern will be computed as LC*CL.  For example, for a TPFA discretization we have dofs located at cell centers,
   * and connected through adjacent faces.  In this example, LC is the cell-to-face connectivity, and LC*CL is the 
   * desired TPFA sparsity pattern.  More generally,
   *
   * Example 1 = ("displacement",NODE,ELEM,3) for a Q1 finite-element interpolation for elasticity
   * Example 2 = ("pressure",ELEM,FACE,1) for a scalar TPFA-type approximation
   * Example 3 = ("pressure",ELEM,NODE,1) for a scalar MPFA-type approximation 
   * Example 4 = ("mass",ELEM,NONE,1) for a diagonal-only sparsity pattern (no connectivitys)
   *
   * When the number of components is greater than one, we always assume they are tightly coupled to one another
   * and form a dense block.  The sparsity pattern LC*CL is then interpreted as the super-node pattern, containing
   * dense sub-blocks.
   */
  void addField(string const & field, 
                Location const location,
                Connectivity const connectivity, 
                localIndex const components = 1,
                string_array const & regions = string_array() );

  /**
   * Add coupling between two fields.
   * The connectivity argument defines how the two fields couple. If the first field has support location A,
   * the second field has support location B, and the connecting object is C, the sparsity pattern will be
   * defined as (AC)(CB).  The final argument indicates if the coupling is symmetric, in the sense that there
   * is a two-way coupling between the fields.  Without this argument, a nonzero block will be added to the 
   * system matrix for block AB, but block BA will remain zero (one-way coupling). 
   *
   * Example 1 = ("node_field","elem_field", ELEM, true) couples all dofs sharing a common element (two-way coupling) 
   * Example 2 = ("node_field_1","node_field_2", NODE, true) couples all dofs sharing a common node (two-way coupling)
   * Example 3 = ("node_field_1","face_field", NODE, false) couples nodal dofs to adjacent faces (one-way coupling)
   */
   
  void addCoupling(string const & rowField, 
                   string const & colField,
                   Connectivity const connectivity,
                   bool const symmetric = false);

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
  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
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
  void getIndices( globalIndex_array & indices,
                   Connectivity const connectivity,
                   localIndex const index,
                   string const & field = "all") const;

private:
  /**
   *  Limit on max number of fields
   */
  localIndex const static MAX_NUM_FIELDS = 10;

  /**
   * Pointer to domain manager
   */
  DomainPartition * m_domain = nullptr;

  /**
   * Pointer to corresponding MeshLevel
   */
  MeshLevel * m_meshLevel = nullptr;

  /**
   * Field description
   */
  struct FieldDescription
  {
    string name;                 //!< field name
    array1d<string> regionNames; //!< active element regions
    array1d<ElementRegion*> regionPtrs; //!< saved pointers to active regions
    Location location;           //!< support location
    localIndex numComponents;    //!< number of vector components
    string key;                  //!< string key for index array
    string docstring;            //!< documentation string
    localIndex numLocalRows;     //!< number of local rows
    globalIndex numGlobalRows;   //!< number of ghost rows
    globalIndex firstLocalRow;   //!< first row on this processor (without field offset)
    globalIndex fieldOffset;     //!< global row offset for multi-field problems

  };

  /**
   * Array of field descriptions
   */
   array1d<FieldDescription> m_fields;

   /*
    * Table of connectivities within and between fields
    */
   array2d<Connectivity> m_connectivity;

   /*
    * Number of MPI ranks
    */
   int mpiSize;

   /*
    * This mpi rank
    */
   int mpiRank;

   /**
    * Check if string key is already being used
    */
   bool keyInUse(string const & key);

   /**
    * Get field index from string key
    */
   localIndex fieldIndex(string const & key);

   /**
    * Create element index array
    */
   void createElemIndexArray(FieldDescription & field);

   /**
    * Create face index array
    */
   void createFaceIndexArray(FieldDescription & field);

   /**
    * Create node index array
    */
   void createNodeIndexArray(FieldDescription & field);
};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_DOFMANAGER_HPP_ */

