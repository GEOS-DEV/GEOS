/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DofManager.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_DOFMANAGER_HPP_
#define GEOSX_LINEARALGEBRA_DOFMANAGER_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/utilities/ComponentMask.hpp"

#include <numeric>

namespace geosx
{

class DomainPartition;
class MeshLevel;
class ObjectManagerBase;
class FluxApproximationBase;

/**
 * @class DofManager
 * @brief The DoFManager is responsible for allocating global dofs, constructing
 * sparsity patterns, and generally simplifying the interaction between
 * PhysicsSolvers and linear algebra operations.
 */
class DofManager
{
public:

  /// Maximum number of components in a field
  static int constexpr MAX_COMP = 32;

  /// Type of component mask used by DofManager
  using CompMask = ComponentMask< MAX_COMP >;

  /**
   * @brief Describes a selection of components from a DoF field.
   *
   * A half-open range [@p loComp, @p hiComp) is selected.
   */
  struct SubComponent
  {
    string fieldName;  ///< Name of the DOF field in DofManager
    CompMask mask;     ///< Mask that defines component selection
  };

  /**
   * @brief Enumeration of geometric objects for support location. Note that this
   * enum is nearly identical to Connectivity, but we keep both for code readability
   * in function calls.
   */
  enum class Location
  {
    Elem, //!< location is element (like pressure in finite volumes)
    Face, //!< location is face (like flux in mixed finite elements)
    Edge, //!< location is edge (like flux between fracture elements)
    Node  //!< location is node (like displacements in finite elements)
  };

  /**
   * @brief Enumeration of geometric objects for connectivity type. Note that this
   * enum is nearly identical to Location, but we keep both for code readability
   * in function calls.
   */
  enum class Connector
  {
    Elem, //!< connectivity is element (like in finite elements)
    Face, //!< connectivity is face (like in finite volumes TPFA)
    Edge, //!< connectivity is edge (like fracture element connectors)
    Node, //!< connectivity is node (like in finite volumes MPFA)
    None,  //!< there is no connectivity (self connected field, like a lumped mass matrix)
    Stencil //!< connectivity is through a (set of) user-provided stencil(s)
  };

  /**
   * @brief Constructor.
   *
   * @param [in] name a unique name for this DoF manager
   */
  explicit DofManager( string name );

  /**
   * @brief Deleted copy constructor.
   */
  DofManager( DofManager const & ) = delete;

  /**
   * @brief Move constructor.
   */
  DofManager( DofManager && ) = default;

  /**
   * @brief Deleted copy assignment.
   * @return
   */
  DofManager & operator=( DofManager const & ) = delete;

  /**
   * @brief Defaulted move assignment.
   * @return
   */
  DofManager & operator=( DofManager && ) = default;

  /**
   * @brief Destructor.
   */
  ~DofManager() = default;

  /**
   * @brief Remove all fields and couplings and re-enable addition of new fields.
   */
  void clear();

  /**
   * @brief Assign a mesh.
   * @param mesh the target mesh (non-const access required for allocating index arrays)
   * @note Calling this function unconditionally removes all previously defined fields and couplings.
   *       The user should only call this when they actually want to replace the mesh with a new one,
   *       or they think the mesh topology might have changed and so the DOFs need to be re-numbered.
   *       They must then re-add all fields and couplings, and call reorderByRank() again.
   */
  void setMesh( MeshLevel & mesh );

  /**
   * @brief Add a new field and enumerate its degrees-of-freedom.
   *
   * @param [in] fieldName the name of the field
   * @param [in] location type of mesh objects the field is defined on
   * @param [in] components number of components
   * @param [in] regions list of region names the field is defined on (if empty, selects all regions)
   */
  void addField( string const & fieldName,
                 Location location,
                 integer components,
                 std::vector< string > const & regions = {} );

  /**
   * @copydoc addField(string const &, Location, integer, std::vector< string > const &)
   *
   * Overload for arrayView1d<string> input used by physics solvers.
   */
  void addField( string const & fieldName,
                 Location location,
                 integer components,
                 arrayView1d< string const > const & regions );

  /**
   * @brief Add coupling between two fields.
   *
   * The connectivity argument defines how the two fields couple. If the first field has support location A,
   * the second field has support location B, and the connecting object is C, the sparsity pattern will be
   * defined as (AC)(CB). The final argument indicates if the coupling is symmetric, in the sense that there
   * is a two-way coupling between the fields. Without this argument, a nonzero block will be added to the
   * system matrix for block AB, but block BA will remain zero (one-way coupling).
   *
   * - Example 1 = ("node_field","elem_field", ELEM, {}, true) couples all dofs sharing a common element (two-way coupling)
   * - Example 2 = ("node_field","node_field", NODE, {}, true) couples all dofs sharing a common node (two-way coupling)
   * - Example 3 = ("node_field","face_field", NODE, {}, false) couples nodal dofs to adjacent faces (one-way coupling)
   *
   * When the number of components is greater than one, we always assume they are tightly coupled to one another
   * and form a dense block. The sparsity pattern LC*CL is then interpreted as the super-node pattern, containing
   * dense sub-blocks.
   *
   * @param [in] rowFieldName string the name of the row field.
   * @param [in] colFieldName string the name of the col field.
   * @param [in] connectivity Connectivity through what they are connected.
   * @param [in] regions names of regions where this coupling is defined.
   * @param [in] symmetric bool is it symmetric, i.e., both row-col and col-row?
   */
  void addCoupling( string const & rowFieldName,
                    string const & colFieldName,
                    Connector connectivity,
                    std::vector< string > const & regions = {},
                    bool symmetric = true );

  /**
   * @copydoc addCoupling(string const &, string const &, Connector, std::vector< string > const &, bool)
   *
   * Overload for arrayView1d<string> input used by physics solvers.
   */
  void addCoupling( string const & rowFieldName,
                    string const & colFieldName,
                    Connector connectivity,
                    arrayView1d< string const > const & regions,
                    bool symmetric = true );

  /**
   * @brief Special interface for self-connectivity through a stencil.
   * @param [in] fieldName name of the field (this is only for diagonal blocks)
   * @param [in] stencils a pointer to FluxApproximation storing the stencils
   *
   * The field must be defined on element support. The set of regions is taken
   * automatically from the field definition.
   */
  void addCoupling( string const & fieldName,
                    FluxApproximationBase const & stencils );

  /**
   * @brief Finish populating fields and apply appropriate dof renumbering.
   *
   * This function must be called after all field and coupling information has been added.
   * It adjusts DoF index arrays to account for presence of other fields (in a global monolithic fashion).
   *
   * @note After DofManager has been closed, new fields and coupling cannot be added, until
   *       @ref clear or @ref setMesh is called.
   *
   * @note After reorderByRank() is called, the meaning of FieldDescription::globalOffset changes from
   *       "global offset of field's block in a global field-wise ordered (block) system" to
   *       "global offset of field's block on current processor in a rank-wise ordered system".
   *       This meaning is consistent with its use throughout. For example, this is the row/col
   *       global offset used to insert the field's sparsity block into a global coupled system.
   */
  void reorderByRank();

  /**
   * @brief Check if string key is already being used.
   * @param name field key to check
   * @return flag true if exists
   */
  bool fieldExists( string const & name ) const;

  /**
   * @brief Return the key used to record the field in the DofManager.
   * @param [in] fieldName string the name of the field.
   * @return string indicating name of the field.
   */
  string const & getKey( string const & fieldName ) const;

  /**
   * @brief @return number of global dofs of a given field.
   * @param fieldName the name of the field
   */
  globalIndex numGlobalDofs( string const & fieldName ) const;

  /**
   * @brief @return total number of dofs across all fields and processors.
   */
  globalIndex numGlobalDofs() const;

  /**
   * @brief @return number of local dofs of a given field.
   * @param fieldName name of the field
   */
  localIndex numLocalDofs( string const & fieldName ) const;

  /**
   * @brief @return total number of dofs across all fields on current processor.
   */
  localIndex numLocalDofs() const;

  /**
   * @brief @return rank offset of given field, i.e. total number of dofs of this field on previous processors.
   * @param fieldName name of the field
   */
  globalIndex rankOffset( string const & fieldName ) const;

  /**
   * @brief @return rank offset of current processor, i.e. total number of dofs on previous processors.
   */
  globalIndex rankOffset() const;

  /**
   * @brief @return number of components in a given field.
   * @param fieldName name of the field
   */
  integer numComponents( string const & fieldName = "" ) const;

  /**
   * @brief @return number of dof components across all fields.
   */
  integer numComponents() const;

  /**
   * @brief Get the support location type of the field.
   * @param [in] fieldName name of the field
   * @return support location type
   */
  Location location( string const & fieldName ) const;

  /**
   * @brief @return The list of region names in field's domain.
   * @param fieldName name of the field
   */
  std::vector< string > const & regions( string const & fieldName ) const;

  /**
   * @brief @return global offset of field's block on current processor in the system matrix.
   * @param [in] fieldName name of the field.
   */
  globalIndex globalOffset( string const & fieldName ) const;

  /**
   * @brief Return an array of number of components per field, sorted by field registration order.
   * @return array of number of components
   */
  array1d< integer > numComponentsPerField() const;

  /**
   * @brief Fill a container with unique dof labels for each local dof.
   * @tparam CONTAINER type of container to fill
   * @param labels the container to fill
   *
   * The labels are assigned starting from zero in increasing order of field registration.
   * Labels are repeated for each support point of a field.
   */
  template< typename CONTAINER >
  void getLocalDofComponentLabels( CONTAINER & labels ) const
  {
    labels.resize( numLocalDofs() );
    typename CONTAINER::value_type labelStart = 0;
    auto it = labels.begin();
    for( FieldDescription const & field : m_fields )
    {
      localIndex const numComp = field.numComponents;
      localIndex const numSupp = field.numLocalDof / numComp;
      for( localIndex i = 0; i < numSupp; ++i, it += numComp )
      {
        std::iota( it, it + numComp, labelStart );
      }
      labelStart += numComp;
    }
  }

  /**
   * @brief Populate sparsity pattern of the entire system matrix.
   * @param [out] pattern the target sparsity pattern
   */
  void setSparsityPattern( SparsityPattern< globalIndex > & pattern ) const;

  /**
   * @brief Copy values from LA vectors to simulation data arrays.
   *
   * @tparam VECTOR type of LA vector
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  template< typename VECTOR >
  void copyVectorToField( VECTOR const & vector,
                          string const & srcFieldName,
                          string const & dstFieldName,
                          real64 scalingFactor,
                          CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Copy values from LA vectors to simulation data arrays.
   *
   * @param localVector source local vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  void copyVectorToField( arrayView1d< real64 const > const & localVector,
                          string const & srcFieldName,
                          string const & dstFieldName,
                          real64 scalingFactor,
                          CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Add values from LA vectors to simulation data arrays.
   *
   * @tparam VECTOR type of LA vector
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  template< typename VECTOR >
  void addVectorToField( VECTOR const & vector,
                         string const & srcFieldName,
                         string const & dstFieldName,
                         real64 scalingFactor,
                         CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Add values from LA vectors to simulation data arrays.
   *
   * @param localVector source local vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  void addVectorToField( arrayView1d< real64 const > const & localVector,
                         string const & srcFieldName,
                         string const & dstFieldName,
                         real64 scalingFactor,
                         CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Copy values from nodes to DOFs.
   *
   * @tparam VECTOR type of LA vector
   * @param vector target LA vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  template< typename VECTOR >
  void copyFieldToVector( VECTOR & vector,
                          string const & srcFieldName,
                          string const & dstFieldName,
                          real64 scalingFactor,
                          CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Copy values from simulation data arrays to vectors.
   *
   * @param localVector target LA vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  void copyFieldToVector( arrayView1d< real64 > const & localVector,
                          string const & srcFieldName,
                          string const & dstFieldName,
                          real64 scalingFactor,
                          CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Add values from a simulation data array to a DOF vector.
   *
   * @tparam VECTOR type of LA vector
   * @param vector target LA vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  template< typename VECTOR >
  void addFieldToVector( VECTOR & vector,
                         string const & srcFieldName,
                         string const & dstFieldName,
                         real64 scalingFactor,
                         CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Add values from a simulation data array to a DOF vector.
   *
   * @param localVector target vector
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param mask component selection mask
   */
  void addFieldToVector( arrayView1d< real64 > const & localVector,
                         string const & srcFieldName,
                         string const & dstFieldName,
                         real64 scalingFactor,
                         CompMask mask = CompMask( MAX_COMP, true ) ) const;

  /**
   * @brief Create a dof selection by filtering out excluded components
   * @param excluded a list of dof components to exclude
   * @return a vector of remaining dof components
   *
   * @note Removed components must not have repeats, and each entry must either have
   *       loComp = 0 or hiComp = numComponents(fieldName) (or both). In other words,
   *       filtered out components must not leave "holes" in DOFs.
   */
  std::vector< SubComponent >
  filterDofs( std::vector< SubComponent > const & excluded ) const;

  /**
   * @brief Populate this manager from another using a sub-selection of fields/components.
   * @param source source dof manager
   * @param selection selection of fields/components
   * @note this will also allocate new dof
   */
  void setupFrom( DofManager const & source,
                  std::vector< SubComponent > const & selection );

  /**
   * @brief Create a matrix that restricts vectors and matrices to a subset of DOFs
   * @tparam MATRIX type of matrix used for restrictor
   * @param selection a list of fields to select; each entry is a struct containing
   *                  the name of the field and low and high selected component indices
   * @param comm the MPI communicator to use in the operator
   * @param transpose if @p true, the transpose (prolongation) operator will be created
   * @param restrictor resulting operator
   *
   * @note Can only be called after reorderByRank(), since global DOF indexing is required
   *       for the restrictor to make sense.
   */
  template< typename MATRIX >
  void makeRestrictor( std::vector< SubComponent > const & selection,
                       MPI_Comm const & comm,
                       bool transpose,
                       MATRIX & restrictor ) const;

  /**
   * @brief Print the summary of declared fields and coupling.
   *
   * @param os output stream
   */
  void printFieldInfo( std::ostream & os = std::cout ) const;

private:

  /**
   * Field description
   */
  struct FieldDescription
  {
    string name;                   ///< field name
    string key;                    ///< string key for index array
    string docstring;              ///< documentation string
    std::vector< string > regions; ///< list of support region names
    Location location;             ///< support location
    integer numComponents = 1;     ///< number of vector components
    localIndex numLocalDof = 0;    ///< number of local rows
    globalIndex numGlobalDof = 0;  ///< number of global rows
    globalIndex blockOffset = 0;   ///< offset of this field's block in a block-wise ordered system
    globalIndex rankOffset = 0;    ///< field's first DoF on current processor (within its block, ignoring other fields)
    globalIndex globalOffset = 0;  ///< global offset of field's DOFs on current processor for multi-field problems
  };

  /**
   * Coupling description
   */
  struct CouplingDescription
  {
    Connector connector = Connector::None;  //!< geometric object defining dof connections
    std::vector< string > regions; //!< list of region names
    FluxApproximationBase const * stencils = nullptr; //!< pointer to flux stencils for stencil based connections
  };

  /**
   * @brief Get field index from string key
   */
  localIndex getFieldIndex( string const & name ) const;

  /**
   * @brief Compute and save dof offsets the field
   * @param fieldIndex index of the field
   */
  void computeFieldDimensions( localIndex fieldIndex );

  /**
   * @brief Create index array for the field
   * @param field the field descriptor
   */
  void createIndexArray( FieldDescription const & field );

  /**
   * @brief Remove an index array for the field
   * @param field the field descriptor
   */
  void removeIndexArray( FieldDescription const & field );

  /**
   * @brief Calculate or estimate the number of nonzero entries in each local row
   * @param rowLengths array of row lengths (values are be incremented, not overwritten)
   * @param rowFieldIndex index of row field (must be non-negative)
   * @param colFieldIndex index of col field (must be non-negative)
   */
  void countRowLengthsOneBlock( arrayView1d< localIndex > const & rowLengths,
                                localIndex rowFieldIndex,
                                localIndex colFieldIndex ) const;

  void countRowLengthsFromStencil( arrayView1d< localIndex > const & rowLengths,
                                   localIndex fieldIndex ) const;

  /**
   * @brief Populate the sparsity pattern for a coupling block between given fields.
   * @param pattern the sparsity to be filled
   * @param rowFieldIndex index of row field (must be non-negative)
   * @param colFieldIndex index of col field (must be non-negative)
   *
   * This private function is used as a building block by higher-level SetSparsityPattern()
   */
  void setSparsityPatternOneBlock( SparsityPatternView< globalIndex > const & pattern,
                                   localIndex rowFieldIndex,
                                   localIndex colFieldIndex ) const;

  void setSparsityPatternFromStencil( SparsityPatternView< globalIndex > const & pattern,
                                      localIndex fieldIndex ) const;

  template< int DIMS_PER_DOF >
  void setFiniteElementSparsityPattern( SparsityPattern< globalIndex > & pattern,
                                        localIndex fieldIndex ) const;

  /**
   * @brief Generic implementation for @ref copyVectorToField and @ref addVectorToField
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param vector source LA vector
   * @param srcFieldName name of the source field (as defined in DofManager)
   * @param scalingFactor a factor to scale vector values by
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param dstFieldName name of the destination field (view wrapper key on the manager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename FIELD_OP, typename POLICY, typename LOCAL_VECTOR >
  void vectorToField( LOCAL_VECTOR localVector,
                      string const & srcFieldName,
                      string const & dstFieldName,
                      real64 scalingFactor,
                      CompMask mask ) const;

  /**
   * @brief Generic implementation for @ref copyFieldToVector and @ref addFieldToVector
   * @tparam FIELD_OP operation to perform (see FieldSpecificationOps.hpp)
   * @tparam POLICY execution policy for the kernel
   * @param manager mesh object manager that contains the target field (subregion for elements)
   * @param srcFieldName name of the source field (view wrapper key on the manager)
   * @param scalingFactor a factor to scale vector values by
   * @param vector ponter to target vector local data (host or device)
   * @param dstFieldName name of the destination field (as defined in DofManager)
   * @param loCompIndex index of starting DoF component (for partial copy)
   * @param hiCompIndex index past the ending DoF component (for partial copy)
   *
   * @note [@p loCompIndex , @p hiCompIndex) form a half-open interval.
   *       Negative value of @p hiCompIndex means use full number of field components
   */
  template< typename FIELD_OP, typename POLICY, typename LOCAL_VECTOR >
  void fieldToVector( LOCAL_VECTOR localVector,
                      string const & srcFieldName,
                      string const & dstFieldName,
                      real64 scalingFactor,
                      CompMask mask ) const;

  /// Name of the manager (unique, for unique identification of index array keys)
  string m_name;

  /// Pointer to corresponding MeshLevel
  MeshLevel * m_mesh = nullptr;

  /// Array of field descriptions
  std::vector< FieldDescription > m_fields;

  /// Table of connector types within and between fields
  std::map< std::pair< localIndex, localIndex >, CouplingDescription > m_coupling;

  /// Flag indicating that DOFs have been reordered rank-wise.
  bool m_reordered;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_DOFMANAGER_HPP_*/
