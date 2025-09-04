/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * SolidMechanicsMortarContact.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_

#define POLYTOP { \
    1e100, -1e100, 1e100, -1e100, 1e100, -1e100, 1e100, -1e100, 1e100, \
    -1e100, 1e100, -1e100, 1e100, -1e100, 1e100, -1e100, 1e100, -1e100  \
}

#define POLYTOP_PRIMITIVES { \
    {1, 0, 0, 1, 1, 0, 1, 1, 0}, \
    {0, 1, 0, 1, 0, 1, -1, 0, 1}, \
    {0, 0, 1, 0, 1, 1, 0, -1, -1} \
}

#define BOUNDING_BOX_EXPANSION 0.025


#include "physicsSolvers/solidMechanics/contact/ContactSolverBase.hpp"

namespace geos
{

class TreeNodeMortar;

class SolidMechanicsMortarContact : public ContactSolverBase
{
public:
  SolidMechanicsMortarContact( const string & name,
                                            Group * const parent );

  ~SolidMechanicsMortarContact() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SolidMechanicsMortarContact";
  }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override {return catalogName(); }

  using FaceTypeMap = std::map< string, std::map< string, array1d< localIndex > > >;

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager,
                          string const & meshSlaveName ) const;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override final;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override final;

  virtual void implicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64 calculateResidualNorm( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs ) override;

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;


  FaceTypeMap createFaceTypeList( MeshLevel const & mesh, 
                                  FaceElementSubRegion const & surface,
                                  string name);

  // bubbles only on the slave side
  void createBubbleCellList( MeshLevel & meshSlave,
                             FaceElementSubRegion const & surfRegion) const;

 

  

private:

  // store pointers to relevant mesh objects of master and slave side
  MeshLevel * m_meshSlave = nullptr;
  MeshLevel * m_meshMaster = nullptr;
  FaceElementSubRegion const* m_surfaceMaster = nullptr;     
  FaceElementSubRegion const* m_surfaceSlave = nullptr;

  // map id of slave elements to id of connected master elements
  ArrayOfArrays< localIndex > m_connectivityMapSlave;  

  // map id of master elements to id of connected master elements
  ArrayOfArrays< localIndex > m_connectivityMapMaster;  

  void addCouplingNumNonzeros( DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

 
  void addCouplingSparsityPattern( DofManager const & dofManager,
                                   SparsityPatternView< globalIndex > const & pattern ) const;

  /// Finite element type to face element index map
  FaceTypeMap m_faceTypesToFaceElementsMaster;
  FaceTypeMap m_faceTypesToFaceElementsSlave;

  /// Finite element type to finite element object map
  std::map< string, std::unique_ptr< geos::finiteElement::FiniteElementBase > > m_faceTypeToFiniteElements;
  
  /// Map gauss point list to master basis functions values
  std::map< string, ArrayOfArrays< real64 > > m_gpToMasterBasis; 

  /// Map gauss point list to corresponding master element id 
  std::map< string, array2d< localIndex > > m_gpToMasterId; 

  /// List of all existing pairs of slave-master elements (derived from m_gpToMasterId)
  std::map< string, array2d< localIndex > > m_mortarPairs;

  ///
  void computeMortarInterpolation();

  void getLocalInterpolationPoints( localIndex nInt, string const & finiteElementName, array2d< real64 > & localCoordsMaster );

  real64 computeRBF( real64 const d, real64 const radius );

  real64 computeRBF( real64 const (&pt1)[3], real64 const (&pt2)[3], real64 const radius );

  void computeRBFweights( arrayView2d<real64> M, arrayView2d<real64> Nm, arrayView2d<real64> wRBF, localIndex nIntPts, localIndex numNodeperElement );

  real64 computeRBFmatrix( arrayView2d<real64> realCoordsInterpolation, arrayView2d<real64> RBFmatrix, localIndex nIntPts);
  
  /// Tandem traversal contact search 
  void contactSearch(std::unique_ptr<TreeNodeMortar> const & nodeMaster,
                     std::unique_ptr<TreeNodeMortar> const & nodeSlave);

  /// Check intersection between two bounding boxes using polytops primitives                   
  bool checkIntersection(std::unique_ptr<TreeNodeMortar> const & nodeMaster,
                         std::unique_ptr<TreeNodeMortar> const & nodeSlave);

  void getConnectivityMap();

};

// binary tree for efficient contact search
class TreeNodeMortar
  {
  public:
    std::unique_ptr<TreeNodeMortar> left = nullptr;
    std::unique_ptr<TreeNodeMortar> right = nullptr;
    bool isLeaf = false;

     // primitives of polytopal bounding box 
    double polytop[18] = POLYTOP; 

     // id of face corresponding to leaf nodes  
    localIndex leafId;  
    
    // populate a tree node (use recursion)
    void createNode(MeshLevel const & mesh, 
                    FaceElementSubRegion const & surf, 
                    array1d<localIndex> & surfList);


    
  };

}; /* namespace geos */



#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_ */
