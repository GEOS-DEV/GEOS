/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2    std::map< ElementShape, array1d< localIndex > > const &
    faceTypesToFaceElementsSlave = m_faceTypeToElementList.at( MortarSide::Slave );

    std::map< ElementShape, array1d< localIndex > > const &
    faceTypesToFaceElementsMaster = m_faceTypeToElementList.at( MortarSide::Master );

    for( const auto & [slaveShape, slaveElementList] : faceTypesToFaceElementsSlave )
    {
      GEOS_UNUSED_VAR( slaveElementList );
      for(const auto & [masterShape, masterElementList] : faceTypesToFaceElementsMaster )
      {
        GEOS_UNUSED_VAR( masterElementList );
        finiteElement::FiniteElementBase const & slaveFE = *(m_faceTypeToFiniteElements.at( slaveShape ));
        finiteElement::FiniteElementBase const & masterFE = *(m_faceTypeToFiniteElements.at( masterShape ));Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
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
#include "mesh/utilities/ComputationalGeometry.hpp"
//#include "finiteElement/FiniteElementDispatch.hpp"

namespace geos
{

class TreeNodeMortar;

enum class ElementShape { Triangle, Quadrilateral };
enum class MortarSide { Slave, Master };

using feTriangleCell = finiteElement::H1_TriangleFace_Lagrange1_Gauss4;

struct MortarSurface
{
  MeshLevel* mesh = nullptr;
  SurfaceElementRegion* surface = nullptr;
};

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

  //using FaceTypeMap = std::map< string, std::map< string, array1d< localIndex > > >;

  using connectivityMapType = std::map< std::pair< ElementShape, ElementShape >, ArrayOfArrays < localIndex > >;

 
  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  void registerMortarDataOnMesh( );

  void setupDofs( DomainPartition const & domain,
                  DofManager & dofManager) const;

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

  real64 calculateContactResidualNorm( DomainPartition const & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< real64 const > const & localRhs );

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    real64 const dt,
                                    DomainPartition & domain ) override;

  void updateState( DomainPartition & domain ) override final;

  void createFaceTypeListMortar( MortarSide side );

  

  void createBubbleCellList( ) const;


  void setMortarSurfaces( DomainPartition & domain);

  

private:

  string m_meshSlaveName;

  std::map< MortarSide, std::map< ElementShape, array1d< localIndex > > > m_faceTypeToElementList;

  // holds pointers to master and slave MeshLevel and SurfaceElementRegion
  std::map< MortarSide, MortarSurface > m_mortarSide;

  // list of pairs of slave/master element ids for each mortar subtriangle
  std::map< MortarSide, std::map< std::pair< ElementShape, ElementShape >, array1d< localIndex > > > m_triCells;

  // list of determinant of each mortar subtriangle
  std::map< std::pair< ElementShape, ElementShape >, array1d< real64 > >  m_triCellsDet;

  // list local coordinates of gauss points of subtriangles on each mortar side
  std::map< MortarSide, std::map< std::pair< ElementShape, ElementShape >, array3d< real64 > > > m_gpLocalCoords;

  // coordinates of triangle gauss points ready to go (they are private in the triangle class)
  constexpr static localIndex nGPtri = feTriangleCell::numQuadraturePoints;
  constexpr static real64 qCoords[nGPtri][2] = {
    { 0.333333333333333, 0.333333333333333 },
    { 0.600000000000000, 0.200000000000000 },
    { 0.200000000000000, 0.600000000000000 },
    { 0.200000000000000, 0.200000000000000 }
  };


  void addBubbleCouplingNumNonzeros( DofManager & dofManager,
                                     arrayView1d< localIndex > const & rowLengths ) const;

  void addMortarCouplingNumNonzeros( DofManager & dofManager,
                                     ElementShape const & slaveShape,
                                     ElementShape const & masterShape,
                                     connectivityMapType const & connectivityMap,
                                     arrayView1d< localIndex > const & rowLengths ) const;

  void addBubbleCouplingSparsityPattern( DofManager const & dofManager,
                                         SparsityPatternView< globalIndex > const & pattern ) const;

  void addMortarCouplingSparsityPattern( DofManager & dofManager,
                                         ElementShape const & slaveShape,
                                         ElementShape const & masterShape,
                                         SparsityPatternView< globalIndex > const & pattern ) const;

  void assembleBubbles( real64 const dt,
                        DomainPartition & domain,
                        DofManager const & dofManager,
                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                        arrayView1d< real64 > const & localRhs );
  
  template< ElementShape shape>
  void assembleMortarBubbles( DofManager const & dofManager,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs );

  void assembleMortar( real64 const dt,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs );

  template< ElementShape shpS, ElementShape shpM >
  void assembleMortar( real64 const dt,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs );

  void computeRotationMatrices( );



  string m_slaveName;
  string m_masterName;

  std::map< ElementShape, std::unique_ptr< geos::finiteElement::FiniteElementBase > > m_faceTypeToMortarFiniteElements;

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {

    constexpr static char const * masterString() { return "master"; }

    constexpr static char const * slaveString() { return "slave"; }

  };

  void computeMortarInterpolation ( connectivityMapType & connectivityMap);

  template< ElementShape slaveShape, ElementShape masterShape >
  void computeMortarInterpolation( ArrayOfArrays<localIndex> const & connections );

  template< ElementShape slaveShape, ElementShape masterShape >
  localIndex processMortarPair( localIndex const slaveFaceId, 
                                localIndex const masterFaceId,
                                arraySlice1d< localIndex const > const & nodesSlave,
                                arraySlice1d< localIndex const > const & nodesMaster, 
                                arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & coordsSlave,
                                arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & coordsMaster,
                                arrayView2d< localIndex > & cellPairs,
                                arrayView1d< real64 > & subTriDeterminants,
                                arrayView3d< real64 > & localCoordsSlave,
                                arrayView3d< real64 > & localCoordsMaster,
                                localIndex const & kPair );
      

  template< MortarSide side >                         
  void projectPointInPlane( real64 const (& coord3d)[3],
                            real64 const (& normal)[3],
                            real64 const (& origin)[3],
                            real64 (& proj2d)[2]);

  template< localIndex sizePoly, localIndex sizeClipper>
  void polygonClipping( array2d<real64> & poly, array2d<real64> & clipPoly);

  void intersect(real64 x1, real64 y1,real64 x2, real64 y2,
                 real64 x3, real64 y3,real64 x4, real64 y4,
                 real64 & xInt, real64 & yInt);

  void clip( array2d<real64> & poly,
             real64 xc1, real64 yc1, real64 xc2, real64 yc2);

  bool validateClip( array2d< real64 > & clipPoly);

  template< ElementShape shape >
  void projectGP( real64 const (& coordsTri)[3][2],
                  arrayView2d<real64 const> const & coordsElem,
                  real64 (& xi)[nGPtri][2]);

  template<ElementShape shape>
  bool checkInFE(real64 xi0, real64 xi1);

  template<localIndex numNodeElement>
  void calcGradN( real64 const (& xi)[2], real64 (& dN)[2][numNodeElement] );

  template<localIndex numNodes>
  void permuteN(real64 (& N)[numNodes]);

  void getConnectivityMap( connectivityMapType & connectivityMap );

  void getMortarConnections( ElementShape slaveShape, ElementShape masterShape, ArrayOfArrays<localIndex> & connections);

  // Tandem traversal contact search between master and slave nodes
  void contactSearch(std::unique_ptr<TreeNodeMortar> const & nodeMaster,
                     std::unique_ptr<TreeNodeMortar> const & nodeSlave,
                     ArrayOfArrays<localIndex> & connections);

  // Check intersection between two bounding boxes using polytops primitives
  bool checkIntersection(std::unique_ptr<TreeNodeMortar> const & nodeMaster,
                         std::unique_ptr<TreeNodeMortar> const & nodeSlave);

  // retrieve finite element type from templated element shape
  template<ElementShape S>
  decltype(auto) getFE() 
  {
    auto & femTypePtr = m_faceTypeToMortarFiniteElements.at(S); // unique_ptr<FiniteElementBase>

    if constexpr (S == ElementShape::Triangle) 
    {
        using femType = finiteElement::H1_TriangleFace_Lagrange1_Gauss4;
        
        return *static_cast<femType*>(femTypePtr.get()); 
    } 
    else if constexpr (S == ElementShape::Quadrilateral) 
    {
        using femType = finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2;
        return *static_cast<femType*>(femTypePtr.get()); 
    } 
    else 
    {
        static_assert(S == ElementShape::Triangle || S == ElementShape::Quadrilateral,
                      "Unsupported ElementShape");
    }
  }

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
    void createNode( MeshLevel const & mesh, 
                     FaceElementSubRegion const & surf, 
                     arrayView1d<localIndex> & surfId, 
                     array1d<localIndex> & surfList );

  };
  

}; /* namespace geos */



#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_ */
