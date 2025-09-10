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
  const SurfaceElementRegion* surface = nullptr;
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

  void createFaceTypeListMortar( MortarSide side );

  

  // bubbles only on the slave side
  // void createBubbleCellList( MeshLevel & meshSlave,
  //                            FaceElementSubRegion const & surfRegion) const;

  void createBubbleCellList( ) const;


  void setMortarSurfaces( DomainPartition & domain);

  
  // call templated lambda with compile time knowledge of the element shapes on the master and the slave side
  // template<typename FUNC>
  // void forMortarSurfaces(FUNC&& func) const
  // {
  //   std::map<ElementShape, array1d<localIndex>> const&
  //   faceTypesToFaceElementsSlave = m_faceTypeToElementList.at(MortarSide::Slave);

  //   std::map<ElementShape, array1d<localIndex>> const&
  //       faceTypesToFaceElementsMaster = m_faceTypeToElementList.at(MortarSide::Master);

  //   for (const auto& [slaveShape, slaveElementList] : faceTypesToFaceElementsSlave) {
  //       GEOS_UNUSED_VAR(slaveElementList);

  //       for (const auto& [masterShape, masterElementList] : faceTypesToFaceElementsMaster) 
  //       {
  //         GEOS_UNUSED_VAR(masterElementList);

  //         switch (slaveShape) 
  //         {
  //         case ElementShape::Triangle:
  //           switch (masterShape) 
  //           {
  //           case ElementShape::Triangle:
  //             func.template operator()<ElementShape::Triangle, ElementShape::Triangle>();
  //             break;
  //           case ElementShape::Quadrilateral:
  //             func.template operator()<ElementShape::Triangle, ElementShape::Quadrilateral>();
  //             break;
  //           }
  //           break;
  //         case ElementShape::Quadrilateral:
  //           switch (masterShape) 
  //           {
  //           case ElementShape::Triangle:
  //             func.template operator()<ElementShape::Quadrilateral, ElementShape::Triangle>();
  //             break;
  //           case ElementShape::Quadrilateral:
  //             func.template operator()<ElementShape::Quadrilateral, ElementShape::Quadrilateral>();
  //             break;
  //           }
  //           break;
  //         }
  //       }
  //   }
  // }

  // template< typename LAMBDA >
  // void forMortarSurfacesOld( LAMBDA && lambda ) const
  // {

  //   std::map< ElementShape, array1d< localIndex > > const &
  //   faceTypesToFaceElementsSlave = m_faceTypeToElementList.at( MortarSide::Slave );

  //   std::map< ElementShape, array1d< localIndex > > const &
  //   faceTypesToFaceElementsMaster = m_faceTypeToElementList.at( MortarSide::Master );

  //   for( const auto & [slaveShape, slaveElementList] : faceTypesToFaceElementsSlave )
  //   {
  //     GEOS_UNUSED_VAR( slaveElementList );
  //     for(const auto & [masterShape, masterElementList] : faceTypesToFaceElementsMaster )
  //     {
  //       GEOS_UNUSED_VAR( masterElementList );
  //       finiteElement::FiniteElementBase const & slaveFE = *(m_faceTypeToMortarFiniteElements.at( slaveShape ));
  //       finiteElement::FiniteElementBase const & masterFE = *(m_faceTypeToMortarFiniteElements.at( masterShape ));

  //       dispatchMortar< MORTAR_FE_TYPES >(  slaveFE, masterFE, 
  //       [&](const auto& finiteElementSlave, const auto& finiteElementMaster)
  //       {
  //         lambda(finiteElementSlave, finiteElementMaster, slaveShape, masterShape);
  //       });
  //     }
  //   }
  // }


  // template< typename LAMBDA >
  // void forMortarSurfaces( LAMBDA && lambda ) const
  // {

  //   std::map< string, array1d< localIndex > > const &
  //   faceTypesToFaceElementsSlave = m_faceTypesToElementList.at( MortarSide::Slave );

  //   std::map< string, array1d< localIndex > > const &
  //   faceTypesToFaceElementsMaster = m_faceTypesToElementList.at( MortarSide::Master );

  //   for( const auto & [slaveShape, slaveElementList] : faceTypesToFaceElementsSlave )
  //   {
  //     for(const auto & [masterShape, masterElementList] : faceTypesToFaceElementsMaster )
  //     {
  //       arrayView1d< localIndex const > const slaveElemList = slaveElementList.toViewConst();
  //       arrayView1d< localIndex const > const masterElemList = masterElementList.toViewConst();

  //       finiteElement::FiniteElementBase const & slaveFE = *(m_faceTypeToFiniteElements.at( slaveShape ));
  //       finiteElement::FiniteElementBase const & masterFE = *(m_faceTypeToFiniteElements.at( masterShape ));

  //       lambda(slaveShape, masterShape, slaveElemList, masterElemList);
  //     }
  //   }
  // }


  // template< typename FETypes, typename LAMBDA >
  // void dispatchMortar( const finiteElement::FiniteElementBase & slaveFE, const finiteElement::FiniteElementBase & masterFE, LAMBDA&& lambda )
  // {
  //   geos::finiteElement::FiniteElementDispatchHandler< FETypes >::dispatch2D( slaveFE,
  //   [&]( const auto & finiteElementSlave )
  //   {
  //     geos::finiteElement::FiniteElementDispatchHandler< FETypes >::dispatch2D( masterFE,
  //     [&]( const auto & finiteElementMaster )
  //     {
  //       lambda(finiteElementSlave, finiteElementMaster);
  //     } );
  //   } );
  // } 


private:

  // Struct to hold the pointers for a single contact surface

  // store pointers to relevant mesh objects of master and slave side
  MeshLevel * m_meshSlave = nullptr;
  MeshLevel * m_meshMaster = nullptr;
  FaceElementSubRegion const* m_surfaceMaster = nullptr;     
  FaceElementSubRegion const* m_surfaceSlave = nullptr;


   std::map< MortarSide, std::map< ElementShape, array1d< localIndex > > > m_faceTypeToElementList;

   // holds pointers to master and slave MeshLevel and SurfaceElementRegion
   std::map< MortarSide, MortarSurface > m_mortarSide;

   // given a certain mortar pair of shapes, holds the id of the element connected by each mortar
   // subtriangle and the corresponding jacobian determinant. The kernel will use these lists
   std::map< std::pair< ElementShape, ElementShape >, array2d< localIndex> > m_triCells;
   std::map< std::pair< ElementShape, ElementShape >, array1d< real64 > >   m_triCellsDet;

   // hold local coordinates of gauss points of each triangle subcell on each mortar side
   std::map< MortarSide, std::map< std::pair< ElementShape, ElementShape >, array2d< real64 > > > m_gpLocalCoords;

  // coordinates of gauss points ready to go
  constexpr static localIndex nGPtri = feTriangleCell::numQuadraturePoints;
  constexpr static real64 qCoords[nGPtri][2] = {
    { 0.333333333333333, 0.333333333333333 },
    { 0.600000000000000, 0.200000000000000 },
    { 0.200000000000000, 0.600000000000000 },
    { 0.200000000000000, 0.200000000000000 }
  };




  // map id of slave elements to id of connected master elements
  //ArrayOfArrays< localIndex > m_connectivityMapSlave;  

  // map id of master elements to id of connected master elements
  //ArrayOfArrays< localIndex > m_connectivityMapMaster;  

  //void addCouplingNumNonzeros( DofManager & dofManager,
  //                            arrayView1d< localIndex > const & rowLengths ) const;

 
  //void addCouplingSparsityPattern( DofManager const & dofManager,
  //                                 SparsityPatternView< globalIndex > const & pattern ) const;

  /// Finite element type to face element index map
  FaceTypeMap m_faceTypesToFaceElementsMaster;
  FaceTypeMap m_faceTypesToFaceElementsSlave;

  string m_slaveName;

  string m_masterName;

  /// Finite element type to finite element object map
  //std::map< string, std::unique_ptr< geos::finiteElement::FiniteElementBase > > m_faceTypeToFiniteElements;

  std::map< ElementShape, std::unique_ptr< geos::finiteElement::FiniteElementBase > > m_faceTypeToMortarFiniteElements;

  /// Map gauss point list to master basis functions values
  std::map< string, ArrayOfArrays< real64 > > m_gpToMasterBasis; 

  

  /// Map gauss point list to corresponding master element id 
  std::map< string, array2d< localIndex > > m_gpToMasterId; 

  /// List of all existing pairs of slave-master elements (derived from m_gpToMasterId)
  std::map< string, array2d< localIndex > > m_mortarPairs;

  struct viewKeyStruct : ContactSolverBase::viewKeyStruct
  {

    constexpr static char const * masterString() { return "master"; }

    constexpr static char const * slaveString() { return "slave"; }

  };

  ///
  //void computeMortarInterpolation();

  template< ElementShape slaveShape, ElementShape masterShape >
  void computeMortarInterpolationNew();

  template< ElementShape slaveShape, ElementShape masterShape >
  void processMortarPair( localIndex const slaveFaceId, 
                          localIndex const masterFaceId,
                          arraySlice1d< localIndex const > const & nodesSlave,
                          arraySlice1d< localIndex const > const & nodesMaster, 
                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & coordsSlave,
                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & coordsMaster );


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

  void getLocalInterpolationPoints( localIndex nInt, string const & finiteElementName, array2d< real64 > & localCoordsMaster );

  real64 computeRBF( real64 const d, real64 const radius );

  real64 computeRBF( real64 const (&pt1)[3], real64 const (&pt2)[3], real64 const radius );

  void computeRBFweights( arrayView2d<real64> M, arrayView2d<real64> Nm, arrayView2d<real64> wRBF, localIndex nIntPts, localIndex numNodeperElement );

  real64 computeRBFmatrix( arrayView2d<real64> realCoordsInterpolation, arrayView2d<real64> RBFmatrix, localIndex nIntPts);
  
  /// Tandem traversal contact search 
  void contactSearch(std::unique_ptr<TreeNodeMortar> const & nodeMaster,
                     std::unique_ptr<TreeNodeMortar> const & nodeSlave,
                     ArrayOfArrays<localIndex> & connectivityMap);

  /// Check intersection between two bounding boxes using polytops primitives                   
  bool checkIntersection(std::unique_ptr<TreeNodeMortar> const & nodeMaster,
                         std::unique_ptr<TreeNodeMortar> const & nodeSlave);

  
  // compute connectivityMap and return total number of connections
  localIndex getConnectivityMap( ElementShape slaveShape, ElementShape masterShape, ArrayOfArrays<localIndex> & connectivityMap);

  template<ElementShape S>
  decltype(auto) getFE() 
  {
    auto & femTypePtr = m_faceTypeToMortarFiniteElements.at(S); // unique_ptr<FiniteElementBase>

    if constexpr (S == ElementShape::Triangle) 
    {
        using femType = finiteElement::H1_TriangleFace_Lagrange1_Gauss6;
        return *static_cast<femType*>(femTypePtr.get()); 
    } 
    else if constexpr (S == ElementShape::Quadrilateral) 
    {
        using femType = finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6;
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
