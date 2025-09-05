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

#define MORTAR_FE_TYPES \
  finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre6, \
  finiteElement::H1_TriangleFace_Lagrange1_Gauss6


#include "physicsSolvers/solidMechanics/contact/ContactSolverBase.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"

namespace geos
{

class TreeNodeMortar;

enum class ElementShape { Triangle, Quadrilateral };
enum class MortarSide { Slave, Master };

using feTriangleCell = finiteElement::H1_TriangleFace_Lagrange1_Gauss6;

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

   std::map< MortarSide, MortarSurface > m_mortarSide;


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

  string m_slaveName;

  string m_masterName;

  /// Finite element type to finite element object map
  std::map< string, std::unique_ptr< geos::finiteElement::FiniteElementBase > > m_faceTypeToFiniteElements;

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
  void computeMortarInterpolation();

  template< ElementShape slaveShape, ElementShape masterShape >
  void computeMortarInterpolationNew();

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
    void createNode(MeshLevel const & mesh, 
                    FaceElementSubRegion const & surf, 
                    array1d<localIndex> & surfList);


    
  };

}; /* namespace geos */



#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSMORTARCONTACT_HPP_ */
