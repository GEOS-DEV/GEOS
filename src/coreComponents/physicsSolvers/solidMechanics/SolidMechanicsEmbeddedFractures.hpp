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

/*
 * SolidMechanicsEmbeddedFractures.hpp
 *
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
using namespace constitutive;

class SolidMechanicsLagrangianFEM;

class SolidMechanicsEmbeddedFractures : public SolverBase
{
public:
  SolidMechanicsEmbeddedFractures( const std::string & name,
                                   Group * const parent );

  ~SolidMechanicsEmbeddedFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "SolidMechanicsEmbeddedFractures";
  }

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

  virtual void SetupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void SetupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            array1d< real64 > & localRhs,
                            array1d< real64 > & localSolution,
                            bool const setSparsity = true ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void ImplicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;


  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition & domain ) override final;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  template< typename CONSTITUTIVE_BASE,
            template< typename SUBREGION_TYPE,
                      typename CONSTITUTIVE_TYPE,
                      typename FE_TYPE > class KERNEL_TEMPLATE,
            typename ... PARAMS >
  void AssemblyLaunch( DomainPartition & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs,
                       PARAMS && ... params );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto solidSolverNameString = "solidSolverName";

    constexpr static auto contactRelationNameString = "contactRelationName";

    constexpr static auto dispJumpString = "displacementJump";

    constexpr static auto deltaDispJumpString = "deltaDisplacementJump";

    constexpr static auto fractureRegionNameString = "fractureRegionName";

  } SolidMechanicsEmbeddedFracturesViewKeys;

protected:

  void AddCouplingNumNonzeros( DomainPartition & domain,
                               DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the perforations
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void AddCouplingSparsityPattern( DomainPartition const & domain,
                                   DofManager const & dofManager,
                                   SparsityPatternView< globalIndex > const & pattern ) const;

  /*
   * @brief Assemble Equilibrium operator
   * @param eqMatrix Equilibrium operator
   * @param embeddedSurfaceSubRegion subRegion
   * @param k cell index
   * @param hInv scaling coefficient
   */
  void AssembleEquilibriumOperator( array2d< real64 > & eqMatrix,
                                    EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion,
                                    const localIndex k,
                                    const real64 hInv );
  /*
   * @brief Assemble Compatibility operator
   * @param compMatrix
   * @param embeddedSurfaceSubRegion
   * @param k cell index
   * @param q quadrature point index
   * @param elemsToNodes element to node map
   * @param nodesCoord nodes coordinates
   * @param embeddedSurfaceToCell embedded surface to cell maps
   * @param numNodesPerElement number of nodes per element
   * @param dNdX shape functions derivatives
   */
  void AssembleCompatibilityOperator( array2d< real64 > & compMatrix,
                                      EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion,
                                      localIndex const k,
                                      localIndex const q,
                                      CellBlock::NodeMapType const & elemsToNodes,
                                      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord,
                                      arrayView1d< localIndex const > const & embeddedSurfaceToCell,
                                      localIndex const numNodesPerElement,
                                      arrayView4d< real64 const > const & dNdX );

  /*
   * @brief Assemble Compatibility operator
   * @param strainMatrix strain matrix (B)
   * @param elIndex element index
   * @param q quadrature point index
   * @param numNodesPerElement number of nodes per element
   * @param dNdX shape functions derivatives
   */
  void AssembleStrainOperator( array2d< real64 > & strainMatrix,
                               localIndex const elIndex,
                               localIndex const q,
                               localIndex const numNodesPerElement,
                               arrayView4d< real64 const > const & dNdX );
  /*
   * @brief Computes traction and derivative on each fracture segment.
   * @param constitutiveManager constant pointer to the constitutive mamanger
   * @param dispJump displacement jump
   * @param tractionVector traction vector
   * @param dTdw Derivative of the traction w.r.t. the jump.
   */
  void ComputeTraction( ConstitutiveManager const * const constitutiveManager,
                        array1d< real64 >  const & dispJump,
                        array1d< real64 > & tractionVector,
                        array2d< real64 > & dTdw );


private:

  /// Solid mechanics solver name
  string m_solidSolverName;

  /// fracture region name
  string m_fractureRegionName;

  /// pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// contact relation name string
  string m_contactRelationName;

};


//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************


template< typename CONSTITUTIVE_BASE,
          template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    typename FE_TYPE > class KERNEL_TEMPLATE,
          typename ... PARAMS >
void SolidMechanicsEmbeddedFractures::AssemblyLaunch( DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs,
                                                      PARAMS && ... params )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = *(domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 ));

  NodeManager const & nodeManager = *(mesh.getNodeManager());
  ElementRegionManager const & elemManager = *(mesh.getElemManager());
  EmbeddedSurfaceRegion * const region = elemManager.GetRegion< EmbeddedSurfaceRegion >( m_fractureRegionName );
  EmbeddedSurfaceSubRegion * const subRegion = region->GetSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  string const dispDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  arrayView1d< globalIndex const > const & jumpDofNumber = subRegion->getReference< globalIndex_array >( jumpDofKey );

  ResetStressToBeginningOfStep( domain );

  real64 const gravityVectorData[3] = { gravityVector().Data()[0],
                                        gravityVector().Data()[1],
                                        gravityVector().Data()[2] };

  real64 maxTraction = finiteElement::
		  regionBasedKernelApplication< parallelDevicePolicy< 32 >,
		                                CONSTITUTIVE_BASE,
		  		                        CellElementSubRegion,
		  		                        KERNEL_TEMPLATE >( mesh,
		  				                targetRegionNames(),
		  				                this->getDiscretizationName(),
		  				                m_solidMaterialNames,
										subRegion,
		  				                dispDofNumber,
		  				                jumpDofNumber,
		  				                dofManager.rankOffset(),
		  				                localMatrix,
		  				                localRhs,
		  				                gravityVectorData,
		  				                std::forward< PARAMS >( params )... );

}

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_ */
