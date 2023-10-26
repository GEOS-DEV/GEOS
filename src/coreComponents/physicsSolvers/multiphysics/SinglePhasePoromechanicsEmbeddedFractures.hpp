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
 * @file SinglePhasePoromechanicsEmbeddedFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"

namespace geos
{

class SinglePhasePoromechanicsEmbeddedFractures : public SinglePhasePoromechanics
{
public:
  SinglePhasePoromechanicsEmbeddedFractures( const std::string & name,
                                             Group * const parent );
  ~SinglePhasePoromechanicsEmbeddedFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsEmbeddedFractures object through the object
   * catalog.
   */
  static string catalogName() { return "SinglePhasePoromechanicsEmbeddedFractures"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override final;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  /**
   * @Brief add extra nnz to each row induced by the coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLengths the number of NNZ of each row
   */
  void addCouplingNumNonzeros( DomainPartition & domain,
                               DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addCouplingSparsityPattern( DomainPartition const & domain,
                                   DofManager const & dofManager,
                                   SparsityPatternView< globalIndex > const & pattern ) const;


  virtual void updateState( DomainPartition & domain ) override final;


  struct viewKeyStruct : SinglePhasePoromechanics::viewKeyStruct
  {
    constexpr static char const * fracturesSolverNameString() { return "fracturesSolverName"; }

    constexpr static char const * dTraction_dPressureString() { return "dTraction_dPressure"; }
  };


protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename EFEM_KERNEL_WRAPPER >
  real64 assemblyLaunch( MeshLevel & mesh,
                         DofManager const & dofManager,
                         arrayView1d< string const > const & regionNames,
                         string const & materialNamesString,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         real64 const & dt );

  string m_fracturesSolverName;

  SolidMechanicsEmbeddedFractures * m_fracturesSolver;

};


template< typename CONSTITUTIVE_BASE,
          typename KERNEL_WRAPPER,
          typename EFEM_KERNEL_WRAPPER >
real64 SinglePhasePoromechanicsEmbeddedFractures::assemblyLaunch( MeshLevel & mesh,
                                                                  DofManager const & dofManager,
                                                                  arrayView1d< string const > const & regionNames,
                                                                  string const & materialNamesString,
                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                  arrayView1d< real64 > const & localRhs,
                                                                  real64 const & dt )
{
  GEOS_MARK_FUNCTION;

  NodeManager const & nodeManager = mesh.getNodeManager();

  ElementRegionManager const & elemManager = mesh.getElemManager();
  SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( m_fracturesSolver->getFractureRegionName() );
  EmbeddedSurfaceSubRegion const & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  string const dofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
  string const jumpDofKey = dofManager.getKey( fields::contact::dispJump::key() );
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
  arrayView1d< globalIndex const > const & jumpDofNumber = subRegion.getReference< globalIndex_array >( jumpDofKey );

  string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  KERNEL_WRAPPER kernelWrapper( dispDofNumber,
                                dofManager.rankOffset(),
                                localMatrix,
                                localRhs,
                                dt,
                                gravityVectorData,
                                flowDofKey,
                                FlowSolverBase::viewKeyStruct::fluidNamesString() );

  real64 const maxForce =
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< >,
                                    CONSTITUTIVE_BASE,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            m_fracturesSolver->getSolidSolver()->getDiscretizationName(),
                                                            materialNamesString,
                                                            kernelWrapper );

  EFEM_KERNEL_WRAPPER EFEMkernelWrapper( subRegion,
                                         dispDofNumber,
                                         jumpDofNumber,
                                         flowDofKey,
                                         dofManager.rankOffset(),
                                         localMatrix,
                                         localRhs,
                                         dt,
                                         gravityVectorData,
                                         FlowSolverBase::viewKeyStruct::fluidNamesString() );

  finiteElement::
    regionBasedKernelApplication< parallelDevicePolicy< >,
                                  CONSTITUTIVE_BASE,
                                  CellElementSubRegion >( mesh,
                                                          regionNames,
                                                          m_fracturesSolver->getSolidSolver()->getDiscretizationName(),
                                                          materialNamesString,
                                                          EFEMkernelWrapper );

  return maxForce;

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_ */
