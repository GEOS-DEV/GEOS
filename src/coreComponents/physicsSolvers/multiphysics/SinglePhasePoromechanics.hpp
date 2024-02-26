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
 * @file SinglePhasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICS_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsSolver.hpp"


namespace geos
{

template< typename FLOW_SOLVER >
class SinglePhasePoromechanics : public PoromechanicsSolver< FLOW_SOLVER >
{
public:

  using Base = PoromechanicsSolver< FLOW_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /**
   * @brief main constructor for SinglePhasePoromechanics objects
   * @param name the name of this instantiation of SinglePhasePoromechanics in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanics
   */
  SinglePhasePoromechanics( const string & name,
                            dataRepository::Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanics() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanics object through the object catalog.
   */
  static string catalogName();
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void postProcessInput() override;

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void updateState( DomainPartition & domain ) override;

  /**@}*/

  struct viewKeyStruct : Base::viewKeyStruct
  {
    // nothing yet here
  };

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  void assembleElementBasedTerms( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

private:

  /**
   * @brief Helper function to recompute the bulk density
   * @param[in] subRegion the element subRegion
   */
  void updateBulkDensity( ElementSubRegionBase & subRegion );

  void createPreconditioner();

  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  real64 assemblyLaunch( MeshLevel & mesh,
                         DofManager const & dofManager,
                         arrayView1d< string const > const & regionNames,
                         string const & materialNamesString,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         real64 const dt,
                         PARAMS && ... params );

};

template< typename FLOW_SOLVER >
template< typename CONSTITUTIVE_BASE,
          typename KERNEL_WRAPPER,
          typename ... PARAMS >
real64 SinglePhasePoromechanics< FLOW_SOLVER >::assemblyLaunch( MeshLevel & mesh,
                                                                DofManager const & dofManager,
                                                                arrayView1d< string const > const & regionNames,
                                                                string const & materialNamesString,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs,
                                                                real64 const dt,
                                                                PARAMS && ... params )
{
  GEOS_MARK_FUNCTION;

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const dofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
  arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( SolverBase::gravityVector() );

  KERNEL_WRAPPER kernelWrapper( dofNumber,
                                dofManager.rankOffset(),
                                localMatrix,
                                localRhs,
                                dt,
                                gravityVectorData,
                                std::forward< PARAMS >( params )... );

  return finiteElement::
           regionBasedKernelApplication< parallelDevicePolicy< >,
                                         CONSTITUTIVE_BASE,
                                         CellElementSubRegion >( mesh,
                                                                 regionNames,
                                                                 this->solidMechanicsSolver()->getDiscretizationName(),
                                                                 materialNamesString,
                                                                 kernelWrapper );
}


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICS_HPP_ */
