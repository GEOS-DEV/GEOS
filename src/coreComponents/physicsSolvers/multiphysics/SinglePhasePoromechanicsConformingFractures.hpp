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

/**
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"

namespace geos
{

template< typename FLOW_SOLVER = SinglePhaseBase >
class SinglePhasePoromechanicsConformingFractures : public PoromechanicsConformingFractures< SinglePhasePoromechanics, FLOW_SOLVER >
{
public:

  using Base = PoromechanicsConformingFractures< SinglePhasePoromechanics, FLOW_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_maxFaceNodes;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFractures"; }

  /**
   * @brief main constructor for SinglePhasePoromechanicsConformingFractures objects
   * @param name the name of this instantiation of SinglePhasePoromechanicsConformingFractures in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanicsConformingFractures
   */
  SinglePhasePoromechanicsConformingFractures( const string & name,
                                               dataRepository::Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanicsConformingFractures() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFractures object through the object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, SinglePhaseBase > )
    {
      return "SinglePhasePoromechanicsConformingFractures";
    }
    else
    {
      return FLOW_SOLVER::catalogName() + "PoromechanicsConformingFractures";
    }
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  GEOS_MGR_STRATEGY_NOT_SUPPORTED()

  virtual void assembleSystem( real64 const time_n,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override
  { Base::assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs ); }

  /**@}*/

protected:

  virtual void assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                                   string_array const & regionNames,
                                                                   DofManager const & dofManager,
                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                   arrayView1d< real64 > const & localRhs ) override;

  virtual integer numFluidComponents() const override { return 1; }

  virtual string getFlowDofKey() const override { return SinglePhaseBase::viewKeyStruct::elemDofFieldString(); }

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_ */
