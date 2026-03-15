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
 * @file MultiphasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanics.hpp"

namespace geos
{

template< typename FLOW_SOLVER = CompositionalMultiphaseBase >
class MultiphasePoromechanicsConformingFractures : public PoromechanicsConformingFractures< MultiphasePoromechanics, FLOW_SOLVER >
{
public:

  using Base = PoromechanicsConformingFractures< MultiphasePoromechanics, FLOW_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_maxFaceNodes;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFractures"; }

  /**
   * @brief main constructor for MultiphasePoromechanicsConformingFractures objects
   * @param name the name of this instantiation of MultiphasePoromechanicsConformingFractures in the repository
   * @param parent the parent group of this instantiation of MultiphasePoromechanicsConformingFractures
   */
  MultiphasePoromechanicsConformingFractures( const string & name,
                                              dataRepository::Group * const parent );

  /// Destructor for the class
  ~MultiphasePoromechanicsConformingFractures() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFractures object through the object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, CompositionalMultiphaseBase > )
    {
      return "MultiphasePoromechanicsConformingFractures";
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

  void assembleSystem( real64 const time_n,
                       real64 const dt,
                       DomainPartition & domain,
                       DofManager const & dofManager,
                       PhysicsSolverBase::MATRIX_VIEW const & localMatrix,
                       arrayView1d< real64 > const & localRhs ) override;

protected:

  virtual void initializePreSubGroups() override;

  virtual void assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                                   string_array const & regionNames,
                                                                   DofManager const & dofManager,
                                                                   PhysicsSolverBase::MATRIX_VIEW const & localMatrix,
                                                                   arrayView1d< real64 > const & localRhs ) override;

  virtual integer numFluidComponents() const override { return this->flowSolver()->numFluidComponents(); }

  virtual string getFlowDofKey() const override { return CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(); }

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_ */
