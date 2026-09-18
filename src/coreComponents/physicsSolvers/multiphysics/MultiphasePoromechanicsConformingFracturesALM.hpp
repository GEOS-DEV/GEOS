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
 * @file MultiphasePoromechanicsConformingFracturesALM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp"

namespace geos
{

template< typename FLOW_SOLVER = CompositionalMultiphaseBase >
class MultiphasePoromechanicsConformingFracturesALM : public PoromechanicsConformingFractures< MultiphasePoromechanics, FLOW_SOLVER, SolidMechanicsAugmentedLagrangianContact >
{
public:

  using Base = PoromechanicsConformingFractures< MultiphasePoromechanics, FLOW_SOLVER, SolidMechanicsAugmentedLagrangianContact >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_maxFaceNodes;

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFracturesALM"; }

  /**
   * @brief main constructor for MultiphasePoromechanicsConformingFracturesALM objects
   * @param name the name of this instantiation of MultiphasePoromechanicsConformingFracturesALM in the repository
   * @param parent the parent group of this instantiation of MultiphasePoromechanicsConformingFracturesALM
   */
  MultiphasePoromechanicsConformingFracturesALM( const string & name,
                                                 dataRepository::Group * const parent );

  /// Destructor for the class
  ~MultiphasePoromechanicsConformingFracturesALM() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new MultiphasePoromechanicsConformingFracturesALM object through the
   * object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, CompositionalMultiphaseBase > )
    {
      return "MultiphasePoromechanicsConformingFracturesALM";
    }
    else
    {
      return FLOW_SOLVER::catalogName() + "PoromechanicsConformingFracturesALM";
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

  GEOS_MGR_STRATEGY_NOT_SUPPORTED()//TODO: no MGR strategy exists yet for multiphase ALM

  virtual void setSparsityPattern( DomainPartition & domain,
                                   DofManager & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                   SparsityPattern< globalIndex > & pattern ) override final;

  /**@}*/

protected:

  virtual void initializePreSubGroups() override
  {
    Base::initializePreSubGroups();

    GEOS_THROW_IF( this->m_isThermal || this->flowSolver()->isThermal(),
                   GEOS_FMT( "{}: thermal coupling is not supported by {}",
                             this->getName(), this->getCatalogName() ),
                   InputError, this->getDataContext() );
  }

  virtual void assembleForceResidualDerivativeWrtPressure( string const & meshName,
                                                            MeshLevel const & mesh,
                                                            string_array const & regionNames,
                                                            DofManager const & dofManager,
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                            arrayView1d< real64 > const & localRhs ) override final;

  virtual void assembleFluidMassResidualDerivativeWrtDisplacement( string const & meshName,
                                                                    MeshLevel const & mesh,
                                                                    string_array const & regionNames,
                                                                    DofManager const & dofManager,
                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                    arrayView1d< real64 > const & localRhs ) override final;

  virtual string getFlowDofKey() const override { return CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(); }

private:

  struct viewKeyStruct : public Base::viewKeyStruct
  {};

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_ */
