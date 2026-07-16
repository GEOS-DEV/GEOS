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
 * @file PhaseFieldFractureSolver.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"

namespace geos
{

class PhaseFieldFractureSolver : public CoupledSolver< SolidMechanicsLagrangianFEM, PhaseFieldDamageFEM >
{
public:

  using Base = CoupledSolver< SolidMechanicsLagrangianFEM, PhaseFieldDamageFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /**
   * @brief main constructor for PhaseFieldFractureSolver objects
   * @param name the name of this instantiation of PhaseFieldFractureSolver in the repository
   * @param parent the parent group of this instantiation of PhaseFieldFractureSolver
   */
  PhaseFieldFractureSolver( const string & name,
                            Group * const parent );

  /// Destructor for the class
  ~PhaseFieldFractureSolver() override = default;

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new PhaseFieldFractureSolver object through the object catalog.
   */
  static string catalogName()
  {
    return "PhaseFieldFracture";
  }
  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "PhaseFieldFracture"; }

  enum class SolverType : integer
  {
    SolidMechanics = 0,
    Damage = 1
  };

  virtual void postInputInitialization() override final;

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the damage solver
   * @return a pointer to the damage solver
   */
  PhaseFieldDamageFEM * damageSolver() const
  {
    return std::get< toUnderlying( SolverType::Damage ) >( m_solvers );
  }

  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType ) override final;

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final {}

private:

  /**
   * @brief interpolate the nodal damage field onto the quadrature points of the solid constitutive model
   * @param[in] domain the domain partition
   */
  void mapDamageToQuadrature( DomainPartition & domain );

};


template< typename FE_TYPE >
struct DamageInterpolationKernel
{
  DamageInterpolationKernel( CellElementSubRegion const & subRegion ):
    m_numElems( subRegion.size() )
  {}

  void interpolateDamage( arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes,
                          arrayView1d< real64 const > const nodalDamage,
                          arrayView2d< real64 > damageFieldOnMaterial )
  {
    forAll< parallelDevicePolicy<> >( m_numElems, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
      constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

      for( localIndex q = 0; q < n_q_points; ++q )
      {
        real64 N[ numNodesPerElement ];
        FE_TYPE::calcN( q, N );

        damageFieldOnMaterial( k, q ) = 0;
        for( localIndex a = 0; a < numNodesPerElement; ++a )
        {
          damageFieldOnMaterial( k, q ) += N[a] * nodalDamage[elemToNodes( k, a )];
        }
      }

    } );
  }

  localIndex m_numElems;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_ */
