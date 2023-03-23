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
 * @file PhaseFieldFractureSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"

namespace geosx
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

  PhaseFieldFractureSolver( const string & name,
                            Group * const parent );

  ~PhaseFieldFractureSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "PhaseFieldFracture";
  }

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "PhaseFieldFracture"; }

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  void imposeFakeBackgroundPressures( DomainPartition & domain );

  void mapDamageAndGradientToQuadrature( DomainPartition & domain );

  enum class SolverType : integer
  {
    SolidMechanics = 0,
    Damage = 1
  };

  virtual void postProcessInput() override final;

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  PhaseFieldDamageFEM * damageSolver() const
  {
    return std::get< toUnderlying( SolverType::Damage ) >( m_solvers );
  }

  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, integer const idx ) override final;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * subcyclingOptionString() { return "subcycling"; }
    constexpr static char const * pressureEffectsString() { return "pressureEffects"; }
  };

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final {}

private:

  integer m_pressureEffects; //only use case is the MultiResolution solver

};

template< typename FE_TYPE >
struct DamageInterpolationKernel
{
  DamageInterpolationKernel( CellElementSubRegion const & subRegion ):
    m_numElems( subRegion.size() )
  {}

  void interpolateDamage( arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes,
                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodes,
                          arrayView1d< real64 const > const nodalDamage,
                          arrayView2d< real64 > damageFieldOnMaterial,
                          arrayView3d< real64 > damageGradOnMaterial )
  {
    forAll< serialPolicy >( m_numElems, [=] ( localIndex const k )
    {
      constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
      constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

      real64 xLocal[ numNodesPerElement ][ 3 ];
      real64 nodalDamageLocal[ numNodesPerElement ];

      for( localIndex a = 0; a < numNodesPerElement; ++a )
      {
        localIndex const localNodeIndex = elemToNodes( k, a );

        for( int dim=0; dim < 3; ++dim )
        {
          xLocal[a][dim] = xNodes[ localNodeIndex ][dim];
        }

        nodalDamageLocal[ a ] = nodalDamage[ localNodeIndex ];
      }

      for( localIndex q = 0; q < n_q_points; ++q )
      {
        real64 N[ numNodesPerElement ];
        FE_TYPE::calcN( q, N );

        real64 dNdX[ numNodesPerElement ][ 3 ];

        real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );

        GEOSX_UNUSED_VAR( detJ );

        real64 qDamage = 0.0;
        real64 qDamageGrad[3] = {0, 0, 0};
        FE_TYPE::valueAndGradient( N, dNdX, nodalDamageLocal, qDamage, qDamageGrad );

        damageFieldOnMaterial( k, q ) = qDamage;

        for( int dim=0; dim < 3; ++dim )
        {
          damageFieldOnMaterial( k, q ) += N[a] * nodalDamage[elemToNodes( k, a )];
          damageGradOnMaterial[k][q][dim] = qDamageGrad[dim];
        }
      }
    } );
  }

  localIndex m_numElems;
};

// template< typename FE_TYPE >
// struct DamageInterpolationKernel
// {
//   DamageInterpolationKernel( CellElementSubRegion const & subRegion ):
//     m_numElems( subRegion.size() )
//   {}

//   void interpolateDamage( arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes,
//                           arrayView1d< real64 const > const nodalDamage,
//                           arrayView2d< real64 > damageFieldOnMaterial )
//   {
//     forAll< parallelDevicePolicy<> >( m_numElems, [=] GEOSX_HOST_DEVICE ( localIndex const k )
//     {
//       constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
//       constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

//       for( localIndex q = 0; q < n_q_points; ++q )
//       {
//         real64 N[ numNodesPerElement ];
//         FE_TYPE::calcN( q, N );

//         damageFieldOnMaterial( k, q ) = 0;
//         for( localIndex a = 0; a < numNodesPerElement; ++a )
//         {
//           damageFieldOnMaterial( k, q ) += N[a] * nodalDamage[elemToNodes( k, a )];
//           //solution is probably not going to work because the solution of the coupled solver
//           //has both damage and displacements. Using the damageResult field from the Damage solver
//           //is probably better
//         }
//       }

//     } );
//   }

//   localIndex m_numElems;
// };

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_PHASEFIELDFRACTURESOLVER_HPP_ */
