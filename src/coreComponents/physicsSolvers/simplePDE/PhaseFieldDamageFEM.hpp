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
 * @file PhaseFieldDamageFEM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/SolverBase.hpp"

struct stabledt
{
  double m_maxdt;
};

namespace geosx
{
namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

class PhaseFieldDamageFEM : public SolverBase
{
public:
  PhaseFieldDamageFEM( const string & name, Group * const parent );

  virtual ~PhaseFieldDamageFEM() override;

  static string catalogName()
  {
    return "PhaseFieldDamageFEM";
  }

  static string coupledSolverAttributePrefix() { return "damage"; }

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived
   * classes
   */
  /**@{*/

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void assembleSystem( real64 const time, real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void applyBoundaryConditions( real64 const time, real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64 calculateResidualNorm( DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs ) override;

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final;

  virtual void
  implicitStepSetup( real64 const &,
                     real64 const &,
                     DomainPartition & ) override {}

  virtual void implicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & ) override {}

  /**@}*/

  void applyDirichletBCImplicit( real64 const time,
                                 DofManager const & dofManager,
                                 DomainPartition & domain,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs );

  void setInitialCrackDamageBCs( DomainPartition & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs );

  void setInitialCrackNodes( arrayView1d< localIndex > const & fracturedNodes );

  void applyIrreversibilityConstraint( DofManager const & dofManager,
                                       DomainPartition & domain,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs );

  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
  #if 0
    static constexpr char const * coeffNameString() { return "coeffField"; }
  #endif
    static constexpr char const * localDissipationOptionString() { return "localDissipation"; }
    static constexpr char const * irreversibilityFlagString() { return "irreversibilityFlag"; }
    static constexpr char const * damageUpperBoundString() { return "damageUpperBound"; }
    static constexpr char const * damageViscosityFlagString() { return "useDamageViscosity"; }
    static constexpr char const * damageViscosityParameterString() { return "damageViscosityCoeff"; }
    static constexpr char const * solidModelNamesString() { return "solidMaterialNames"; }

    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey damageVarName = { "Damage" };
  } PhaseFieldDamageFEMViewKeys;

  inline ParallelVector const * getSolution() const
  {
    return &m_solution;
  }

  inline globalIndex getSize() const
  {
    return m_matrix.numGlobalRows();
  }

  string const getFieldName() const
  {
    return m_damageName;
  }

  void setPressureEffects()
  {
    m_pressureEffectsFlag = 1;
  }

protected:
  virtual void postProcessInput() override final;

private:
  string m_damageName;
  stabledt m_stabledt;
  timeIntegrationOption m_timeIntegrationOption;
  string m_localDissipationOption;
  integer m_pressureEffectsFlag;
  integer m_irreversibilityFlag;
  real64 m_damageUpperBound;
  integer m_damageViscosityFlag;
  real64 m_damageViscosityCoeff;
  array1d< localIndex > m_initialCrack;
  //SortedArrayView< localIndex const > const & m_subdomainElems;

  PhaseFieldDamageFEM();
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_ */
