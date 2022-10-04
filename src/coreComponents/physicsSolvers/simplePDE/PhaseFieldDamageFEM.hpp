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

/**
 * @brief class that implements the solver for the damage equation in a phase-field solver of fracture
 *
 * This class implements the an FEM solver for the reaction-diffusion equation that governs damage
 * evolution in a phase-field description of brittle fracture. The drivinf force for damage growth
 * is the strain energy density, that is stored in the solid constitutive option and retrived in the
 * assembly routine of this solver.
 */

class PhaseFieldDamageFEM : public SolverBase
{
public:
  PhaseFieldDamageFEM( const string & name, Group * const parent );

  virtual ~PhaseFieldDamageFEM() override;

  static string catalogName()
  {
    return "PhaseFieldDamageFEM";
  }

  virtual void registerDataOnMesh( Group & meshBodies ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions implement all the standard solution steps of an FEM solver
   * class
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

  void applyDirichletBCImplicit( real64 const time,
                                 DofManager const & dofManager,
                                 DomainPartition & domain,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs );

  /**@}*/

  ///enum with possible time integration options, but Explicit is not implemented
  enum class timeIntegrationOption
  {
    SteadyState,
    ImplicitTransient,
    ExplicitTransient
  };

  ///struct with names used in the xml file
  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * coeffNameString() { return "coeffField"; }
    static constexpr char const * localDissipationOptionString() { return "localDissipation"; }
    static constexpr char const * solidModelNamesString() { return "solidMaterialNames"; }

    dataRepository::ViewKey timeIntegrationOption = { "timeIntegrationOption" };
    dataRepository::ViewKey fieldVarName = { "fieldName" };
  } PhaseFieldDamageFEMViewKeys;

  inline ParallelVector const * getSolution() const
  {
    return &m_solution;
  }

  inline globalIndex getSize() const
  {
    return m_matrix.numGlobalRows();
  }

  string const & getFieldName() const
  {
    return m_fieldName;
  }

protected:
  virtual void postProcessInput() override final;

private:
  ///name of the primary variable, will be switched to Damage in future PR
  string m_fieldName;
  ///the stable time step
  stabledt m_stabledt;
  ///the time integration option, can probably be set to SteadyState in all cases
  timeIntegrationOption m_timeIntegrationOption;
  ///the type of local dissipation function, can either be Linear (AT1 model) or Quadratic (AT2 model)
  string m_localDissipationOption;
  ///this was used for debugging in an older version of the code, it is already set to be removed in another PR
  array1d< real64 > m_coeff;

  PhaseFieldDamageFEM();
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGE_HPP_ */
