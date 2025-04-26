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
 * @file Electrostatics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_ELECTROSTATICS_HPP_
#define GEOS_PHYSICSSOLVERS_ELECTROSTATICS_HPP_


#include "common/format/EnumStrings.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "physicsSolvers/PhysicsSolverBase.hpp"

namespace geos
{

/**
 * @class Electrostatics
 *
 * This class implements a finite element solution to the charge balance equation with interface reaction.
 */
class Electrostatics : public PhysicsSolverBase
{
public:
  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  Electrostatics(const string& name, Group* const parent);

  /// Destructor
  virtual ~Electrostatics() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "Electrostatics"; }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  // virtual void initializePreSubGroups() override;

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh(Group & meshBodies) override final;

  virtual real64 solverStep(real64 const& time_n, real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition& domain) override;

  virtual void implicitStepSetup(real64 const& time_n, real64 const& dt,
                                 DomainPartition& domain) override;

  virtual void setupDofs(DomainPartition const& domain, DofManager& dofManager) const override;

  virtual void setupSystem(DomainPartition& domain, DofManager& dofManager,
                           CRSMatrix<real64, globalIndex>& localMatrix,
                           ParallelVector& rhs, ParallelVector& solution,
                           bool const setSparsity = false) override;

  virtual void assembleSystem(real64 const time, real64 const dt,
                              DomainPartition& domain, DofManager const& dofManager,
                              CRSMatrixView<real64, globalIndex const> const& localMatrix,
                              arrayView1d<real64> const& localRhs) override;

  virtual void applyBoundaryConditions(real64 const time, real64 const dt,
                                       DomainPartition& domain, DofManager const& dofManager,
                                       CRSMatrixView<real64, globalIndex const> const& localMatrix,
                                       arrayView1d<real64> const& localRhs) override;

  // virtual real64 calculateResidualNorm(real64 const& time_n, real64 const& dt,
  //                                      DomainPartition const& domain, DofManager const& dofManager,
  //                                      arrayView1d<real64 const> const& localRhs) override;

  virtual void applySystemSolution(DofManager const& dofManager,
                                   arrayView1d<real64 const> const& localSolution,
                                   real64 const scalingFactor, real64 const dt,
                                   DomainPartition& domain) override;

  virtual void updateState(DomainPartition& domain) override final;

  virtual void resetStateToBeginningOfStep(DomainPartition& GEOS_UNUSED_PARAM(domain)) override;

  virtual void implicitStepComplete(real64 const& time, real64 const& dt,
                                    DomainPartition& domain) override;

  void applyPotentialBC(real64 const time, DofManager const& dofManager, DomainPartition& domain,
                        CRSMatrixView<real64, globalIndex const> const& localMatrix,
                        arrayView1d<real64> const& localRhs);

  void applyCurrentBC(real64 const time, DofManager const& dofManager,
                      DomainPartition& domain, arrayView1d<real64> const& localRhs);

  enum class TimeIntegrationOption : integer
  {
    QuasiStatic,
    ImplicitTransient
  };

  struct viewKeyStruct : public PhysicsSolverBase::viewKeyStruct
  {
    static constexpr char const* timeIntegrationOption() { return "timeIntegrationOption"; }
    static constexpr char const* fieldVarName() { return "fieldName"; }
    static constexpr char const* electroMaterialNamesString() {return "electroMaterialNames";}
  };

private:
  string m_fieldName;
  TimeIntegrationOption m_timeIntegrationOption;
};

ENUM_STRINGS(Electrostatics::TimeIntegrationOption,
             "QuasiStatic", "ImplicitTransient");
} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_ELECTROSTATICS_HPP_ */