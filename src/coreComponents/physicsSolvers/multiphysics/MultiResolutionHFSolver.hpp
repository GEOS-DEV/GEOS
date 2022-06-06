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
 * @file MultiResolutionHFSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIRESOLUTIONHFSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIRESOLUTIONHFSOLVER_HPP_

#include "physicsSolvers/multiphysics/PhaseFieldFractureSolver.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"

namespace geosx
{

class SurfaceGenerator;


class MultiResolutionHFSolver : public SolverBase
{
public:
  MultiResolutionHFSolver( const string & name,
                           Group * const parent );

  ~MultiResolutionHFSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "MultiResolutionHF";
  }

  virtual void RegisterDataOnMesh( Group & MeshBodies );

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  //virtual void updateNodeMaps( MeshLevel const & base, MeshLevel const & patch );

  virtual void setInitialCrackDamageBCs( MeshLevel const & patch, MeshLevel const & base );

  virtual void prepareSubProblemBCs( MeshLevel const & base,
                                     MeshLevel & patch );

  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  virtual void setNextDt( real64 const & currentDt,
                          real64 & nextDt ) override;


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * baseSolverNameString() { return "baseSolverName"; }

    constexpr static char const * patchSolverNameString() { return "patchSolverName"; }

    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * surfaceGeneratorNameString() { return "surfaceGeneratorName"; }

    constexpr static char const * maxNumResolvesString() { return "maxNumResolves"; }

  };

protected:
  virtual void postProcessInput() override final;

  virtual void
  initializePostInitialConditionsPreSubGroups() override final;

private:

  // name of the base efem solver
  string m_baseSolverName;

  //name of patch phase-field fracture solver
  string m_patchSolverName;

  // list of damage dofs to be fixed in the subdomain boundary
  array1d< localIndex > m_nodeFixDamage;

  // list of disp dofs to be fixed in the subdomain boundary
  array1d< localIndex > m_nodeFixDisp;

  // list of displacement values to be prescribed in the subdomain boundary
  array2d< real64 > m_fixedDispList;

  // pointer to base efem solver
  SolidMechanicsEmbeddedFractures * m_baseSolver;

  // pointer to patch phase-field fracture solver
  PhaseFieldFractureSolver * m_patchSolver;

  integer m_maxNumResolves;
  integer m_numResolves[2];

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIRESOLUTIONHFSOLVER_HPP_ */
