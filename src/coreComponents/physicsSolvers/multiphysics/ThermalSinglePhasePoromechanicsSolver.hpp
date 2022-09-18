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
 * @file ThermalSinglePhasePoromechanicsSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_THERMALSINGLEPHASEPOROMECHANICSSOLVER_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class ThermalSinglePhasePoromechanicsSolver : public SolverBase
{
public:
  ThermalSinglePhasePoromechanicsSolver( const string & name,
                                         Group * const parent );
  ~ThermalSinglePhasePoromechanicsSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "ThermalSinglePhasePoromechanics";
  }

  virtual void registerDataOnMesh( Group & MeshBodies ) override final;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  enum class CouplingTypeOption : integer
  {
    FixedStress,
    TightlyCoupled
  };

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * couplingTypeOptionString() { return "couplingTypeOption"; }

    constexpr static char const * totalMeanStressString() { return "totalMeanStress"; }
    constexpr static char const * oldTotalMeanStressString() { return "oldTotalMeanStress"; }

    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
    constexpr static char const * fluidSolverNameString() { return "fluidSolverName"; }
    constexpr static char const * subcyclingOptionString() { return "subcycling"; }
  };

protected:
  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:

  string m_solidSolverName;
  string m_fluidSolverName;
  CouplingTypeOption m_couplingTypeOption;
  integer m_subcyclingOption;

};

ENUM_STRINGS( ThermalSinglePhasePoromechanicsSolver::CouplingTypeOption,
              "FixedStress",
              "TightlyCoupled" );

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_ThermalSinglePhasePoromechanicsSOLVER_HPP_ */
