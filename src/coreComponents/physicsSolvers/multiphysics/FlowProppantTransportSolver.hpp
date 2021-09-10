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
 * @file FlowProppantTransportSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_FLOWPROPPANTTRANSPORTSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_FLOWPROPPANTTRANSPORTSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class ProppantTransport;
class FlowSolverBase;

class FlowProppantTransportSolver : public SolverBase
{
public:
  FlowProppantTransportSolver( const string & name,
                               Group * const parent );
  ~FlowProppantTransportSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "FlowProppantTransport"; }

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * proppantSolverNameString() { return "proppantSolverName"; }
    constexpr static char const * flowSolverNameString() { return "flowSolverName"; }
  };

  void preStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      DomainPartition & domain );

  void postStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       DomainPartition & domain );

protected:

  virtual void postProcessInput() override final;

private:

  string m_proppantSolverName;
  string m_flowSolverName;

  FlowSolverBase * m_flowSolver;
  ProppantTransport * m_proppantSolver;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_FLOWPROPPANTTRANSPORTSOLVER_HPP_ */
