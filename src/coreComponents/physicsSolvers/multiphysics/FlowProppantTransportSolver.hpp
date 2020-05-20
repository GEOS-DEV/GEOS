/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowProppantTransportSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_FLOWPROPPANTTRANSPORTSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_FLOWPROPPANTTRANSPORTSOLVER_HPP_

#include "../SolverBase.hpp"

namespace geosx
{

class FlowSolverBase;
class ProppantTransport;

class FlowProppantTransportSolver : public SolverBase
{
public:
  FlowProppantTransportSolver( const std::string & name,
                               Group * const parent );
  ~FlowProppantTransportSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "FlowProppantTransport"; }

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override final;

  virtual void
  ImplicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition * const domain ) override final;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition * const domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto proppantSolverNameString = "proppantSolverName";
    constexpr static auto flowSolverNameString = "flowSolverName";

  } flowProppantTransportSolverViewKeys;


  void PreStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      integer const cycleNumber,
                      DomainPartition * domain );

  void PostStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition * domain );

protected:
  virtual void PostProcessInput() override final;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;


private:

  string m_proppantSolverName;
  string m_flowSolverName;

  FlowSolverBase * m_flowSolver;
  ProppantTransport * m_proppantSolver;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_FLOWPROPPANTTRANSPORTSOLVER_HPP_ */
