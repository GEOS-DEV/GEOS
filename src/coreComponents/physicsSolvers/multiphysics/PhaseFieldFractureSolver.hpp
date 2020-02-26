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
 * @file PhaseFieldFractureSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_PhaseFieldFractureSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_PhaseFieldFractureSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class PhaseFieldFractureSolver : public SolverBase
{
public:
  PhaseFieldFractureSolver( const std::string& name,
                     Group * const parent );
  ~PhaseFieldFractureSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "PhaseFieldFracture"; }

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


  real64 SplitOperatorStep( real64 const& time_n,
                            real64 const& dt,
                            integer const cycleNumber,
                            DomainPartition * const domain);

  virtual void ApplySystemSolution(  DofManager const &dofManager,
                                     ParallelVector const &solution,
                                     real64 const scalingFactor,
                                     DomainPartition *const domain) override;


  enum class couplingTypeOption : int
  {
    FixedStress,
    TightlyCoupled
  };



  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto couplingTypeOptionString = "couplingTypeOptionEnum";
    constexpr static auto couplingTypeOptionStringString = "couplingTypeOption";

    constexpr static auto totalMeanStressString = "totalMeanStress";
    constexpr static auto oldTotalMeanStressString = "oldTotalMeanStress";

    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto damageSolverNameString = "damageSolverName";
  } PhaseFieldFractureSolverViewKeys;


protected:
  virtual void PostProcessInput() override final;

  virtual void InitializePostInitialConditions_PreSubGroups(dataRepository::Group * const problemManager) override final;


private:

  string m_solidSolverName;
  string m_damageSolverName;
  string m_couplingTypeOptionString;
  couplingTypeOption m_couplingTypeOption;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_PhaseFieldFractureSOLVER_HPP_ */
