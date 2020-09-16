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
 * @file HydrofractureSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

#define Trilinos 1
#define Hypre 2
#define Petsc 3
#undef USING_TRILINOS
#if ( GEOSX_LA_INTERFACE == Trilinos )
#define USING_TRILINOS
#endif

namespace geosx
{

class FlowSolverBase;
class SolidMechanicsLagrangianFEM;

class HydrofractureSolver : public SolverBase
{
public:
  HydrofractureSolver( const std::string & name,
                       Group * const parent );

  ~HydrofractureSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "Hydrofracture";
  }

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;
#endif

  virtual void SetupDofs( DomainPartition const * const domain,
                          DofManager & dofManager ) const override;

  virtual void SetupSystem( DomainPartition * const domain,
                            DofManager & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition * const domain,
                     DofManager & dofManager,
                     ParallelMatrix & matrix,
                     ParallelVector & rhs,
                     ParallelVector & solution ) override final;

  virtual void ImplicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition * const domain ) override final;

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition * const domain,
                               DofManager const & dofManager,
                               ParallelMatrix & matrix,
                               ParallelVector & rhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition * const domain,
                                        DofManager const & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void SolveSystem( DofManager const & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual real64
  ScalingForSystemSolution( DomainPartition const * const domain,
                            DofManager const & dofManager,
                            ParallelVector const & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * const domain ) override;

  virtual void SetNextDt( real64 const & currentDt,
                          real64 & nextDt ) override;


  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition * const domain ) override;

  void UpdateDeformationForCoupling( DomainPartition * const domain );

//  void ApplyFractureFluidCoupling( DomainPartition * const domain,
//                                   systemSolverInterface::EpetraBlockSystem & blockSystem );

  void AssembleForceResidualDerivativeWrtPressure( DomainPartition * const domain,
                                                   ParallelMatrix * const matrix01,
                                                   ParallelVector * const rhs0 );

  void AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const * const domain,
                                                           ParallelMatrix * const matrix10,
                                                           ParallelVector * const rhs0,
							   real64 const dt);


  real64 SplitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition * const domain );

  //TJ: getter for the tip location after each time step
  inline real64 getConvergedTipLoc() const {return m_convergedTipLoc;}

  //TJ: getter for the total time
  inline real64 getTotalTime() const {return m_totalTime;}

  void initializeNewFaceElements( DomainPartition const & domain );

  enum class couplingTypeOption : int
  {
    FIM,
    SIM_FixedStress
  };

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto couplingTypeOptionString = "couplingTypeOptionEnum";
    constexpr static auto couplingTypeOptionStringString = "couplingTypeOption";

    constexpr static auto totalMeanStressString = "totalMeanStress";
    constexpr static auto oldTotalMeanStressString = "oldTotalMeanStress";

    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto fluidSolverNameString = "fluidSolverName";

    constexpr static auto contactRelationNameString = "contactRelationName";
    constexpr static auto maxNumResolvesString = "maxNumResolves";

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    constexpr static auto separationCoeff0String = "separationCoeff0";
    constexpr static auto apertureAtFailureString = "apertureAtFailure";
#endif

  } HydrofractureSolverViewKeys;

protected:
  virtual void PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;

private:

  string m_solidSolverName;
  string m_flowSolverName;
  string m_contactRelationName;
  string m_couplingTypeOptionString;

  couplingTypeOption m_couplingTypeOption;

  SolidMechanicsLagrangianFEM * m_solidSolver;
  FlowSolverBase * m_flowSolver;

  //TJ: a flag to indicate of finding the tip location
  //    through fixed-pointer iterations
  //TJ: Since we have globalSolverIterFlag, we don't need this private member anymore
//  bool m_tipIterationFlag;

  //TJ: old tip location
  real64 m_oldTipLocation;

  //TJ: new tip location
  real64 m_newTipLocation;

  //TJ: tip location after each time step
  real64 m_convergedTipLoc;

  //TJ: an array to store the tip locations during the iteration
  array1d<real64> m_tipLocationHistory;

  //TJ: tip element, the element that is partially opened
  localIndex m_tipElement;

  //TJ: element on the boundary of the channel region. The element
  //    itself is fully opened, and its aperture is used to locate
  //    the tip
  localIndex m_channelElement;

  //TJ: total_time, needed for evaluating velocity
  real64 m_totalTime;

#ifdef USING_TRILINOS
  real64 m_densityScaling;
  real64 m_pressureScaling;
#endif

  std::unique_ptr< ParallelMatrix > m_blockDiagUU;

  ParallelMatrix m_matrix01;
  ParallelMatrix m_matrix10;

  ParallelMatrix m_permutationMatrix0; // it's used to have the output based on global ordering
  ParallelMatrix m_permutationMatrix1; // it's used to have the output based on global ordering

  integer m_maxNumResolves;
  integer m_numResolves[2];
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_ */
