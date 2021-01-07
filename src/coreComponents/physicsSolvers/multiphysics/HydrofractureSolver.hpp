/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

#include "common/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include <unordered_set>

namespace geosx
{

class FlowSolverBase;
class SolidMechanicsLagrangianFEM;
//TJ: add the surfaceGenerator solver;
class SurfaceGenerator;

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

  virtual void SetupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void SetupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            array1d< real64 > & localRhs,
                            array1d< real64 > & localSolution,
                            bool const setSparsity = true ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void ImplicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override final;

  void EikonalEquationSolver( DomainPartition & domain,
			      FaceElementSubRegion & subRegion,
			      std::unordered_set<localIndex> const & partiallyOpenFaceElmts );

  real64 CalculateSignedDistance1stOrder( localIndex const node,
                                          std::array<std::unordered_set<localIndex>, 2> const & neighbors,
                                          std::unordered_map<localIndex, int> & nodeStatus,
                                          NodeManager * const nodeManager);

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual real64
  ScalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void SetNextDt( real64 const & currentDt,
                          real64 & nextDt ) override;


  virtual real64 ExplicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  void UpdateDeformationForCoupling( DomainPartition & domain );

//  void ApplyFractureFluidCoupling( DomainPartition * const domain,
//                                   systemSolverInterface::EpetraBlockSystem & blockSystem );

  void AssembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix );


  real64 SplitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  void initializeNewFaceElements( DomainPartition const & domain );



  enum class CouplingTypeOption : integer
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
    constexpr static auto surfaceGeneratorSolverNameString = "surfaceGeneratorName";

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

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void addFluxApertureCouplingNNZ( DomainPartition & domain,
                                   DofManager & dofManager,
                                   arrayView1d< localIndex > const & rowLengths ) const;


  /**
   * @Brief add the sparsity pattern induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                               DofManager & dofManager,
                                               SparsityPatternView< globalIndex > const & pattern ) const;

private:

  string m_solidSolverName;
  string m_flowSolverName;
  string m_contactRelationName;
  //TJ: name of the surfaceGenerator solver
  string m_surfaceGeneratorName;

  CouplingTypeOption m_couplingTypeOption;

  SolidMechanicsLagrangianFEM * m_solidSolver;
  FlowSolverBase * m_flowSolver;
  //TJ: surfaceGenerator solver to access tip element related
  //    quantities, such as m_tipNodes and m_trailingFaces
  SurfaceGenerator * m_surfaceGeneratorSolver;

  std::unique_ptr< ParallelMatrix > m_blockDiagUU;

  integer m_maxNumResolves;
  integer m_numResolves[2];
};

ENUM_STRINGS( HydrofractureSolver::CouplingTypeOption, "FIM", "SIM_FixedStress" )

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_ */
