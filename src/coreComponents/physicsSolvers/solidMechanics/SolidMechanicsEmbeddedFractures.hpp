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

/*
 * SolidMechanicsEmbeddedFractures.hpp
 *
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

class SolidMechanicsLagrangianFEM;

class SolidMechanicsEmbeddedFractures : public SolverBase
{
public:
  SolidMechanicsEmbeddedFractures( const std::string& name,
                                   Group * const parent );

  ~SolidMechanicsEmbeddedFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "SolidMechanics_FEM_AES";
  }

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

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

  virtual void ImplicitStepComplete( real64 const& time_n,
                                     real64 const& dt,
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

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override final;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition * const domain ) override;


  void AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const * const domain,
                                                           ParallelMatrix * const matrix10,
                                                           ParallelVector * const rhs0 );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto totalMeanStressString = "totalMeanStress";
    constexpr static auto oldTotalMeanStressString = "oldTotalMeanStress";

    constexpr static auto solidSolverNameString = "solidSolverName";

    constexpr static auto contactRelationNameString = "contactRelationName";

    constexpr static auto dispJumpString = "displacementJump";

    constexpr static auto deltaDispJumpString = "deltaDisplacementJump";

  } SolidMechanicsEmbeddedFracturesViewKeys;

protected:
  virtual void PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;

private:

  string m_solidSolverName;

  SolidMechanicsLagrangianFEM * m_solidSolver;

  string m_contactRelationName;

  ParallelMatrix m_matrix11;
  ParallelMatrix m_matrix01;
  ParallelMatrix m_matrix10;

  ParallelVector m_residual1; // Traction balance on the fracture

  ParallelMatrix m_permutationMatrix0; // it's used to have the output based on global ordering
  ParallelMatrix m_permutationMatrix1; // it's used to have the output based on global ordering

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_ */
