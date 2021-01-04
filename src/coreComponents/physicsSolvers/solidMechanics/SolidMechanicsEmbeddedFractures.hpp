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

/*
 * SolidMechanicsEmbeddedFractures.hpp
 *
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
using namespace constitutive;

class SolidMechanicsLagrangianFEM;

class SolidMechanicsEmbeddedFractures : public SolverBase
{
public:
  SolidMechanicsEmbeddedFractures( const std::string & name,
                                   Group * const parent );

  ~SolidMechanicsEmbeddedFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  {
    return "SolidMechanicsEmbeddedFractures";
  }

  virtual void RegisterDataOnMesh( dataRepository::Group * const MeshBodies ) override final;

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

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void ResetStateToBeginningOfStep( DomainPartition & domain ) override final;

  void updateState( DomainPartition & domain );

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;


  void AddCouplingNumNonzeros( DomainPartition & domain,
                               DofManager & dofManager,
                               arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void AddCouplingSparsityPattern( DomainPartition const & domain,
                                   DofManager const & dofManager,
                                   SparsityPatternView< globalIndex > const & pattern ) const;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto solidSolverNameString = "solidSolverName";

    constexpr static auto contactRelationNameString = "contactRelationName";

    constexpr static auto dispJumpString = "displacementJump";

    constexpr static auto deltaDispJumpString = "deltaDisplacementJump";

    constexpr static auto fractureRegionNameString = "fractureRegionName";

    constexpr static auto fractureTractionString = "fractureTraction";

    constexpr static auto dTraction_dJumpString = "dTraction_dJump";

  } SolidMechanicsEmbeddedFracturesViewKeys;

  string & getContactRelationName() { return m_contactRelationName; };

  string const & getContactRelationName() const { return m_contactRelationName; };

  string & getFractureRegionName() { return m_fractureRegionName; };

  string const & getFractureRegionName() const { return m_fractureRegionName; };

  void applyTractionBC( real64 const time_n,
                        real64 const dt,
                        DomainPartition & domain );

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( Group * const problemManager ) override final;

  virtual void PostProcessInput() override final;



private:

  /// Solid mechanics solver name
  string m_solidSolverName;

  /// fracture region name
  string m_fractureRegionName;

  /// pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// contact relation name string
  string m_contactRelationName;
};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_ */
