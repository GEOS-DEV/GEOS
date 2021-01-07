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

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_ELASTICSTATICFEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_ELASTICSTATICFEM_HPP_

#include "common/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

namespace geosx
{

class ElasticStaticFEM : public SolverBase
{
public:
  ElasticStaticFEM() = delete;

  ElasticStaticFEM( std::string const & name,
                    Group * const parent );

  virtual ~ElasticStaticFEM() override;

  static string CatalogName() { return "ElasticStaticFEM"; }
  
  virtual void 
  RegisterDataOnMesh( Group * const MeshBodies ) override final;

  virtual real64 
  SolverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual void 
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void 
  ImplicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override;


  virtual void 
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;


  virtual void 
  SetupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = false ) override;

  virtual void
  AssembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void 
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  void ApplyDisplacementBC_implicit( real64 const time,
                                     DofManager const & dofManager,
                                     DomainPartition & domain,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs );

  void CRSApplyTractionBC( real64 const time,
                           DofManager const & dofManager,
                           DomainPartition & domain,
                           arrayView1d< real64 > const & localRhs );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr auto solidMaterialNamesString = "solidMaterialNames";
  };

protected:

  array1d< string > m_solidMaterialNames;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_ELASTICSTATICFEM_HPP_ */
