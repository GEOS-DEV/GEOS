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
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_

#include "SolidMechanicsLagrangianFEM.hpp"
#include "SolidMechanicsLagrangianSSLEKernels.hpp"

namespace geosx
{

/**
 * @class SolidMechanicsLagrangianSSLE
 *
 * This class contains an implementation of a small strain linear elastic solution to the equations of motion which are
 * called through the interface in SolidMechanicsLagrangianFEM.
 */
class SolidMechanicsLagrangianSSLE : public SolidMechanicsLagrangianFEM
{
public:
  SolidMechanicsLagrangianSSLE( string const & name,
                                Group * const parent );
  virtual ~SolidMechanicsLagrangianSSLE() override;

  static string CatalogName() { return "SolidMechanicsLagrangianSSLE"; }

  virtual void
  updateStress( DomainPartition * const domain ) override;


  virtual void ApplySystemSolution( DofManager const & dofManager,
                                    ParallelVector const & solution,
                                    real64 const scalingFactor,
                                    DomainPartition * const domain  ) override;

  virtual real64
  ExplicitElementKernelLaunch( localIndex NUM_NODES_PER_ELEM,
                               localIndex NUM_QUADRATURE_POINTS,
                               constitutive::ConstitutiveBase * const constitutiveRelation,
                               set<localIndex> const & elementList,
                               arrayView2d<localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM> const & elemsToNodes,
                               arrayView3d< R1Tensor const> const & dNdX,
                               arrayView2d<real64 const> const & detJ,
                               arrayView1d<R1Tensor const> const & u,
                               arrayView1d<R1Tensor const> const & vel,
                               arrayView1d<R1Tensor> const & acc,
                               arrayView2d<R2SymTensor> const & stress,
                               real64 const dt ) const override
  {
    using ExplicitKernel = SolidMechanicsLagrangianSSLEKernels::ExplicitKernel;
    return SolidMechanicsLagrangianFEMKernels::
           ElementKernelLaunchSelector<ExplicitKernel>( NUM_NODES_PER_ELEM,
                                                        NUM_QUADRATURE_POINTS,
                                                        constitutiveRelation,
                                                        elementList,
                                                        elemsToNodes,
                                                        dNdX,
                                                        detJ,
                                                        u,
                                                        vel,
                                                        acc,
                                                        stress,
                                                        dt );
  }

  virtual real64
  ImplicitElementKernelLaunch( localIndex NUM_NODES_PER_ELEM,
                               localIndex NUM_QUADRATURE_POINTS,
                               constitutive::ConstitutiveBase * const constitutiveRelation,
                               localIndex const numElems,
                               real64 const dt,
                               arrayView3d<R1Tensor const> const & dNdX,
                               arrayView2d<real64 const > const& detJ,
                               FiniteElementBase const * const fe,
                               arrayView1d< integer const > const & elemGhostRank,
                               arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes,
                               arrayView1d< globalIndex const > const & globalDofNumber,
                               arrayView1d< R1Tensor const > const & disp,
                               arrayView1d< R1Tensor const > const & uhat,
                               arrayView1d< R1Tensor const > const & vtilde,
                               arrayView1d< R1Tensor const > const & uhattilde,
                               arrayView2d< real64 const > const & density,
                               arrayView1d< real64 const > const & fluidPressure,
                               arrayView1d< real64 const > const & deltaFluidPressure,
                               arrayView1d< real64 const > const & biotCoefficient,
                               timeIntegrationOption const tiOption,
                               real64 const stiffnessDamping,
                               real64 const massDamping,
                               real64 const newmarkBeta,
                               real64 const newmarkGamma,
                               R1Tensor const & gravityVector,
                               DofManager const * const dofManager,
                               ParallelMatrix * const matrix,
                               ParallelVector * const rhs ) const override
  {
    GEOSX_MARK_FUNCTION;
    using ImplicitKernel = SolidMechanicsLagrangianSSLEKernels::ImplicitKernel;
    return SolidMechanicsLagrangianFEMKernels::
           ElementKernelLaunchSelector<ImplicitKernel>( NUM_NODES_PER_ELEM,
                                                        NUM_QUADRATURE_POINTS,
                                                        constitutiveRelation,
                                                        numElems,
                                                        dt,
                                                        dNdX,
                                                        detJ,
                                                        fe,
                                                        elemGhostRank,
                                                        elemsToNodes,
                                                        globalDofNumber,
                                                        disp,
                                                        uhat,
                                                        vtilde,
                                                        uhattilde,
                                                        density,
                                                        fluidPressure,
                                                        deltaFluidPressure,
                                                        biotCoefficient,
                                                        tiOption,
                                                        stiffnessDamping,
                                                        massDamping,
                                                        newmarkBeta,
                                                        newmarkGamma,
                                                        gravityVector,
                                                        dofManager,
                                                        matrix,
                                                        rhs );
  }
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_ */
