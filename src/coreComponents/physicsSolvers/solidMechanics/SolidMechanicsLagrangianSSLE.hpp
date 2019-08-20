/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#ifndef CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_
#define CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_

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
                                ManagedGroup * const parent );
  virtual ~SolidMechanicsLagrangianSSLE() override;

  static string CatalogName() { return "SolidMechanicsLagrangianSSLE"; }

  virtual real64
  ExplicitElementKernelLaunch( localIndex NUM_NODES_PER_ELEM,
                               localIndex NUM_QUADRATURE_POINTS,
                               constitutive::ConstitutiveBase * const constitutiveRelation,
                               set<localIndex> const & elementList,
                               arrayView2d<localIndex const> const & elemsToNodes,
                               arrayView3d< R1Tensor const> const & dNdX,
                               arrayView2d<real64 const> const & detJ,
                               arrayView1d<R1Tensor const> const & u,
                               arrayView1d<R1Tensor const> const & vel,
                               arrayView1d<R1Tensor> const & acc,
                               arrayView2d<real64> const & meanStress,
                               arrayView2d<R2SymTensor> const & devStress,
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
                                                        meanStress,
                                                        devStress,
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
                               arrayView2d< localIndex const > const & elemsToNodes,
                               arrayView1d< globalIndex const > const & globalDofNumber,
                               arrayView1d< R1Tensor const > const & disp,
                               arrayView1d< R1Tensor const > const & uhat,
                               arrayView1d< R1Tensor const > const & vtilde,
                               arrayView1d< R1Tensor const > const & uhattilde,
                               arrayView1d< real64 const > const & density,
                               arrayView1d< real64 const > const & fluidPressure,
                               arrayView1d< real64 const > const & deltaFluidPressure,
                               arrayView1d< real64 const > const & biotCoefficient,
                               timeIntegrationOption const tiOption,
                               real64 const stiffnessDamping,
                               real64 const massDamping,
                               real64 const newmarkBeta,
                               real64 const newmarkGamma,
                               DofManager const * const dofManager,
                               ParallelMatrix * const matrix,
                               ParallelVector * const rhs ) const override
  {
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
                                                        dofManager,
                                                        matrix,
                                                        rhs );
  }
};

} /* namespace geosx */

#endif /* CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_ */
