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

namespace geosx
{

class SolidMechanicsLagrangianSSLE : public SolidMechanicsLagrangianFEM
{
public:
  SolidMechanicsLagrangianSSLE( string const & name,
                                ManagedGroup * const parent );
  virtual ~SolidMechanicsLagrangianSSLE() override;

  static string CatalogName() { return "SolidMechanicsLagrangianSSLE"; }

  virtual real64
  ExplicitElementKernelLaunchSelector( localIndex const er,
                                       localIndex const esr,
                                       set<localIndex> const & elementList,
                                       arrayView2d<localIndex> const & elemsToNodes,
                                       arrayView3d< R1Tensor > const & dNdX,
                                       arrayView2d<real64> const & detJ,
                                       arrayView1d<R1Tensor> const & u,
                                       arrayView1d<R1Tensor> const & vel,
                                       arrayView1d<R1Tensor> & acc,
                                       ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase>& constitutiveRelations,
                                       ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                                       ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                                       real64 const dt,
                                       localIndex NUM_NODES_PER_ELEM,
                                       localIndex NUM_QUADRATURE_POINTS ) override
  {
    return ExplicitElementKernelLaunchSelectorT<SolidMechanicsLagrangianSSLE>( er,
                                                                              esr,
                                                                              elementList,
                                                                              elemsToNodes,
                                                                              dNdX,
                                                                              detJ,
                                                                              u,
                                                                              vel,
                                                                              acc,
                                                                              constitutiveRelations,
                                                                              meanStress,
                                                                              devStress,
                                                                              dt,
                                                                              NUM_NODES_PER_ELEM,
                                                                              NUM_QUADRATURE_POINTS );
  }

  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
  static real64
  ExplicitElementKernelLaunch( localIndex const er,
                               localIndex const esr,
                               set<localIndex> const & elementList,
                               arrayView2d<localIndex> const & elemsToNodes,
                               arrayView3d< R1Tensor > const & dNdX,
                               arrayView2d<real64> const & detJ,
                               arrayView1d<R1Tensor> const & u,
                               arrayView1d<R1Tensor> const & vel,
                               arrayView1d<R1Tensor> & acc,
                               ElementRegionManager::ConstitutiveRelationAccessor<constitutive::ConstitutiveBase> constitutiveRelations,
                               ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & meanStress,
                               ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const & devStress,
                               real64 const dt );
};

} /* namespace geosx */

#endif /* CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANSSLE_HPP_ */
