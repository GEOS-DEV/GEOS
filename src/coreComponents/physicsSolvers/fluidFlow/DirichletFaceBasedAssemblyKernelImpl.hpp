/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DirichletFaceBasedAssemblyKernelImpl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_DIRICHLETFACEBASEDASSEMBLYKERNELIMPL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_DIRICHLETFACEBASEDASSEMBLYKERNELIMPL_HPP

#include "DirichletFaceBasedAssemblyKernel.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

template< integer NUM_COMP, integer NUM_DOF, typename FLUIDWRAPPER >
DirichletFaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, FLUIDWRAPPER >::
DirichletFaceBasedAssemblyKernel( integer const numPhases,
                                  globalIndex const rankOffset,
                                  FaceManager const & faceManager,
                                  BoundaryStencilWrapper const & stencilWrapper,
                                  FLUIDWRAPPER const & fluidWrapper,
                                  DofNumberAccessor const & dofNumberAccessor,
                                  CompFlowAccessors const & compFlowAccessors,
                                  MultiFluidAccessors const & multiFluidAccessors,
                                  CapPressureAccessors const & capPressureAccessors,
                                  PermeabilityAccessors const & permeabilityAccessors,
                                  real64 const dt,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs,
                                  BitFlags< FaceBasedAssemblyKernelFlags > kernelFlags )
  : Base( numPhases,
          rankOffset,
          stencilWrapper,
          dofNumberAccessor,
          compFlowAccessors,
          multiFluidAccessors,
          capPressureAccessors,
          permeabilityAccessors,
          dt,
          localMatrix,
          localRhs,
          kernelFlags ),
  m_facePres( faceManager.getField< fields::flow::facePressure >() ),
  m_faceTemp( faceManager.getField< fields::flow::faceTemperature >() ),
  m_faceCompFrac( faceManager.getField< fields::flow::faceGlobalCompFraction >() ),
  m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
  m_fluidWrapper( fluidWrapper )
{}

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_DIRICHLETFACEBASEDASSEMBLYKERNELIMPL_HPP
