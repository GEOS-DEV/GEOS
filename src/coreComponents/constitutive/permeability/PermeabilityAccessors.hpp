/*	
 * ------------------------------------------------------------------------------------------------------------	
 * SPDX-License-Identifier: LGPL-2.1-only	
 *	
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC	
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University	
 * Copyright (c) 2018-2020 Total, S.A	
 * Copyright (c) 2020-     GEOSX Contributors	
 * All right reserved	
 *	
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.	
 * ------------------------------------------------------------------------------------------------------------	
 */

#ifndef GEOSX_PROPPANTTRANSPORTKERNELS_PERMEABILITYACCESSORS_HPP
#define GEOSX_PROPPANTTRANSPORTKERNELS_PERMEABILITYACCESSORS_HPP

#include "PermeabilityExtrinsicData.hpp"

#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

class PermeabilityAccessorsImpl
{
public:
  PermeabilityAccessorsImpl( ElementRegionManager const & elemManager,
                             string const & solverName )
    : m_impl( elemManager, solverName )
  { }

  auto permeability()
  {
    return m_impl.get< extrinsicMeshData::permeability::permeability >();
  }

  auto permeabilityMultiplier()
  {
    return m_impl.get< extrinsicMeshData::permeability::permeabilityMultiplier >();
  }

private:
  StencilMaterialAccessors< constitutive::PermeabilityBase,
                            extrinsicMeshData::permeability::permeability,
                            extrinsicMeshData::permeability::permeabilityMultiplier > m_impl;
};

} // end of namespace ProppantTransportKernels

namespace CompositionalMultiphaseFVMKernels
{

class PermeabilityAccessorsImpl
{
public:
  PermeabilityAccessorsImpl( ElementRegionManager const & elemManager,
                             string const & solverName )
    : m_impl( elemManager, solverName )
  { }

  auto permeability()
  {
    return m_impl.get< extrinsicMeshData::permeability::permeability >();
  }

  auto permeability() const
  {
    return m_impl.get< extrinsicMeshData::permeability::permeability >();
  }

  auto dPerm_dPressure()
  {
    return m_impl.get< extrinsicMeshData::permeability::dPerm_dPressure >();
  }

  auto dPerm_dPressure() const
  {
    return m_impl.get< extrinsicMeshData::permeability::dPerm_dPressure >();
  }

private:
  StencilMaterialAccessors< constitutive::PermeabilityBase,
                            extrinsicMeshData::permeability::permeability,
                            extrinsicMeshData::permeability::dPerm_dPressure > m_impl;
};

} // end of namespace CompositionalMultiphaseFVMKernels

namespace SinglePhaseFVMKernels
{

using PermeabilityAccessorsImpl = CompositionalMultiphaseFVMKernels::PermeabilityAccessorsImpl;

} // end of namespace SinglePhaseFVMKernels

} // end of namespace geosx

#endif //GEOSX_PROPPANTTRANSPORTKERNELS_PERMEABILITYACCESSORS_HPP
