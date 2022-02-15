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

#ifndef GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXSLURRYFLUIDACCESSORS_HPP
#define GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXSLURRYFLUIDACCESSORS_HPP

#include "SlurryFluidExtrinsicData.hpp"

#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geosx
{

namespace ProppantTransportKernels
{

class CellBasedFluxSlurryFluidAccessorsImpl
{
public:
  CellBasedFluxSlurryFluidAccessorsImpl( ElementRegionManager const & elemManager,
                                         string const & solverName )
    : m_impl( elemManager, solverName )
  { }

  auto density()
  {
    return m_impl.get< extrinsicMeshData::singlefluid::density >();
  }

  auto viscosity()
  {
    return m_impl.get< extrinsicMeshData::singlefluid::viscosity >();
  }

private:
  StencilMaterialAccessors< constitutive::SlurryFluidBase,
                            extrinsicMeshData::singlefluid::density,
                            extrinsicMeshData::singlefluid::viscosity > m_impl;
};

} // end of namespace ProppantTransportKernels

} // end of namespace geosx

#endif //GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXSLURRYFLUIDACCESSORS_HPP
