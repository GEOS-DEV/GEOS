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
                                         string const & solverName,
                                         arrayView1d< string const > const & regionNames,
                                         arrayView1d< string const > const & materialNames )
    : m_impl( elemManager, solverName, regionNames, materialNames )
  { }

  auto density()
  {
    return m_impl.get< extrinsicMeshData::slurryfluid::density >();
  }

  auto viscosity()
  {
    return m_impl.get< extrinsicMeshData::slurryfluid::viscosity >();
  }

private:
  StencilAccessors< extrinsicMeshData::slurryfluid::density,
                    extrinsicMeshData::slurryfluid::viscosity > m_impl;
};

} // end of namespace ProppantTransportKernels

} // end of namespace geosx

#endif //GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXSLURRYFLUIDACCESSORS_HPP
