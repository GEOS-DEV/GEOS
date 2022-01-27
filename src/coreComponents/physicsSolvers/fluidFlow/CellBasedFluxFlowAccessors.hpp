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

#ifndef GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXFLOWACCESSORS_HPP
#define GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXFLOWACCESSORS_HPP

#include "FlowSolverBaseExtrinsicData.hpp"


namespace geosx
{

namespace ProppantTransportKernels
{

class CellBasedFluxFlowAccessorsImpl
{
public:
  CellBasedFluxFlowAccessorsImpl( ElementRegionManager const & elemManager,
                                  string const & solverName )
    : m_impl( elemManager, solverName )
  { }

  auto pressure()
  {
    return m_impl.get< extrinsicMeshData::flow::pressure >();
  }

  auto gravityCoefficient()
  {
    return m_impl.get< extrinsicMeshData::flow::gravityCoefficient >();
  }

  auto elementAperture()
  {
    return m_impl.get< extrinsicMeshData::elementAperture >();
  }

private:
  StencilAccessors< extrinsicMeshData::flow::pressure,
                    extrinsicMeshData::flow::gravityCoefficient,
                    extrinsicMeshData::elementAperture> m_impl;
};

} // end of namespace ProppantTransportKernels

} // end of namespace geosx

#endif //GEOSX_PROPPANTTRANSPORTKERNELS_CELLBASEDFLUXFLOWACCESSORS_HPP
