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
 * @file CellBase.cpp
 */

#include "ElementSubRegionBase.hpp"

namespace geosx
{
using namespace dataRepository;

ElementSubRegionBase::ElementSubRegionBase( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_constitutiveModels( groupKeyStruct::constitutiveModelsString(), this ),
  m_numNodesPerElement(),
  m_numEdgesPerElement(),
  m_numFacesPerElement(),
  m_elementCenter(),
  m_elementVolume()
{
  registerGroup( groupKeyStruct::constitutiveModelsString(), &m_constitutiveModels ).
    setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::numNodesPerElementString(), &m_numNodesPerElement );

  registerWrapper( viewKeyStruct::numEdgesPerElementString(), &m_numEdgesPerElement );

  registerWrapper( viewKeyStruct::numFacesPerElementString(), &m_numFacesPerElement );

  registerWrapper( viewKeyStruct::elementCenterString(), &m_elementCenter ).
    setPlotLevel( PlotLevel::LEVEL_1 ).
    reference().resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::elementVolumeString(), &m_elementVolume ).
    setPlotLevel( PlotLevel::LEVEL_1 );
}

ElementSubRegionBase::~ElementSubRegionBase()
{}

} /* namespace geosx */
