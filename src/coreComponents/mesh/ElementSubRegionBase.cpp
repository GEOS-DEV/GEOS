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
* @file CellBase.cpp
*/

#include "ElementSubRegionBase.hpp"

namespace geosx
{
using namespace dataRepository;

ElementSubRegionBase::ElementSubRegionBase( string const & name, Group * const parent ):
ObjectManagerBase(name,parent),
m_constitutiveModels(groupKeyStruct::constitutiveModelsString,this),
m_numNodesPerElement(),
m_numEdgesPerElement(),
m_numFacesPerElement(),
m_elementCenter(),
m_elementVolume()
{
RegisterGroup( groupKeyStruct::constitutiveModelsString, &m_constitutiveModels, 0 );
}

ElementSubRegionBase::~ElementSubRegionBase()
{
}

void ElementSubRegionBase::SetElementType( string const & elementType )
{
m_elementTypeString = elementType;
m_elementType =FiniteElementBase::StringToElementType( elementType );
}



} /* namespace geosx */
