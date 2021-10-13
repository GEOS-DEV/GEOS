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
 * @file NumericalMethodsManager.cpp
 */

#include "NumericalMethodsManager.hpp"


namespace geosx
{
using namespace dataRepository;

NumericalMethodsManager::NumericalMethodsManager( string const & name, Group * const parent ):
  Group( name, parent ),
  m_finiteElementDiscretizationManager( groupKeysStruct::finiteElementDiscretizationsString(), this ),
  m_finiteVolumeManager( groupKeysStruct::finiteVolumeManagerString(), this )
{
  setInputFlags( InputFlags::OPTIONAL );

  this->registerGroup( groupKeysStruct::finiteElementDiscretizationsString(), &m_finiteElementDiscretizationManager );
  this->registerGroup( groupKeysStruct::finiteVolumeManagerString(), &m_finiteVolumeManager );
}

NumericalMethodsManager::~NumericalMethodsManager()
{}

Group * NumericalMethodsManager::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                              string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}



} /* namespace geosx */
