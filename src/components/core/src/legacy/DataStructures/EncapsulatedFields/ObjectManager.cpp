// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * ObjectManager.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "ObjectManager.h"

/**
 * @author R.R. Settgast
 * @return
 */
ObjectManager::ObjectManager():
  m_objects(),
  m_numObjects(0),
  m_objectMemory(),
  m_requestedFields(),
  m_fieldRegistry()
{}


/**
 * @author R.R. Settgast
 * @return
 */
ObjectManager::ObjectManager( const ObjectManager& init )
{}


/**
 * @author R.R. Settgast
 * @return
 */
ObjectManager::~ObjectManager()
{}


void ObjectManager::RegisterFields( std::vector< std::pair< FieldKey, FieldType > > requestedFields )
{
  for(std::vector< std::pair< FieldKey, FieldType > >::const_iterator field = requestedFields.begin() ;
      field != requestedFields.end() ; ++field )
  {
    m_requestedFields.insert( *field );
  }

}

void ObjectManager::OrganizeFields(  )
{
  m_fieldRegistry.SetRegistry(m_requestedFields);
}


void ObjectManager::AllocateObjectFields()
{
  m_objectMemory.resize( m_fieldRegistry.m_objectSize * m_numObjects );
}

void ObjectManager::resize( const size_t size )
{
  m_numObjects = size;
  m_objects.resize( size );


  this->AllocateObjectFields();

  realT* p_object = &(m_objectMemory[0]);
  for( std::vector<Object*>::iterator object=m_objects.begin() ; object!=m_objects.end() ; ++object )
  {
    (*object) = (Object*)p_object;
    (*object)->m_fieldRegistry = &(this->m_fieldRegistry);
    p_object += m_fieldRegistry.m_objectSize;
  }

}
