// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file Object.cpp
 * @author settgast1
 * @date Nov 5, 2010
 */

#include "Object.h"

Object::Object():
  m_fieldData(NULL),
  m_fieldRegistry(NULL)
{
  // TODO Auto-generated constructor stub

}

Object::Object(const Object& init)
{
  operator=(init);
}

Object::~Object()
{
  if( m_fieldData != NULL )
    delete[] m_fieldData;
}



void Object::resize( const size_t size )
{
  if( m_fieldRegistry->m_objectSize == size )
  {
    if( m_fieldData != NULL )
      delete[] m_fieldData;

    m_fieldRegistry->m_objectSize = size;
    m_fieldData = new realT[size];

  }
}
