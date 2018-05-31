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
 * @file FieldRegistry.cpp
 * @author settgast1
 * @date Nov 10, 2010
 */

#include "FieldRegistry.h"

FieldRegistry::FieldRegistry():
  m_objectSize(0),
  m_numFields(0),
  m_Registry(),
  m_fieldEnumToOffset()
{
  // TODO Auto-generated constructor stub

}

FieldRegistry::~FieldRegistry()
{
  // TODO Auto-generated destructor stub
}

void FieldRegistry::SetRegistry( const std::map< FieldKey, FieldType >& requestedFields )
{

  // set number of fields
  m_numFields = requestedFields.size();

  FieldRecord TempFieldRecord;

  // loop over requested fields
  for(std::map< FieldKey, FieldType >::const_iterator field = requestedFields.begin() ;
      field != requestedFields.end() ;
      ++field )
  {

    TempFieldRecord.m_fieldKey = field->first;
    TempFieldRecord.m_fieldType = field->second;
    TempFieldRecord.m_fieldSize = FieldSize(TempFieldRecord.m_fieldType);

    m_Registry.push_back(TempFieldRecord);
  }


  int offset = 1;
  for( std::vector<FieldRecord>::iterator field=m_Registry.begin() ; field!= m_Registry.end() ; ++field )
  {
    field->m_fieldOffset = offset;
    offset += field->m_fieldSize;
  }
  m_objectSize = offset;

  for( int i=0 ; i<m_numFields ; ++i )
  {
    m_fieldEnumToOffset[ m_Registry[i].m_fieldKey ] = m_Registry[i].m_fieldOffset;
  }

}
