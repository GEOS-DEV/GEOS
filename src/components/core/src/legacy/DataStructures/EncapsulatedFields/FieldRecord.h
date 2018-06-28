/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file FieldRecord.h
 * @author settgast1
 * @date Nov 16, 2010
 */

#ifndef FIELDRECORD_H_
#define FIELDRECORD_H_

#include "Common/Common.h"

class FieldRecord
{
public:
  FieldRecord():
    m_fieldOffset(0),
    m_fieldSize(0),
    m_fieldKey(0),
    m_fieldType(FieldInfo::integerField)
  { }

  FieldRecord(const FieldRecord& init)
  {
    this->operator=(init);
  }

  FieldRecord& operator=(const FieldRecord& init)
  {
    this->m_fieldOffset = init.m_fieldOffset;
    this->m_fieldSize   = init.m_fieldSize;
    this->m_fieldKey   = init.m_fieldKey;
    this->m_fieldType   = init.m_fieldType;

    return (*this);
  }

  ~FieldRecord()
  { }


  int m_fieldOffset;
  int m_fieldSize;
  FieldKey m_fieldKey;
  FieldType m_fieldType;


};

#endif /* FIELDRECORD_H_ */
