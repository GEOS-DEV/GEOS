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
 * @file FieldRegistry.h
 * @author settgast1
 * @date Nov 10, 2010
 */

#ifndef FIELDREGISTRY_H_
#define FIELDREGISTRY_H_



#include <vector>
#include <map>

#include "FieldRecord.h"

class FieldRegistry
{
public:
  FieldRegistry();
  ~FieldRegistry();

  size_t m_objectSize;
  int m_numFields;

  std::vector<FieldRecord> m_Registry;


  int m_fieldEnumToOffset[FieldInfo::numFieldEnums];



  void SetRegistry( const std::map< FieldKey, FieldType >& requestedFields );

};



#endif /* FIELDREGISTRY_H_ */
