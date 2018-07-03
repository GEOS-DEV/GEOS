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
 * @file ObjectManager.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */


#ifndef OBJECTMANAGER_H_
#define OBJECTMANAGER_H_

#include "Common/Common.h"
#include "Object.h"
#include <map>
#include "FieldRegistry.h"


/**
 * @author Randolph Settgast
 *
 * The ObjectManager
 */
class ObjectManager
{
public:

  /// default constructor
  ObjectManager();

  /// copy constructor
  ObjectManager( const ObjectManager& init );

  /// default destructor
  ~ObjectManager();


  std::vector< Object* > m_objects;



  size_t m_numObjects;
  dvector m_objectMemory;


  std::map< FieldKey, FieldType > m_requestedFields;

  FieldRegistry m_fieldRegistry;

  virtual void resize( const size_t size );

  void RegisterFields( std::vector< std::pair< FieldKey, FieldType > > requestedFields );

  void OrganizeFields(  );

  void AllocateObjectFields();



protected:

private:
  ObjectManager& operator=( const ObjectManager&);

};



#endif /* ObjectManager_H_ */
