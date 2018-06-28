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

/*
 * NodeManager.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "NodeManager.h"
#include <fstream>

/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager():
  ObjectManager(),
  m_numNodes(ObjectManager::m_numObjects)
{

  std::vector< std::pair< FieldKey, FieldType > > fields;
  fields.resize(7);

  fields[0] = std::pair< FieldKey, FieldType >(FieldInfo::referencePosition,FieldInfo::R1TensorField);
  fields[1] = std::pair< FieldKey, FieldType >(FieldInfo::displacement,FieldInfo::R1TensorField);
  fields[2] = std::pair< FieldKey, FieldType >(FieldInfo::incrementalDisplacement,FieldInfo::R1TensorField);
  fields[3] = std::pair< FieldKey, FieldType >(FieldInfo::velocity,FieldInfo::R1TensorField);
  fields[4] = std::pair< FieldKey, FieldType >(FieldInfo::acceleration,FieldInfo::R1TensorField);
  fields[5] = std::pair< FieldKey, FieldType >(FieldInfo::force,FieldInfo::R1TensorField);
  fields[6] = std::pair< FieldKey, FieldType >(FieldInfo::mass,FieldInfo::realField);

  this->RegisterFields(fields);
  this->OrganizeFields();
  this->AllocateObjectFields();
}


/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager( const NodeManager& init ):
  ObjectManager(init),
  m_numNodes(ObjectManager::m_numObjects)
{}


/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::~NodeManager()
{}

/**
 * @author Settgast
 * @param geometryStream open file stream to read nodal data from
 *
 * Needs to be replaced once we settle on a file format
 */
void NodeManager::ReadAsciiNodeInput( std::ifstream& geometryStream )
{

  int globalNodeIndex;
  int junk;

  const int xOffset = m_fieldRegistry.m_fieldEnumToOffset[FieldInfo::referencePosition];

  for( std::vector<Object*>::iterator pnode = m_objects.begin() ; pnode != m_objects.end() ; ++pnode )
  {
    Object* const node = *pnode;
    R1Tensor& m_refposition = node->GetFieldFromOffset<R1Tensor>(xOffset);
    geometryStream>>globalNodeIndex>>m_refposition(0)>>m_refposition(1)>>m_refposition(2)>>junk;

  }

}



/**
 * @author R. Settgast
 * @param destination local node number of destination node
 * @param source local node number of source node
 */
void NodeManager::CopyNode( const int destination, const int source )
{
  *(m_objects[destination]) = *(m_objects[source]);

}

/**
 * @author R. Settgast
 * @param nodenum local node number of node to tranlate
 * @param offset amount to translate node
 */
void NodeManager::TranslateNode( const int nodenum, const R1Tensor& offset  )
{
  m_objects[nodenum]->GetFieldFromOffset<R1Tensor>(m_fieldRegistry.m_fieldEnumToOffset[FieldInfo::referencePosition]);

}
