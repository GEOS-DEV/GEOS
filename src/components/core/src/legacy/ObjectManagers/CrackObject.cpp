/*
 * CrackSurface.cpp
 *
 *  Created on: Nov 4, 2014
 *      Author: annavarapusr1, rrsettgast
 */

#include "CrackObject.h"





CrackObjectManager::CrackObjectManager():
ObjectDataStructureBaseT( ObjectDataStructureBaseT::CrackSurfaceManager ),
m_crackSurface1ToCutFaces(m_VariableOneToManyMaps["crackSurface1ToCutFaces"]),
m_crackSurface2ToCutFaces(m_VariableOneToManyMaps["crackSurface2ToCutFaces"]),
m_crackToElements(m_VariableOneToManyMaps["crackToElements"])
{
//  // For now use this definition to define a crack
//  this->AddKeylessDataField<R1Tensor>("crackPositionStart");
//  this->AddKeylessDataField<R1Tensor>("crackPositionEnd");
//  m_crackFaceManager.resize(1);

}

CrackObjectManager::~CrackObjectManager()
{
  // TODO Auto-generated destructor stub
}

