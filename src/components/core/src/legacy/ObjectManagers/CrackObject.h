/*
 * CrackSurface.h
 *
 *  Created on: Nov 4, 2014
 *      Author: annavarapusr1, rrsettgast
 */

//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "../../dataRepository/Group.hpp"

#ifndef CRACKOBJECT_H_
#define CRACKOBJECT_H_


class CrackFaceManager: public ObjectDataStructureBaseT
{
public:
  CrackFaceManager():
    ObjectDataStructureBaseT( ObjectDataStructureBaseT::CrackFaceManager )
  {
    this->AddKeyedDataField<FieldInfo::referencePosition>();
    this->AddKeylessDataField<R1Tensor>("normal");
  }
  ~CrackFaceManager();

};

class CrackObjectManager: public ObjectDataStructureBaseT
{
public:
  CrackObjectManager();
  virtual ~CrackObjectManager();



  virtual void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                 Array1dT<gArray1d>& objectToCompositionObject )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("CrackObject::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }

  virtual void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL) {}

  virtual void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL) {}

  virtual void Initialize(  ) {}

//  UnorderedVariableOneToManyRelation& m_crackSurface1ToCutFaces;
//  UnorderedVariableOneToManyRelation& m_crackSurface2ToCutFaces;

  OrderedVariableOneToManyRelation& m_crackSurface1ToCutFaces;
  OrderedVariableOneToManyRelation& m_crackSurface2ToCutFaces;
  OrderedVariableOneToManyRelation& m_crackToElements;

  Array1dT<rArray2d> m_crackNodeCoords;
  Array1dT<iArray2d> m_crackConn;

//  UnorderedVariableOneToManyRelation& m_crackSurfaceToCrackFaces;

};

#endif /* CRACKOBJECT_H_ */
