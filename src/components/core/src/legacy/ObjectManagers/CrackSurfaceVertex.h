/*
 * CrackSurfaceVertex.h
 *
 *  Created on: Nov 4, 2014
 *      Author: rrsettgast
 */

#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#ifndef CRACKSURFACEVERTEX_H_
#define CRACKSURFACEVERTEX_H_

class CrackSurfaceVertex : public ObjectDataStructureBaseT
{
public:
  CrackSurfaceVertex();
  virtual ~CrackSurfaceVertex();

  virtual void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                 array<gArray1d>& objectToCompositionObject )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("CrackSurfaceVertex::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }

  virtual void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL) {}

  virtual void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL) {}

  virtual void Initialize(  ) {}

  array<R1Tensor> Vertices;

};

#endif /* CRACKSURFACEVERTEX_H_ */
