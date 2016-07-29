/*
 * TempODS.h
 *
 *  Created on: Jan 17, 2013
 *      Author: settgast1
 */

#ifndef TEMPODS_H_
#define TEMPODS_H_

#include "ObjectDataStructureBaseT.h"

class TempODS: public ObjectDataStructureBaseT
{
public:
  TempODS();
  virtual ~TempODS();

  void Initialize(  ) {}


   void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                 Array1dT<gArray1d>& objectToCompositionObject ){}

  /// pure virtual function that sets what objects are on the boundary of the domain
   void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL){}

  /// pure virtual function that sets what objects are external
   void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL){}

};

#endif /* TEMPODS_H_ */
