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
 * TempODS.h
 *
 *  Created on: Jan 17, 2013
 *      Author: settgast1
 */

#ifndef TEMPODS_H_
#define TEMPODS_H_

#include "ObjectDataStructureBaseT.h"

class TempODS : public ObjectDataStructureBaseT
{
public:
  TempODS();
  virtual ~TempODS();

  void Initialize(  ) {}


  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<gArray1d>& objectToCompositionObject ){}

  /// pure virtual function that sets what objects are on the boundary of the
  // domain
  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL){}

  /// pure virtual function that sets what objects are external
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL){}

};

#endif /* TEMPODS_H_ */
