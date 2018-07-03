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
