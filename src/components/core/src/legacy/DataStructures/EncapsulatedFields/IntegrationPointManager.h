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
 * @file IntegrationPointManager.h
 * @author settgast1
 * @date Dec 3, 2010
 */

#ifndef INTEGRATIONPOINTMANAGER_H_
#define INTEGRATIONPOINTMANAGER_H_

#include "ObjectManager.h"

class IntegrationPointManager : public ObjectManager
{
public:
  IntegrationPointManager();
  ~IntegrationPointManager();

  Array2dT< array<R1Tensor> >  m_dNdX;
  rArray2d m_detJ;

};

#endif /* INTEGRATIONPOINTMANAGER_H_ */
