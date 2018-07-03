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
 * @file ElementRegion.h
 * @author settgast1
 * @date Dec 2, 2010
 */

#ifndef ELEMENTREGION_H_
#define ELEMENTREGION_H_

#include "ObjectManager.h"



class ElementRegion : public ObjectManager
{
public:
  ElementRegion();
  ~ElementRegion();

  const localIndex& m_numElements;

  array<R1Tensor> xi;

  array<R1Tensor>  m_dNdXi;


  Array2dT< array<R1Tensor> >  m_dNdX;

  rArray2d detJ;

  ivector m_MapLocalToGlobal;


  std::vector< std::vector<Object*> > m_Nodes;


};

#endif /* ELEMENTREGION_H_ */
