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
 * @file LinearTriangle.h
 * @author Fu, Pengcheng
 * @date July 3, 2012
 */

#include "legacy/ElementLibrary/FiniteElement.h"

#ifndef LINE_H_
#define LINE_H_

class Line : public FiniteElement<2>
{
public:
  Line();
  virtual ~Line();
  void reinit(const std::vector<R1TensorT<3> > &mapped_support_points);

};

#endif /* LINE_H_ */
