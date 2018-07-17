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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TriangleShell.h
 * @This element is implemented so that flow network has something to attach on.
 * @No mechanical behavior of real shell is implemented.
 * @m_elementGeometryID is "TRSH".
 * @author Fu, Pengcheng
 * @date Aug 12, 2012
 */

#include "legacy/ElementLibrary/FiniteElement.h"

#ifndef TRIANGLESHELL_H_
#define TRIANGLESHELL_H_

class TriangleShell : public FiniteElement<3>
{
public:
  TriangleShell();
  virtual ~TriangleShell();
  void reinit(const std::vector<R1TensorT<3> > &mapped_support_points);

};

#endif /* TRIANGLESHELL_H_ */
