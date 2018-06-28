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
 * File: ShapeFunctionBaseT.h
 * Shape function Class
 *
 * created : RRS (09/14/2010)
 */
#ifndef _SHAPE_T_H_
#define _SHAPE_T_H_


#include "../../Common/Common.h"

class ShapeFunctionBaseT
{
public:

  ShapeFunctionBaseT(void);

  virtual ~ShapeFunctionBaseT(void);

  void Calc_Shape_Deriv(const realT fac);


  realT ShapeFunctionValue( const R1Tensor& Xi,
                            const R1Tensor& Xi_node );

  int CalculateJacobian( const int elem );



  //***** Data Member Accessors **********************************************
public:



};



#endif
