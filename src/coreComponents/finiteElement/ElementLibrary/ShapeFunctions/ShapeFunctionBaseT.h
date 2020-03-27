/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
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

  ShapeFunctionBaseT( void );

  virtual ~ShapeFunctionBaseT( void );

  void Calc_Shape_Deriv( const realT fac );


  realT ShapeFunctionValue( const R1Tensor & Xi,
                            const R1Tensor & Xi_node );

  int CalculateJacobian( const int elem );



  //***** Data Member Accessors **********************************************
public:



};



#endif
