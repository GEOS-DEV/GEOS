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
 * aperturatorlib.h
 *
 *  Created on: Jun 5, 2012
 *      Author: johnson346
 */

#ifndef FRACTALSURFACELIB_H_
#define FRACTALSURFACELIB_H_

#include "FractalSurface.h"

int GenerateFractalSurface(double loweru, double lowerv,
                           double upperu, double upperv,
                           double hurst, double hfct,
                           double mean, double stdev,
                           int nlevels, int n0, int n1);

double ValueFractalSurface(double u, double v);

#endif /* FRACTALSURFACELIB_H_ */
