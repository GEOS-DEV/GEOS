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
