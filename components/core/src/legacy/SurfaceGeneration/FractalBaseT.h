/**
 * @file FractalBaseT.h
 *
 *  Created on: June 23, 2012
 *  Author: scottjohnson
 */

#ifndef FRACTALBASET_H_
#define FRACTALBASET_H_

#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"

/**
 * @author Scott Johnson
 * @brief FractalBaseT creates self-similar distributions of properties
 */
class FractalBaseT
{
public:
  FractalBaseT();

protected:
  void FillFractalParameters(const realT mean, const realT stdev,
                             const realT hurst, const localIndex nlevels,
                             Array2dT<realT>& parameters);

  ///Smoothing length multiplication factor
  realT m_hfct;

  ///Initial values of n0, n1, and n2
  localIndex m_n0, m_n1;
};

#endif /* FRACTALBASET_H_ */
