/**
 * @file FractalVolume.h
 *
 *  Created on: June 23, 2012
 *  Author: scottjohnson
 */

#ifndef FRACTALVOLUME_H_
#define FRACTALVOLUME_H_

#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"
#include "TensorT/R1TensorT.h"
#include "StatisticalDistributionBaseT.h"
#include "Interpolation/VolumeKernel.h"
#include "FractalBaseT.h"

/**
 * @author Scott Johnson
 * @brief FractalVolume creates self-similar distributions of apertures
 */
class FractalVolume : public FractalBaseT
{
public:
  FractalVolume();

  unsigned InitializeFractal(const realT mean, const realT stdev,
                             const R1Tensor& lower, const R1Tensor& upper,
                             const realT hurst = 1.3, const localIndex nlevels = 6,
                             const localIndex n0 = 1, const localIndex n1 = 1,
                             const localIndex n2 = 1, const realT hfct = 1.2,
                             const unsigned seed = 0);

  unsigned Initialize(const Array2dT<realT>& parameters,
                      const R1Tensor& lower, const R1Tensor& upper,
                      const localIndex n0 = 1, const localIndex n1 = 1,
                      const localIndex n2 = 1, const realT hfct = 1.2,
                      const unsigned seed = 0);

  realT Value(const R1Tensor& position) const;

  localIndex Positions(const realT dx,
                       Array1dT<R1Tensor>& positions,
                       const realT weight = 1.0) const;

  localIndex Positions(const R1Tensor& min,
                       const R1Tensor& max,
                       const realT dx,
                       Array1dT<R1Tensor>& positions,
                       const realT weight = 1.0) const;

private:
  void InitializeSumMW(const int nj, const int nk, const int nl, const int ioffset,
                       const localIndex ilevel, Array3dT<VolumeKernel*>& curr);

  void InitializeLevel(const int nj, const int nk, const int nl, const int ioffset,
                       const realT dx0, const realT dx1, const realT dx2,
                       const realT mean, const realT stdev, const realT h,
                       Array3dT<VolumeKernel*>& curr, localIndex& icurr);

  void FillValues(const int ioffset, const localIndex nlevels,
                  Array1dT< Array3dT<VolumeKernel*>* >& vals);

  ///lower 3D coordinates of the points that will be queried
  R1Tensor m_lower;

  ///upper 3D coordinates of the points that will be queried
  R1Tensor m_upper;

  ///Initial value of n2
  localIndex m_n2;

  ///For each cell at the finest level, holds, for each level, the kernels applicable to the cell
  Array3dT<Array1dT<Array1dT<VolumeKernel*> > > m_values;

  ///Kernels referenced in values
  Array1dT<VolumeKernel> m_kernels;
};

#endif /* FRACTALVOLUME_H_ */
