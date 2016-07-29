/**
 * @file VolumetricInterpolant.h
 *
 *  Created on: September 4, 2013
 *  Author: scottjohnson
 */

#ifndef VOLUMEKERNEL_H_
#define VOLUMEKERNEL_H_

#include "Interpolation/InterpolantBaseT.h"
#include "TensorT/R1TensorT.h"

/**
 * @author Scott Johnson
 * @brief VolumetricInterpolant holds the definition of a 3D interpolation kernel with compact support
 */
class VolumeKernel : public InterpolantBaseT
{
public:
  VolumeKernel(): m_x(0.0)
  {
  };

  void Initialize(const R1Tensor& x, const realT h, const realT mean, const realT stdev);

  virtual realT SetInverseSmoothingLengthFactor(const realT h);

  void Evaluate(const R1Tensor& x, realT& fmW4, realT& mW4) const;

  inline realT mW(const realT dxd2) const
  {
    return m_value * W(dxd2);
  }

private:
  inline realT W(const realT dxd2) const
  {
    return M4(dxd2 * m_invh2, Alpha()) * m_invhd;
  }

  static inline realT Alpha()
  {
    static realT alpha = InterpolantBaseT::Alpha(3);
    return alpha;
  }

public:
  R1Tensor m_x;

};

#endif /* VOLUMEKERNEL_H_ */
