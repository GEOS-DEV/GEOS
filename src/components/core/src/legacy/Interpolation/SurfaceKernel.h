/**
 * @file SurfaceInterpolant.h
 *
 *  Created on: September 4, 2013
 *  Author: scottjohnson
 */

#ifndef SURFACEKERNEL_H_
#define SURFACEKERNEL_H_

#include "Interpolation/InterpolantBaseT.h"
#include "TensorT/R1TensorT.h"

/**
 * @author Scott Johnson
 * @brief SurfaceInterpolant holds the definition of a 2D interpolation kernel with compact support
 */
class SurfaceKernel : public InterpolantBaseT
{
public:
  SurfaceKernel(): m_x(0.0) {}

  void Initialize(const R1TensorT<2>& x, const realT h, const realT mean, const realT stdev);

  virtual realT SetInverseSmoothingLengthFactor(const realT h);

  void Evaluate(const R1TensorT<2>& x, realT& fmW4, realT& mW4) const;

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
    static realT alpha = InterpolantBaseT::Alpha(2);
    return alpha;
  }

public:
  R1TensorT<2> m_x;

};

#endif /* SURFACEKERNEL_H_ */
