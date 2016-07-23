/**
 * @file InterpolantBaseT.h
 *
 *  Created on: September 3, 2013
 *  Author: scottjohnson
 */

#ifndef INTERPOLANTBASET_H_
#define INTERPOLANTBASET_H_

#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"
#include "SurfaceGeneration/StatisticalDistributionBaseT.h"

/**
 * @author Scott Johnson
 * @brief InterpolantBaseT holds the definition of a generic interpolation kernel with compact support
 */
class InterpolantBaseT
{
public:
  InterpolantBaseT() : m_value(0.0), m_sum_mW(0.0), m_invhd(0.0), m_invh2(0.0)
  {
  }

  inline realT Value() const
  {
    return this->m_value;
  }

  inline realT Mass() const
  {
    return this->m_mass;
  }

  inline realT CurrentSum() const
  {
    return this->m_sum_mW;
  }

  inline void ZeroSum()
  {
    this->m_sum_mW = 0.0;
  }

  inline realT IncrementSum(const realT mwij)
  {
    this->m_sum_mW += mwij;
    return this->m_sum_mW;
  }

  virtual realT SetInverseSmoothingLengthFactor(const realT h) = 0;

protected:

  static inline realT M4(const realT s, const realT alpha)
  {
    if(s > 2)
      return 0.0;
    const static realT two_thirds = 2.0/3.0;
    const static realT one_sixth = 1.0/6.0;
    if(s >= 1)
    {
      realT w = 2.0 - s;
      w *= alpha * w * w * one_sixth;
      return w;
    }
    else
    {
      const realT w = alpha * (two_thirds - s * s + 0.5 * s * s * s);
      return w;
    }
  }

  static inline realT dM4(const realT s, const realT alpha)
  {
    if(s > 2)
      return 0.0;
    const realT two_min_s = 2 - s;
    const realT ret = alpha * (s < 1 ? -2.0 + 1.5 * s : -0.5 * two_min_s * two_min_s / s);
    return ret;
  }

public:
  inline realT W4(const realT dxd2, const realT alpha) const
  {
    const realT s = sqrt(dxd2 * m_invh2);
    return M4(s, alpha) * m_invhd;
  }

  inline realT dW4(const realT dxd2, const realT alpha) const
  {
    const realT s = sqrt(dxd2 * m_invh2);
    return dM4(s, alpha) * m_invhd * m_invh2;
  }

protected:
  realT Evaluate(const realT dxd2, const realT alpha) const
  {
    const realT s2 = dxd2 * m_invh2;
    const realT w = W4(s2, alpha);
    return w * (this->m_value / this->m_sum_mW);//fW/rho
  }

  static inline realT Alpha(const localIndex dim) {
    if(dim > 3 || dim < 1)
      throw GPException("No interpolant normalization factor for dimension > 3");
    return dim == 1 ? 1.0/6.0 : 0.31830988618 * (dim==2 ? 15.0/14.0 : 0.25);
  }

protected:
  realT m_mass;
  realT m_value;
  realT m_sum_mW;
  realT m_invhd;
  realT m_invh2;
};

#endif /* INTERPOLANTBASET_H_ */
