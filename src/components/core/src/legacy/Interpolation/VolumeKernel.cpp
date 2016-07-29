/*
 * SurfaceKernel.cpp
 *
 *  Created on: Sep 5, 2013
 *  Author: scottjohnson
 */

#include "Interpolation/VolumeKernel.h"

void VolumeKernel::Evaluate(const R1Tensor& x,
                            realT& fmW4, realT& mW4) const
{
  R1Tensor dx(x);
  dx -= m_x;
  const realT w = InterpolantBaseT::W4( Dot(dx,dx), Alpha());
  fmW4 += m_value * m_mass * w;
  mW4 += m_mass * w;
}

realT VolumeKernel::SetInverseSmoothingLengthFactor(const realT h)
{
  m_invh2 = 1.0 / (h * h);
  m_invhd = m_invh2 / h;
  return m_invhd;
}

void VolumeKernel::Initialize(const R1Tensor& x, const realT h, const realT mean, const realT stdev)
{
  m_mass = h*h*h;
  m_x = x;
  SetInverseSmoothingLengthFactor(h);
  m_value = StatisticalDistributionBaseT::NormalSample(mean, stdev);
  m_sum_mW = 0.0;
}
