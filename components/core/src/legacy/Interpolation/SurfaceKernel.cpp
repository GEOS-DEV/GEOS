/*
 * SurfaceKernel.cpp
 *
 *  Created on: Sep 5, 2013
 *  Author: scottjohnson
 */

#include "Interpolation/SurfaceKernel.h"

void SurfaceKernel::Evaluate(const R1TensorT<2>& x,
                             realT& fmW4, realT& mW4) const
{
  R1TensorT<2> dx(x);
  dx(0) -= m_x(0);
  dx(1) -= m_x(1);
  const realT w = InterpolantBaseT::W4( Dot(dx,dx), Alpha());
  fmW4 += m_value * m_mass * w;
  mW4 += m_mass * w;
}

realT SurfaceKernel::SetInverseSmoothingLengthFactor(const realT h)
{
  m_invh2 = 1.0 / (h * h);
  m_invhd = m_invh2;
  return m_invhd;
}

void SurfaceKernel::Initialize(const R1TensorT<2>& x, const realT h, const realT mean, const realT stdev)
{
  m_mass = h*h;
  m_x = x;
  SetInverseSmoothingLengthFactor(h);
  m_value = StatisticalDistributionBaseT::NormalSample(mean, stdev);
  m_sum_mW = 0.0;
}
