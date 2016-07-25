/*
 * FractalVolume.cpp
 *
 *  Created on: June 23, 2012
 *  Author: scottjohnson
 */

#include "FractalVolume.h"
#include <limits.h>

/**
 * @brief Default constructor for the aperture generator class
 * @author Scott Johnson
 */
FractalVolume::FractalVolume():
m_lower(std::numeric_limits<realT>::max()),
m_upper(-std::numeric_limits<realT>::max()),
m_n2(0),
m_values(),
m_kernels()
{
}

/**
 * @brief Creates an oct-tree estimate of a fractally self-similar distribution of volumetric parameters
 * @author Scott Johnson
 * @param mean Mean of the parameter distribution
 * @param stdev Standard deviation of the parameter distribution
 * @hurst hurst Hurst exponent for the fractal distribution $1<H<2$ where $D=3-H$
 * @param lower Lower bound of the point values to be queried
 * @param upper Upper bound of the point values to be queried
 * @param nlevels Number of oct-tree levels to use to represent the multi-scale distribution
 * @param n0 Number of oct-tree elements in the x-direction at level 0
 * @param n1 Number of oct-tree elements in the y-direction at level 0
 * @param n2 Number of oct-tree elements in the z-direction at level 0
 * @param hfct Smoothing length multiplication factor
 * @return return
 */
unsigned FractalVolume::InitializeFractal(const realT mean, const realT stdev,
                                          const R1Tensor& lower, const R1Tensor& upper,
                                          const realT hurst, const localIndex nlevels,
                                          const localIndex n0, const localIndex n1,
                                          const localIndex n2, const realT hfct,
                                          const unsigned seed)
{
  Array2dT<realT> parameters(nlevels, 2);
  FillFractalParameters(mean, stdev, hurst, nlevels, parameters);

  //call the atomic initialize
  return Initialize(parameters, lower, upper, n0, n1, n2, hfct, seed);
}

unsigned FractalVolume::Initialize(const Array2dT<realT>& parameters,
                                   const R1Tensor& lower, const R1Tensor& upper,
                                   const localIndex n0, const localIndex n1,
                                   const localIndex n2, const realT hfct,
                                   const unsigned seed)
{
  if(parameters.Dimension(1) < 2)
    throw GPException("For normally distributed values, you must have at least two parameters: mean and standard deviation respectively");

  //initialize member variables
  m_hfct = hfct;
  m_lower = lower;
  m_upper = upper;
  m_n0 = n0;
  m_n1 = n1;
  m_n2 = n2;

  //make sure the random seed is appropriately set
  const unsigned sret = StatisticalDistributionBaseT::InitializeRandom(seed);

  //for the cubic interpolant, the compact support is hfct*2 cells
  const int ioffset = (int)(hfct * 2) + 1;
  const localIndex nlevels = parameters.Dimension(0);

  //resize the kernels array to its final size, so we can reference the memory without worrying!
  {
    size_t sz = 0;

    int aa = 1;
    for(localIndex i = 0; i < nlevels; i++, aa *= 2)
    {
      const int nj = 2*ioffset + n0 * aa;
      const int nk = 2*ioffset + n1 * aa;
      const int nl = 2*ioffset + n2 * aa;
      sz += nj * nk * nl;
    }
    m_kernels.resize(sz);
  }

  //fill the intermediate value hierarchical grid
  Array1dT< Array3dT<VolumeKernel*>* > vals;
  vals.resize(nlevels);
  {
    realT dx0 = 1.0 / n0;
    realT dx1 = 1.0 / n1;
    realT dx2 = 1.0 / n2;

    int aa = 1;
    localIndex icurr = 0;
    for(localIndex i = 0; i < nlevels; i++, aa *= 2, dx0 *= 0.5, dx1 *= 0.5, dx2 *= 0.5)
    {
      vals[i] = new Array3dT<VolumeKernel*>();
      const int nj = 2*ioffset + n0 * aa;
      const int nk = 2*ioffset + n1 * aa;
      const int nl = 2*ioffset + n2 * aa;

      //set the temporary statistical parameters to those for the current level
      const realT mean = parameters(i,0);
      const realT stdev = parameters(i,1);

      //set the smoothing length to that of the current level
      const realT h = (dx0 < dx1 ? (dx2 > dx1 ? dx2 : dx1) : (dx2 > dx0 ? dx2 : dx0)) * hfct;

      //setup new 3D array of surface kernel references
      Array3dT<VolumeKernel*>& curr = *vals[i];
      curr.Allocate(nj, nk, nl);
      curr = 0;

      //initialize the kernels and calculate the mW sums
      InitializeLevel(nj, nk, nl, ioffset, dx0, dx1, dx2, mean, stdev, h, curr, icurr);
    }
  }

  //fill the final list for each cell
  FillValues(ioffset, nlevels, vals);

  for(localIndex i = 0; i < nlevels; i++)
  {
    delete vals[i];
  }

  return sret;
}

void FractalVolume::InitializeLevel(const int nj, const int nk, const int nl, const int ioffset,
                                    const realT dx0, const realT dx1, const realT dx2,
                                    const realT mean, const realT stdev, const realT h,
                                    Array3dT<VolumeKernel*>& curr, localIndex& icurr)
{
  //calculate initial values and set position
  {
    R1Tensor xx(0.0);
    for(int j = 0; j < nj; j++)
    {
      xx(0) = dx0 * ((j - ioffset) + 0.5);
      for(int k = 0; k < nk; k++)
      {
        xx(1) = dx1 * ((k - ioffset) + 0.5);
        for(int l = 0; l < nl; l++)
        {
          xx(2) = dx2 * ((l - ioffset) + 0.5);
          m_kernels[icurr].Initialize(xx, h, mean, stdev);
          curr(j, k, l) = &m_kernels[icurr++];
        }
      }
    }
  }

  //Get point-wise density
  //InitializeSumMW(nj, nk, nl, ioffset, ilevel, curr);
}


void FractalVolume::InitializeSumMW(const int nj, const int nk, const int nl, const int ioffset,
                                    const localIndex, Array3dT<VolumeKernel*>& curr)
{
  //calculate sum_mW ... only need to do this once, since it's a regular grid
  R1Tensor dx(0.0);
  for(int j = 0; j < nj; j++)
  {
    for(int k = 0; k < nk; k++)
    {
      for(int l = 0; l < nl; l++)
      {
        R1Tensor& xjkl = curr(j, k, l)->m_x;
        for(int jj = j-ioffset; jj <= j+ioffset; jj++)
        {
          for(int kk = k-ioffset; kk <= k+ioffset; kk++)
          {
            for(int ll = l-ioffset; ll <= l+ioffset; ll++)
            {
              if(jj < 0 || jj >= nj || kk < 0 || kk >= nk || ll < 0 || ll >= nl)
                continue;
              dx = xjkl;
              dx -= curr(jj, kk, ll)->m_x;
              const realT dxd2 = Dot(dx, dx);
              curr(j, k, l)->IncrementSum(curr(jj, kk, ll)->mW(dxd2));
            }
          }
        }
        //std::cout << "kernel " << curr(j,k)->m_x[0] << " " << curr(j,k)->m_x[1] << " " << ilevel << " " << curr(j,k)->IncrementSum(0) << std::endl;
      }
    }
  }
}


void FractalVolume::FillValues(const int ioffset, const localIndex nlevels,
                               Array1dT< Array3dT<VolumeKernel*>* >& vals)
{
  //note: values will only hold cells away from the compact support
  //informed boundary, so make sure to remove the indices of
  //boundary cells
  const int nj = vals.back()->Dimension(0)-(2*ioffset);
  const int nk = vals.back()->Dimension(1)-(2*ioffset);
  const int nl = vals.back()->Dimension(2)-(2*ioffset);

  //set the size of the values array
  {
    int nAll = 2*ioffset + 1;
    nAll *= (nAll * nAll);
    m_values.clear();
    m_values.Allocate(nj, nk, nl);
    for(localIndex i = 0; i < m_values.size(); i++)
    {
      m_values[i].resize(nlevels);
      for(localIndex j = 0; j < nlevels; j++)
      {
          m_values[i][j].reserve(nAll);
      }
    }
  }

  for( int j = 0; j < nj; j++)  //for each cell at the finest level - j
  {
    for( int k = 0; k < nk; k++)  //for each cell at the finest level - k
    {
      for( int l = 0; l < nl; l++)  //for each cell at the finest level - l
      {
        Array1dT<Array1dT<VolumeKernel*> >& currentVal = m_values(j, k, l); // get the list of kernels for each level
        int nper = 1;
        for (localIndex i = 1; i < nlevels; i++)
          nper *= 2;
        for (localIndex i = 0; i < nlevels; i++, nper /= 2)
        {
          Array3dT<VolumeKernel*>& curr = *vals[i];
          //get the coordinates of the parent cell at level i
          int jj0 = j / nper;
          int jj1 = jj0 + 2 * ioffset;
          int kk0 = k / nper;
          int kk1 = kk0 + 2 * ioffset;
          int ll0 = l / nper;
          int ll1 = ll0 + 2 * ioffset;
          for (int jj = jj0; jj <= jj1; jj++) {
            for (int kk = kk0; kk <= kk1; kk++) {
              for (int ll = ll0; ll <= ll1; ll++) {
                currentVal[i].push_back(curr(jj, kk, ll));
              }
            }
          }
        }
      }
    }
  }
}

realT FractalVolume::Value(const R1Tensor& position) const
{
  //get the normalized position to evaluate in the hierarchical grid
  R1Tensor xt(position);
  {
    R1Tensor tmp2(m_upper);
    tmp2 -= m_lower;
    for (localIndex i=0; i<nsdof; ++i)
    {
      if (isEqual(tmp2[i],0.0))
      {
        xt[i] = 0.5;
      }
      else
      {
        xt[i] -= m_lower[i];
        xt[i] /= tmp2[i];
      }
    }
  }

  //get the bin of the normalized position and then evaluate the position for the kernels applicable at each level
  int coords[nsdof];
  for(localIndex i = 0; i < nsdof; i++)
    coords[i] = (int)(m_values.Dimension(i) * (xt(i) >= 1 ? (1.0-1.0e-10) : xt(i)) );
  const Array1dT<Array1dT<VolumeKernel*> >& valArray = m_values(coords[0], coords[1], coords[2]);

  //for each refinement level and then for each applicable kernel at that refinement level add the contribution
  {
    realT sum = 0.0;
    for(localIndex i = 0; i < valArray.size(); i++)
    {
      realT fmW = 0.0;
      realT mW = 0.0;
      for(localIndex j = 0; j < valArray[i].size(); j++)
      {
        valArray[i][j]->Evaluate(xt, fmW, mW); //add kernel contribution to total
      }
      sum += fmW / mW;
    }
    return sum;
  }
}

localIndex FractalVolume::Positions(const realT dx,
                                 Array1dT<R1Tensor>& positions,
                                 const realT weight) const
{
  return Positions(m_lower, m_upper, dx, positions, weight);
}

localIndex FractalVolume::Positions(const R1Tensor& min,
                                    const R1Tensor& max,
                                    const realT dx,
                                    Array1dT<R1Tensor>& positions,
                                    const realT weight) const
{
  //------GET INCREMENTAL DISTANCE ALONG EACH DIRECTION
  R1Tensor dxc;
  lArray1d nn(nsdof);
  {
    R1Tensor dxv(m_upper);
    dxv -= m_lower;
    dxc = dxv;
    dxc /= dx;
    for(localIndex i = 0; i < nsdof; ++i)
    {
      nn(i) = (localIndex)dxc(i);
      if(nn(i) == 0)
        ++nn(i);
      dxc(i) = dxv(i) / nn(i);
    }
  }

  //------GET INCREMENTAL VOLUME
  const realT dV = dxc(0) * dxc(1) * dxc(2);
  R1Tensor dxch(dxc);
  dxch *= 0.5;

  //------GET DISTRIBUTION
  R1Tensor x;
  for(localIndex ii = 0; ii < nn(0); ii++)
  {
    x(0) = m_lower(0) + dxc(0) * (ii + 0.5);
    for(localIndex jj = 0; jj < nn(1); jj++)
    {
      x(1) = m_lower(1) + dxc(1) * (jj + 0.5);
      for(localIndex kk = 0; kk < nn(2); kk++)
      {
        x(2) = m_lower(2) + dxc(2) * (kk + 0.5);
        const realT val = weight * dV * Value(x);
        const realT mval = fmod(val, 1.0);
        const localIndex njoints = ((localIndex) (val - mval)) + (mval > StatisticalDistributionBaseT::UniformSample(0, 1.0) ? 1 : 0);

        R1Tensor xmax(x);
        xmax += dxch;
        R1Tensor xmin(x);
        xmin -= dxch;

        for(localIndex k = 0; k < njoints; k++)
        {
          R1Tensor xx;
          for(localIndex iii = 0; iii < nsdof; ++iii)
            xx(iii) = StatisticalDistributionBaseT::UniformSample(xmin(iii), xmax(iii));
          positions.push_back(xx);
        }
      }
    }
  }
  return positions.size();
}
