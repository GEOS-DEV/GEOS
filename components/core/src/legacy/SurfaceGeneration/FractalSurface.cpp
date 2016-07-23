/*
 * FractalSurface.cpp
 *
 *  Created on: Sep 27, 2011
 *  Author: scottjohnson
 */

#include "FractalSurface.h"
#include <limits.h>

/**
 * @brief Default constructor for the aperture generator class
 * @author Scott Johnson
 */
FractalSurface::FractalSurface():
m_lower(std::numeric_limits<realT>::max()),
m_upper(-std::numeric_limits<realT>::max()),
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
 * @param hfct Smoothing length multiplication factor
 * @return return
 */
unsigned FractalSurface::InitializeFractal(const realT mean, const realT stdev,
                                          const R1TensorT<2>& lower, const R1TensorT<2>& upper,
                                          const realT hurst, const localIndex nlevels,
                                          const localIndex n0, const localIndex n1,
                                          const realT hfct,
                                          const unsigned seed)
{
  Array2dT<realT> parameters(nlevels, 2);
  FillFractalParameters(mean, stdev, hurst, nlevels, parameters);

  //call the atomic initialize
  return Initialize(parameters, lower, upper, n0, n1, hfct, seed);
}

unsigned FractalSurface::Initialize(const Array2dT<realT>& parameters,
                                   const R1TensorT<2>& lower, const R1TensorT<2>& upper,
                                   const localIndex n0, const localIndex n1,
                                   const realT hfct,
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
      sz += nj * nk;
    }
    m_kernels.resize(sz);
  }

  //fill the intermediate value hierarchical grid
  Array1dT< Array2dT<SurfaceKernel*> > vals;
  vals.resize(nlevels);
  {
    realT dx0 = 1.0 / n0;
    realT dx1 = 1.0 / n1;

    int aa = 1;
    localIndex icurr = 0;
    for(localIndex i = 0; i < nlevels; i++, aa *= 2, dx0 *= 0.5, dx1 *= 0.5)
    {
      const int nj = 2*ioffset + n0 * aa;
      const int nk = 2*ioffset + n1 * aa;

      //set the temporary statistical parameters to those for the current level
      const realT mean = parameters(i,0);
      const realT stdev = parameters(i,1);

      //set the smoothing length to that of the current level
      const realT h = (dx0 < dx1 ? dx1 : dx0) * hfct;

      //setup new 3D array of surface kernel references
      Array2dT<SurfaceKernel*>& curr = vals[i];
      curr.resize2(nj, nk);
      for(localIndex ii = 0; ii < curr.Dimension(0); ii++)
        for(localIndex jj = 0; jj < curr.Dimension(1); jj++)
          curr(ii,jj) = 0;

      //initialize the kernels and calculate the mW sums
      InitializeLevel(nj, nk, ioffset, dx0, dx1, mean, stdev, h, i, curr, icurr);
    }
  }

  //fill the final list for each cell
  FillValues(ioffset, nlevels, vals);

  return sret;
}

void FractalSurface::InitializeLevel(const int nj, const int nk, const int ioffset,
                                     const realT dx0, const realT dx1,
                                     const realT mean, const realT stdev, const realT h,
                                     const localIndex,
                                     Array2dT<SurfaceKernel*>& curr, localIndex& icurr)
{
  //calculate initial values and set position
  R1TensorT<2> xx(0.0);
  for(int j = 0; j < nj; j++)
  {
    xx(0) = dx0 * ((j - ioffset) + 0.5);
    for(int k = 0; k < nk; k++)
    {
      xx(1) = dx1 * ((k - ioffset) + 0.5);
      m_kernels[icurr].Initialize(xx, h, mean, stdev);
      curr(j, k) = &m_kernels[icurr++];
      //std::cout << "kernel " << curr(j,k)->m_x[0] << " " << curr(j,k)->m_x[1] << " " << ilevel << " " << curr(j,k)->Value() << std::endl;
    }
  }

  //Get point-wise density
  //InitializeSumMW(nj, nk, ioffset, ilevel, curr);
}

void FractalSurface::InitializeSumMW(const int nj, const int nk, const int ioffset,
                                     const localIndex, Array2dT<SurfaceKernel*>& curr)
{
  //calculate sum_mW ... only need to do this once, since it's a regular grid
  R1TensorT<2> dx(0.0);
  for (int j = 0; j < nj; j++)
  {
    for (int k = 0; k < nk; k++)
    {
      curr(j, k)->ZeroSum();
      const R1TensorT<2>& xjk = curr(j, k)->m_x;
      for (int jj = j - ioffset; jj <= j + ioffset; jj++)
      {
        for (int kk = k - ioffset; kk <= k + ioffset; kk++)
        {
          if (jj < 0 || jj >= nj || kk < 0 || kk >= nk)
            continue;
          dx = xjk;
          dx -= curr(jj, kk)->m_x;
          const realT dxd2 = Dot(dx, dx);
          curr(j, k)->IncrementSum(curr(jj, kk)->mW(dxd2));
        }
      }
      //std::cout << "kernel " << curr(j,k)->m_x[0] << " " << curr(j,k)->m_x[1] << " " << ilevel << " " << curr(j,k)->IncrementSum(0) << std::endl;
    }
  }
}

void FractalSurface::FillValues(const int ioffset, const localIndex nlevels,
                                Array1dT< Array2dT<SurfaceKernel*> >& vals)
{
  //note: values will only hold cells away from the compact support
  //informed boundary, so make sure to remove the indices of
  //boundary cells
  const int nj = vals.back().Dimension(0) - (2 * ioffset);
  const int nk = vals.back().Dimension(1) - (2 * ioffset);

  //set the size of the values array
  {
    int nAll = 2 * ioffset + 1;
    nAll *= nAll;
    m_values.clear();
    m_values.resize2(nj, nk);
    for (localIndex i = 0; i < m_values.Dimension(0); i++)
    {
      for (localIndex j = 0; j < m_values.Dimension(1); j++)
      {
        m_values(i,j).resize(nlevels);
        for (localIndex k = 0; k < nlevels; k++)
        {
          m_values(i,j)[k].reserve(nAll);
        }
      }
    }
  }

  for (int j = 0; j < nj; j++) //for each cell at the finest level - j
  {
    for (int k = 0; k < nk; k++) //for each cell at the finest level - k
    {
      Array1dT<Array1dT<SurfaceKernel*> >& currentVal = m_values(j, k); // get the list of kernels for each level
      int nper = 1;
      for (localIndex i = 1; i < nlevels; i++)
        nper *= 2;
      for (localIndex i = 0; i < nlevels; i++, nper /= 2)
      {
        Array2dT<SurfaceKernel*>& curr = vals[i];
        //get the coordinates of the parent cell at level i
        int jj0 = j / nper;
        int jj1 = jj0 + 2 * ioffset;
        int kk0 = k / nper;
        int kk1 = kk0 + 2 * ioffset;
        for (int jj = jj0; jj <= jj1; jj++)
        {
          for (int kk = kk0; kk <= kk1; kk++)
          {
            currentVal[i].push_back(curr(jj, kk));
          }
        }
      }
    }
  }
}

realT FractalSurface::Value(const R1TensorT<2>& position) const
{
  //get the normalized position to evaluate in the hierarchical grid
  R1TensorT<2> xt(position);
  {
    R1TensorT<2> tmp2(m_upper);
    tmp2 -= m_lower;
    for (localIndex i=0; i<2; ++i)
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
  int coords[2];
  for(localIndex i = 0; i < 2; i++)
    coords[i] = (int)(m_values.Dimension(i) * (xt(i) >= 1 ? (1.0-1.0e-10) : xt(i)) );
  const Array1dT<Array1dT<SurfaceKernel*> >& valArray = m_values(coords[0], coords[1]);

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
