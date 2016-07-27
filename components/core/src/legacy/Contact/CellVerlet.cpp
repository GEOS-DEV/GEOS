//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************
 * @file CellVerlet.cpp
 * @date Nov 11, 2011
 * @author Scott Johnson
 *
 * @brief Cell-Verlet spatial sorting class
 ************************************************************/

#include "CellVerlet.h"
#include "Utilities/Utilities.h"

#include "SpatialSorterFactory.h"

namespace SpatialSorting
{

  CellVerlet::CellVerlet() :
      binFct(1.5),
      bins(),
      neighborBins(),
      xmin(std::numeric_limits<realT>::max()),
      xmax(-std::numeric_limits<realT>::max()),
      dx(0.0)
  {
  }

  CellVerlet::~CellVerlet()
  {
  }

  void CellVerlet::Clear()
  {
    bins.clear();
    neighborBins.clear();
  }

  bool CellVerlet::UpdateMinMaxDimension(const rArray1d& radii, const Array1dT<R1Tensor>& x,
                                         const int* const excludeFromContact)
  {
    //initialize the other variables
    realT dx_new = 0.;
    R1Tensor xx1 = static_cast<R1Tensor>( -std::numeric_limits<realT>::max() );
    R1Tensor xx0 = static_cast<R1Tensor>( std::numeric_limits<realT>::max() );

    const localIndex num = radii.size();
    if (excludeFromContact != 0)
    {
      for (localIndex a = 0; a < num; ++a)
      {
        if (excludeFromContact[a] > 0)
          continue;
        dx_new = radii[a] > dx_new ? radii[a] : dx_new;
        xx0.SetMin(x[a]);
        xx1.SetMax(x[a]);
      }
    }
    else
    {
      for (localIndex a = 0; a < num; ++a)
      {
        dx_new = radii[a] > dx_new ? radii[a] : dx_new;
        xx0.SetMin(x[a]);
        xx1.SetMax(x[a]);
      }
    }

    //currently holds max radius; multiply by 2 for the diameter
    dx_new *= 2.;

    bool ret = false;
    if (dx_new > this->dx)
    {
      this->dx = dx_new;
      ret = true;
    }
    for (localIndex i = 0; i < nsdof; ++i)
    {
      if (this->xmin(i) > xx0(i))
      {
        this->xmin(i) = xx0(i);
        ret = true;
      }
      if (this->xmax(i) < xx1(i))
      {
        this->xmax(i) = xx1(i);
        ret = true;
      }
    }
    return ret;
  }

  bool CellVerlet::UpdateMinMaxDimension(const rArray1d& radii, const Array1dT<R1Tensor>& x,
                                         const lSet& toResort)
  {
    //initialize the other variables
    realT dx_new = 0.;
    R1Tensor xx1 = static_cast<R1Tensor>(-std::numeric_limits<realT>::max());
    R1Tensor xx0 = static_cast<R1Tensor>(std::numeric_limits<realT>::max());

    for (lSet::const_iterator iter = toResort.begin(); iter != toResort.end(); ++iter)
    {
      const localIndex a = *iter;
      dx_new = radii[a] > dx_new ? radii[a] : dx_new;
      xx0.SetMin(x[a]);
      xx1.SetMax(x[a]);
    }

    //currently holds max radius; multiply by 2 for the diameter
    dx_new *= 2.;

    bool ret = false;
    if (dx_new > this->dx)
    {
      this->dx = dx_new;
      ret = true;
    }
    for (localIndex i = 0; i < nsdof; ++i)
    {
      if (this->xmin(i) > xx0(i))
      {
        this->xmin(i) = xx0(i);
        ret = true;
      }
      if (this->xmax(i) < xx1(i))
      {
        this->xmax(i) = xx1(i);
        ret = true;
      }
    }
    return ret;
  }


  bool CellVerlet::Update(const rArray1d& radii, const Array1dT<R1Tensor>& x, const lSet& toResort,
                          Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                          const int* const excludeFromSorting)
  {
    if (this->UpdateMinMaxDimension(radii, x, toResort))
      return SortSub(true, radii, x, neighborList, neighborListInverse, excludeFromSorting);

    //things are not so bad ... everything has stayed within the bounds of the assumed domain
    //first, let's remove neighbor list entries involving those that need to be resorted
    SpatialSorterBase::Remove(toResort, neighborList, neighborListInverse); //remove the moved from the neighbor lists

    //second, let's re-insert the entries
    iArray1d index(nsdof);
    for (lSet::const_iterator iter = toResort.begin(); iter != toResort.end();
        ++iter)
    {
      BinIndices(x[*iter], index);
      Remove(*iter, index);
    }
    for (lSet::const_iterator iter = toResort.begin(); iter != toResort.end();
        ++iter)
    {
      BinIndices(x[*iter], index);
      Add(radii, x, *iter, index, neighborList, neighborListInverse, 0);
    }
    return true;
  }

  void CellVerlet::RemoveToCheck(const std::map<localIndex, iArray1d>& toCheckFurther)
  {
    //remove those to check from the current bins
    for (int i = 0; i < this->bins.Dimension(0); ++i)
    {
      for (int j = 0; j < this->bins.Dimension(1); ++j)
      {
        for (int k = 0; k < this->bins.Dimension(2); ++k)
        {
          lArray1d& current = this->bins(i, j, k);
          lArray1d::iterator iter = current.begin();
          while (iter != current.end())
          {
            std::map<localIndex, iArray1d>::const_iterator it = toCheckFurther.find(*iter);
            if (it != toCheckFurther.end())
              iter = current.erase(iter);
            else
              ++iter;
          }
        }
      }
    }
  }

  /**
   * @brief Brief
   * @author Scott Johnson
   * Description
   * @param[in] radii Radii vector
   * @param[in] x Centers vector
   * @param[out] neighborList Neighbor
   * @param[out] neighborListInverse Neighbor inverse
   * @param[in] excludeFromSorting Exclusion list with integer binary values
   */
  bool CellVerlet::Sort(const rArray1d& radii, const Array1dT<R1Tensor>& x,
                        Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                        const int* const excludeFromSorting)
  {
    //update problem dimensions
    const bool reset = this->UpdateMinMaxDimension(radii, x, excludeFromSorting);
    return SortSub(reset, radii, x, neighborList, neighborListInverse, excludeFromSorting);

  } //Sort

  bool CellVerlet::Remove(const localIndex iRemove, const iArray1d& guess)
  {
    //check last known location
    {
      lArray1d& binCurrent = this->bins(guess[0], guess[1], guess[2]);
      lArray1d::iterator it = binCurrent.begin();
      while(it != binCurrent.end())
      {
        if(*it == iRemove)
        {
          binCurrent.erase(it);
          return true;
        }
        ++it;
      }
    }

    //check neighbors
    {
      const Array1dT<lArray1d*>& neighborBinsCurrent = this->neighborBins(guess[0], guess[1],
                                                                          guess[2]);
      for (Array1dT<lArray1d*>::const_iterator i0 = neighborBinsCurrent.begin();
          i0 != neighborBinsCurrent.end(); ++i0)
      {
        lArray1d& binCurrent = *(*i0);
        lArray1d::iterator it = binCurrent.begin();
        while(it != binCurrent.end())
        {
          if(*it == iRemove)
          {
            binCurrent.erase(it);
            return true;
          }
          ++it;
        }
      }
    }
    return false;
  }



  void CellVerlet::Add(const rArray1d& radii, const Array1dT<R1Tensor>& x, const localIndex ixfc0,
                       const iArray1d& index, Array1dT<lArray1d>& neighborList,
                       Array1dT<lSet>& neighborListInverse, const int* const excludeFromSorting)
  {
    //check neighbors
    {
      const Array1dT<lArray1d*>& neighborBinsCurrent = this->neighborBins(index[0], index[1],
                                                                          index[2]);
      for (Array1dT<lArray1d*>::const_iterator i0 = neighborBinsCurrent.begin();
          i0 != neighborBinsCurrent.end(); ++i0)
      {
        const lArray1d& neighbor = *(*i0);
        for (lArray1d::const_iterator i1 = neighbor.begin(); i1 != neighbor.end(); ++i1)
          SpatialSorterBase::AddIfClose(ixfc0, *i1, radii, x, neighborList, neighborListInverse);
      }
    }

    //check home bin and add
    {
      lArray1d& binCurrent = this->bins(index[0], index[1], index[2]);
      for (lArray1d::const_iterator i1 = binCurrent.begin(); i1 != binCurrent.end(); ++i1)
        SpatialSorterBase::AddIfClose(ixfc0, *i1, radii, x, neighborList, neighborListInverse);
      binCurrent.push_back(ixfc0);
    }
  }

  /**
   * @brief Brief
   * @author Scott Johnson
   * Description
   * @param[in] radii Radii vector
   * @param[in] x Centers vector
   * @param[out] neighborList Neighbor
   * @param[out] neighborListInverse Neighbor inverse
   * @param[in] excludeFromSorting Exclusion list with integer binary values
   */
  bool CellVerlet::SortSub(const bool reset, const rArray1d& radii, const Array1dT<R1Tensor>& x,
                           Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                           const int* const excludeFromSorting)
  {
    //---------------------------------------------------
    // DEAL WITH A TOTAL REBUILDING OF THE BIN STRUCTURE
    //---------------------------------------------------
    iArray1d index(nsdof);
    if (reset)
    {
      //std::cout << "CELL VERLET: rebuilding bin structure" << std::endl;

      this->Clear();

      //give it a little buffer, so you're not always updating ...
      this->xmin -= this->dx;
      this->xmax += this->dx;

      //return if the faces are dimensionless ... should this be an exception?
      if (isZero(this->dx))
        return false;

      //allocate the bins
      {
        R1Tensor xx(this->xmax);
        BinIndices(xx, index);
        ++index[0];
        ++index[1];
        ++index[2];
        this->bins.Allocate(index[0], index[1], index[2]);
        this->neighborBins.Allocate(index[0], index[1], index[2]);
      }

      //assign neighbor bin lists
      for (long ix = 0; ix < this->bins.Dimension(0); ++ix)
      {
        for (long iy = 0; iy < this->bins.Dimension(1); ++iy)
        {
          for (long iz = 0; iz < this->bins.Dimension(2); ++iz)
          {
            for (int i = -1; i < 2; i++)
            {
              long ii = ix + i;
              if (ii >= 0 && ii < this->bins.Dimension(0))
              {
                for (int j = -1; j < 2; j++)
                {
                  long jj = iy + j;
                  if (jj >= 0 && jj < this->bins.Dimension(1))
                  {
                    for (int k = -1; k < 2; k++)
                    {
                      long kk = iz + k;
                      if (kk >= 0 && kk < this->bins.Dimension(2))
                      {
                        //std::cout << "   trying " << ix << " " << iy << " " << iz << " <-- " << ii << " " << jj << " " << kk;
                        const bool notok = (i == 0) && (j == 0) && (k == 0);
                        if (notok)
                        {
                          //std::cout << " NOT OK " << i << " " << j << " " << k << "\n";
                          continue;
                        }
                        lArray1d& neighbor = this->bins(ii, jj, kk);
                        this->neighborBins(ix, iy, iz).push_back(&neighbor);
                        //std::cout << " ok\n";
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else
    {
      //just clear the entries in the existing bin structures
      //std::cout << "CELL VERLET: clearing bin structure's lists" << std::endl;
      for (long ix = 0; ix < this->bins.Dimension(0); ++ix)
      {
        for (long iy = 0; iy < this->bins.Dimension(1); ++iy)
        {
          for (long iz = 0; iz < this->bins.Dimension(2); ++iz)
          {
            this->bins(ix, iy, iz).clear();
            this->neighborBins(ix, iy, iz).clear();
          }
        }
      }
    }

    //--------------------------------------------------------
    // DEAL WITH REFILLING THE NEIGHBOR LIST AND ITS INVERSE
    //--------------------------------------------------------

    //clear neighbor lists
    //std::cout << "CELL VERLET: rebuilding neighbor lists" << std::endl;
    for (localIndex ixfc0 = 0; ixfc0 < neighborList.size(); ++ixfc0)
    {
      neighborList[ixfc0].clear();
      neighborListInverse[ixfc0].clear();
    }

    //now, bin the faces and simultaneously fill the neighbor list
    const localIndex num = radii.size();
    for (localIndex ixfc0 = 0; ixfc0 < num; ++ixfc0)
    {
      if (excludeFromSorting != 0 && excludeFromSorting[ixfc0] > 0)
        continue;
      BinIndices(x[ixfc0], index);
      //std::cout << ixfc0 << ": index = " << index[0] << " " << index[1] << " " << index[2] << "\n";
      Add(radii, x, ixfc0, index, neighborList, neighborListInverse, excludeFromSorting);
    }
    return true;
  }
}

/// Register spatial sorter in the spatial sorter factory
REGISTER_SPATIALSORTER( CellVerlet )
