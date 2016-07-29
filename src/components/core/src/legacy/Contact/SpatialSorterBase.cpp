/*
 * SpatialSorter.cpp
 *
 *  Created on: Nov 11, 2011
 *      Author: scottjohnson
 */

#include "SpatialSorterBase.h"

namespace SpatialSorting
{

  SpatialSorterBase::SpatialSorterBase()
  {
  }

  SpatialSorterBase::~SpatialSorterBase()
  {
  }

  void SpatialSorterBase::Clear()
  {
  }

  void SpatialSorterBase::Remove(const lSet& toRemove,
                             Array1dT<lArray1d>& neighborList,
                             Array1dT<lSet>& neighborListInverse)
  {
    for(lSet::const_iterator iter = toRemove.begin(); iter != toRemove.end(); ++iter)
    {
      lArray1d& current = neighborList[*iter];
      for(lArray1d::iterator it = current.begin(); it != current.end(); ++it)
      {
        lSet& tmp = neighborListInverse[*it];
        lSet::iterator itset = tmp.find(*iter);
        if(itset != tmp.end())
          tmp.erase(itset);
      }
      current.clear();
    }
  }

  void SpatialSorterBase::RemoveDuplicates(Array1dT<lArray1d>& neighborList)
  {
    lArray1d::iterator it;
    for(Array1dT<lArray1d>::iterator iter = neighborList.begin(); iter != neighborList.end(); ++iter)
    {
      lArray1d& curr = *iter;
      if(curr.size() < 2)
        continue;
      std::sort(curr.begin(), curr.end());
      it = curr.begin();
      localIndex last = *it;
      ++it;
      while(it != curr.end())
      {
        if(*it == last)
        {
          it = curr.erase(it);
        }
        else
        {
          last = *it;
          ++it;
        }
      }
    }
  }

  void SpatialSorterBase::insert(const localIndex i, const localIndex count,
                             Array1dT<lArray1d>& neighborList,
                             Array1dT<lSet>& neighborListInverse)
  {
    if(i >= neighborList.size())
      throw GPException("SpatialSorter::insert out of range");

    //update neighborList
    {
      for(localIndex ii = 0; ii < i; ++ii)
      {
        lArray1d& tmp = neighborList[ii];
        for(lArray1d::iterator it = tmp.begin(); it != tmp.end(); ++it)
            if(*it >= i)
              (*it) += count;
      }
      for(localIndex ii = i; ii < neighborList.size(); ++ii)
      {
        lArray1d& tmp = neighborList[ii];
        for(lArray1d::iterator it = tmp.begin(); it != tmp.end(); ++it)
          (*it) += count;
      }
      for(localIndex ii = 0; ii < count; ii++)
        neighborList.insert(neighborList.begin() + i, lArray1d());
    }

    //update inverse
    {
      for(localIndex ii = i+1; ii < neighborListInverse.size(); ++ii)
      {
        const lSet& tmp = neighborListInverse[ii];
        lSet newList;
        for(lSet::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
          newList.insert(*it >= i ? *it + count : *it);
        neighborListInverse[ii] = newList;
      }
      for(localIndex ii = 0; ii < count; ii++)
        neighborListInverse.insert(neighborListInverse.begin() + i, lSet());
    }
  }

  void SpatialSorterBase::Add(const localIndex kf0, const localIndex kf1,
                          Array1dT<lArray1d>& neighborList,
                          Array1dT<lSet>& neighborListInverse)
  {
    if(kf0 < kf1)
    {
      neighborList[kf0].push_back(kf1);
      neighborListInverse[kf1].insert(kf0);
      //std::cout << "add " << kf0 << " " << kf1 << "\n";
    }
    else
    {
      neighborList[kf1].push_back(kf0);
      neighborListInverse[kf0].insert(kf1);
      //std::cout << "add " << kf1 << " " << kf0 << "\n";
    }
  }

  void SpatialSorterBase::AddIfClose(const localIndex kf0, const localIndex kf1,
                                 const rArray1d& radii, const Array1dT<R1Tensor>& centers,
                                 Array1dT<lArray1d>& neighborList,
                                 Array1dT<lSet>& neighborListInverse)
  {
    if(Close(kf0, kf1, radii, centers))
      Add(kf0, kf1, neighborList, neighborListInverse);
  }
}
