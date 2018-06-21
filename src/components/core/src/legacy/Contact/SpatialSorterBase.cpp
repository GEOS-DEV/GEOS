// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
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
{}

SpatialSorterBase::~SpatialSorterBase()
{}

void SpatialSorterBase::Clear()
{}

void SpatialSorterBase::Remove(const lSet& toRemove,
                               array<lArray1d>& neighborList,
                               array<lSet>& neighborListInverse)
{
  for(lSet::const_iterator iter = toRemove.begin() ; iter != toRemove.end() ; ++iter)
  {
    lArray1d& current = neighborList[*iter];
    for(lArray1d::iterator it = current.begin() ; it != current.end() ; ++it)
    {
      lSet& tmp = neighborListInverse[*it];
      lSet::iterator itset = tmp.find(*iter);
      if(itset != tmp.end())
        tmp.erase(itset);
    }
    current.clear();
  }
}

void SpatialSorterBase::RemoveDuplicates(array<lArray1d>& neighborList)
{
  lArray1d::iterator it;
  for(array<lArray1d>::iterator iter = neighborList.begin() ; iter != neighborList.end() ; ++iter)
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
                               array<lArray1d>& neighborList,
                               array<lSet>& neighborListInverse)
{
  if(i >= neighborList.size())
    throw GPException("SpatialSorter::insert out of range");

  //update neighborList
  {
    for(localIndex ii = 0 ; ii < i ; ++ii)
    {
      lArray1d& tmp = neighborList[ii];
      for(lArray1d::iterator it = tmp.begin() ; it != tmp.end() ; ++it)
        if(*it >= i)
          (*it) += count;
    }
    for(localIndex ii = i ; ii < neighborList.size() ; ++ii)
    {
      lArray1d& tmp = neighborList[ii];
      for(lArray1d::iterator it = tmp.begin() ; it != tmp.end() ; ++it)
        (*it) += count;
    }
    for(localIndex ii = 0 ; ii < count ; ii++)
      neighborList.insert(neighborList.begin() + i, lArray1d());
  }

  //update inverse
  {
    for(localIndex ii = i+1 ; ii < neighborListInverse.size() ; ++ii)
    {
      const lSet& tmp = neighborListInverse[ii];
      lSet newList;
      for(lSet::const_iterator it = tmp.begin() ; it != tmp.end() ; ++it)
        newList.insert(*it >= i ? *it + count : *it);
      neighborListInverse[ii] = newList;
    }
    for(localIndex ii = 0 ; ii < count ; ii++)
      neighborListInverse.insert(neighborListInverse.begin() + i, lSet());
  }
}

void SpatialSorterBase::Add(const localIndex kf0, const localIndex kf1,
                            array<lArray1d>& neighborList,
                            array<lSet>& neighborListInverse)
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
                                   const array<real64>& radii, const array<R1Tensor>& centers,
                                   array<lArray1d>& neighborList,
                                   array<lSet>& neighborListInverse)
{
  if(Close(kf0, kf1, radii, centers))
    Add(kf0, kf1, neighborList, neighborListInverse);
}
}
