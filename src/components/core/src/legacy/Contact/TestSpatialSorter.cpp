#include "SurfaceGeneration/StatisticalDistributionBaseT.h"

#include "SpatialSorter.h"
#include "CellVerlet.h"
#include "CGRID.h"
#include "DESS.h"
#include "NBS.h"
#include "NSquared.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

template<int T>
localIndex Fill(const localIndex num, rArray1d& radii, Array1dT<R1Tensor>& x)
{
  StatisticalDistributionBaseT::InitializeRandom();

  radii.reserve(num);
  x.reserve(num);
  const realT radius = 0.6 * pow(num, -1.0 / T);
  for (localIndex i = 0; i < num; i++)
  {
    const realT r = radius * (1.0 + StatisticalDistributionBaseT::UniformSample(-0.05, 0.05));//radius +/- 5%
    R1Tensor xx;
    for (localIndex j = 0; j < T; j++)
      xx(j) = StatisticalDistributionBaseT::UniformSample(-0.5, 0.5);
    radii.push_back(r);
    x.push_back(xx);

    std::cout << i << ": " << r << " ( " << xx(0) << " , " << xx(1) << " , " << xx(2) << " )\n";
  }
  return num;
}

template<int T>
localIndex FillRegular(const localIndex numPerDim, rArray1d& radii, Array1dT<R1Tensor>& x)
{
  StatisticalDistributionBaseT::InitializeRandom();

  localIndex num = numPerDim;
  for(localIndex i = 1; i < T; i++)
    num *= numPerDim;

  radii.reserve(num);
  x.reserve(num);
  const realT radius = 1.0;
  for (localIndex i = 0; i < numPerDim; i++)
  {
    R1Tensor xx = 0.0;
    xx(0) = i * radius * 1.8;
    if(T > 1)
    {
      for(localIndex j = 0; j < numPerDim; j++)
      {
        xx(1) = j * radius * 1.8;
        if(T > 2)
        {
          for(localIndex k = 0; k < numPerDim; k++)
          {
            xx(2) = k * radius * 1.8;
            R1Tensor tmp(xx);
            x.push_back(tmp);
            radii.push_back(radius);
            std::cout << "[" << i << ", " << j << ", " << k << "] : r=" << radius << " x=( " << xx(0) << " , " << xx(1) << " , " << xx(2) << " )\n";
          }
        }
        else
        {
          R1Tensor tmp(xx);
          x.push_back(tmp);
          radii.push_back(radius);
          std::cout << "[" << i << ", " << j << "] : r=" << radius << " x=( " << xx(0) << " , " << xx(1) << " , " << xx(2) << " )\n";
        }
      }
    }
    else
    {
      R1Tensor tmp(xx);
      x.push_back(tmp);
      radii.push_back(radius);
      std::cout << i << ": r=" << radius << " x=( " << xx(0) << " , " << xx(1) << " , " << xx(2) << " )\n";
    }
  }
  return num;
}

bool Compare(const Array1dT<lArray1d>& truth, const Array1dT<lArray1d>& current)
{
  if (truth.size() != current.size())
    return false;

  for (localIndex i = 0; i < truth.size(); i++)
  {
    const lArray1d& truthsub = truth[i];
    const lArray1d& currentsub = current[i];

    std::cout << i << ":";

    for (localIndex j = 0; j < truthsub.size(); j++)
    {
      std::cout << " " << truthsub[j];
      bool found = false;
      for (localIndex k = 0; k < currentsub.size(); k++)
      {
        if (truthsub[j] == currentsub[k])
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        std::cout << " <--NOT FOUND (current sub size is " << currentsub.size() << " :";
        for(localIndex k = 0; k < currentsub.size(); k++)
          std::cout << " " << currentsub[k];
        std::cout << ") \n";
        return false;
      }
    }
    if(truthsub.size() < currentsub.size())
    {
      std::cout << " (";
      for (localIndex j = 0; j < currentsub.size(); j++)
      {
        std::cout << currentsub[j] << ", ";
      }
      std::cout << ")";
//      std::cout << "new neighbor list too long. extra entries: \n";
//      lArray1d tmp(currentsub);
//      lArray1d::iterator it;
//      for (localIndex j = 0; j < truthsub.size(); j++)
//      {
//        it = tmp.begin();
//        while(it != tmp.end())
//        {
//          if(*it == truthsub[j])
//          {
//            it = tmp.erase(it);
//            break;
//          }
//          else
//            ++it;
//        }
//      }
//      for(localIndex j = 0; j < tmp.size(); j++)
//        std::cout << tmp[j] << ", ";
    }
    std::cout << std::endl;
  }
  return true;
}

void Sort(const localIndex type, const rArray1d& radii, const Array1dT<R1Tensor>& x,
          const lSet& toResort,
          Array1dT<lArray1d>& current, Array1dT<lSet>& currentInverse)
{
  static SpatialSorting::CellVerlet cellverlet;
  static SpatialSorting::CGRID cgrid;
  static SpatialSorting::DESS dess;
  static SpatialSorting::NBS nbs;
  static SpatialSorting::NSquared nsquared;
  const bool fullSort = toResort.size() == 0;

  switch (type)
  {
    case 0:
      if(fullSort)
        nsquared.Sort(radii, x, current, currentInverse);
      else
        nsquared.Update(radii, x, toResort, current, currentInverse);
      break;
    case 1:
      if(fullSort)
        cellverlet.Sort(radii, x, current, currentInverse);
      else
        cellverlet.Update(radii, x, toResort, current, currentInverse);
      break;
    case 2:
      if(fullSort)
        dess.Sort(radii, x, current, currentInverse);
      else
        dess.Update(radii, x, toResort, current, currentInverse);
      break;
    case 3:
      if(fullSort)
        nbs.Sort(radii, x, current, currentInverse);
      else
        nbs.Update(radii, x, toResort, current, currentInverse);
      break;
    case 4:
      if(fullSort)
        cgrid.Sort(radii, x, current, currentInverse);
      else
        cgrid.Update(radii, x, toResort, current, currentInverse);
      break;
    default:
      break;
  }
}

int main(int* argc, char** argv)
{
  const localIndex typeLast = 2;

  std::cout << "+++++++++++++++++++++++++++++++++++++++++\n";
  std::cout << "+       TESTING SPATIAL SORTING         +\n";
  std::cout << "+++++++++++++++++++++++++++++++++++++++++\n";

  //setup the problem with a random distribution
  rArray1d radii;
  Array1dT<R1Tensor> x;
  std::cout << "\n-> FILLING -- \n";

  const localIndex num = Fill<3>(80, radii, x);
  std::cout << "   filled " << num << "\n";
  //exit(1);

  //check each sorter
  Array1dT<lArray1d> current;
  Array1dT<lSet> currentInverse;
  current.resize(num);
  currentInverse.resize(num);
  lSet toResort, toResortUpdate;

  //do initial sort and check
  {
    //get the truth
    std::cout << "\n-> GETTING THE TRUTH\n";
    Array1dT<lArray1d> truth;
    {
      Array1dT<lSet> truthInverse;
      truth.resize(num);
      truthInverse.resize(num);
      Sort(0, radii, x, toResort, truth, truthInverse);

      localIndex totalEntries = 0;
      for(localIndex i = 0; i < truth.size(); i++)
        totalEntries += truth[i].size();
      std::cout << "   total number of entries: " << totalEntries << "\n";
      //return 1;
    }

    //swap positions 1 & 2
    {
      const realT r0 = radii[0];
      radii[0] = radii[1];
      radii[1] = r0;

      R1Tensor x0 = x[0];
      x[0] = x[1];
      x[1] = x0;
    }

    //get the truth
    std::cout << "\n-> GETTING THE TRUTH FOR UPDATE\n";
    Array1dT<lArray1d> truthUpdate;
    {
      Array1dT<lSet> truthInverseUpdate;
      truthUpdate.resize(num);
      truthInverseUpdate.resize(num);
      Sort(0, radii, x, toResort, truthUpdate, truthInverseUpdate);

      localIndex totalEntries = 0;
      for(localIndex i = 0; i < truthUpdate.size(); i++)
        totalEntries += truthUpdate[i].size();
      std::cout << "   total number of entries: " << totalEntries << "\n";
      //return 1;

      toResortUpdate.insert(0);
      toResortUpdate.insert(1);
    }

    for (localIndex type = 0; type <= typeLast; type++)
    {
      std::cout << "+++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "+++   " << type << "   +++" << std::endl;
      std::cout << "+++++++++++++++++++++++++++++++++++++++" << std::endl;
      //swap positions 1 & 2 ... again
      {
        const realT r0 = radii[0];
        radii[0] = radii[1];
        radii[1] = r0;

        R1Tensor x0 = x[0];
        x[0] = x[1];
        x[1] = x0;
      }
      Sort(type, radii, x, toResort, current, currentInverse);
      if (Compare(truth, current))
      {
        std::cout << "SUCCESS: " << type << " ";
      }
      else
      {
        std::cout << "***** failure *****: " << type << " ";
      }
      localIndex totalEntries = 0;
      for(localIndex i = 0; i < current.size(); i++)
        totalEntries += current[i].size();
      std::cout << "(" << totalEntries << ") ";

      switch (type)
      {
        case 0:
          std::cout << "N^2\n";
          break;
        case 1:
          std::cout << "CellVerlet\n";
          break;
        case 2:
          std::cout << "DESS\n";
          break;
        case 3:
          std::cout << "NBS\n";
          break;
        case 4:
          std::cout << "CGRID\n";
          break;
        default:
          std::cout << "--unknown--\n";
          break;
      }

      //swap positions 1 & 2 ... again
      {
        const realT r0 = radii[0];
        radii[0] = radii[1];
        radii[1] = r0;

        R1Tensor x0 = x[0];
        x[0] = x[1];
        x[1] = x0;
      }
      Sort(type, radii, x, toResortUpdate, current, currentInverse);
      if (Compare(truthUpdate, current))
      {
        std::cout << "SUCCESS (UPDATE): " << type << " ";
      }
      else
      {
        std::cout << "***** failure (update) *****: " << type << " ";
      }
      totalEntries = 0;
      for(localIndex i = 0; i < current.size(); i++)
        totalEntries += current[i].size();
      std::cout << "(" << totalEntries << ") ";

      switch (type)
      {
        case 0:
          std::cout << "N^2\n";
          break;
        case 1:
          std::cout << "CellVerlet\n";
          break;
        case 2:
          std::cout << "DESS\n";
          break;
        case 3:
          std::cout << "NBS\n";
          break;
        case 4:
          std::cout << "CGRID\n";
          break;
        default:
          std::cout << "--unknown--\n";
          break;
      }
      std::cout << "+++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "+++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << std::endl;
    }
  }
  return 0;
}
