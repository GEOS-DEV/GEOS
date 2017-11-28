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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************
* @file NSquared.cpp
* @date Nov 11, 2011
* @author Scott Johnson
*
* @brief Order N^2 spatial sorting class
************************************************************/

#include "NSquared.h"

#include "SpatialSorterFactory.h"

namespace SpatialSorting
{

NSquared::NSquared()
{
  // TODO Auto-generated constructor stub

}

NSquared::~NSquared()
{
  // TODO Auto-generated destructor stub
}

bool NSquared::Update(const array<real64>& radii,
                      const array<R1Tensor>& x,
                      const lSet& toResort,
                      array<lArray1d>& neighborList,
                      array<lSet>& neighborListInverse,
                      const int* const excludeFromSorting)
{
  return Sort(radii, x, neighborList, neighborListInverse, excludeFromSorting);
}

/**
 * @brief Brief
 * @author Scott Johnson
 * Description
 * @param[in] radii Radii vector
 * @param[in] centers Centers vector
 * @param[out] neighborList Neighbor
 */
bool NSquared::Sort(const array<real64>& radii,
                    const array<R1Tensor>& centers,
                    array<lArray1d>& neighborList,
                    array<lSet>& neighborListInverse,
                    const int* const excludeFromContact)
{
  localIndex num = radii.size();

  //here, we at least limit by bounding sphere interaction (2)
  R1Tensor dd(0.);
  for (localIndex kf0 = 0 ; kf0 < num ; ++kf0)
  {
    neighborList[kf0].clear();
    neighborListInverse[kf0].clear();
  }

  //iterate through first
  if(excludeFromContact != 0)
  {
    for (localIndex kf0 = 0 ; kf0 < num ; ++kf0)
    {
      if(excludeFromContact[kf0]>0)
        continue;

      //iterate through second
      for (localIndex kf1 = kf0 + 1 ; kf1 < num ; ++kf1)
      {
        if(excludeFromContact[kf1]>0)
          continue;
        SpatialSorterBase::AddIfClose(kf0, kf1, radii, centers, neighborList, neighborListInverse);
      }
    }
  }
  else
  {
    for (localIndex kf0 = 0 ; kf0 < num ; ++kf0)
      for (localIndex kf1 = kf0 + 1 ; kf1 < num ; ++kf1)
        SpatialSorterBase::AddIfClose(kf0, kf1, radii, centers, neighborList, neighborListInverse);
  }
  return true;
}  //Sort
}

/// Register spatial sorter in the spatial sorter factory
REGISTER_SPATIALSORTER( NSquared )
