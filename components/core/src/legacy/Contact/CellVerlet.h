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
 * @file CellVerlet.h
 * @date Nov 11, 2011
 * @author Scott Johnson
 *
 * @brief Cell-Verlet spatial sorting class
 ************************************************************/

#ifndef CELLVERLET_H_
#define CELLVERLET_H_

#include "SpatialSorterBase.h"

namespace SpatialSorting
{

  class CellVerlet: public SpatialSorterBase
  {

  public:
    CellVerlet();
    virtual ~CellVerlet();

    static std::string SpatialSorterName() { return "CellVerlet"; }

    virtual bool Sort(const rArray1d& radii, const Array1dT<R1Tensor>& x,
                      Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                      const int* const excludeFromSorting = 0);

    virtual bool Update(const rArray1d& radii, const Array1dT<R1Tensor>& x, const lSet& toResort,
                        Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                        const int* const excludeFromSorting = 0);

    virtual void Clear();

  private:

    void RemoveToCheck(const std::map<localIndex, iArray1d>& toCheckFurther);
    void Add(const rArray1d& radii, const Array1dT<R1Tensor>& x, const localIndex ixfc0, const iArray1d& index,
             Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
             const int* const excludeFromSorting);

    bool UpdateMinMaxDimension(const rArray1d& radii,
                               const Array1dT<R1Tensor>& x,
                               const int* const excludeFromContact);

    bool UpdateMinMaxDimension(const rArray1d& radii,
                               const Array1dT<R1Tensor>& x,
                               const lSet& toResort);

    bool SortSub(const bool reset, const rArray1d& radii, const Array1dT<R1Tensor>& x,
                 Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                 const int* const excludeFromSorting);

    bool Remove(const localIndex, const iArray1d& guess);

    realT binFct;
    Array3dT< lArray1d > bins;
    Array3dT< Array1dT<lArray1d*> > neighborBins;
    R1Tensor xmin, xmax;
    realT dx;

    inline void BinIndices(const R1Tensor& x, iArray1d& bin)
    {
      R1Tensor xx(x);
      xx -= this->xmin; //adjust by the minimum
      xx *= 1.0 / this->dx;
      for (unsigned int i = 0; i < nsdof; i++)
        bin[i] = static_cast<int>(xx(i));
    };
  };

}

#endif /* CELLVERLET_H_ */
