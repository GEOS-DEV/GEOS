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
#ifndef STRUCTURED_GRID_UTILITIES_H
#define STRUCTURED_GRID_UTILITIES_H

/**
 * @file StructuredGridUtilities.h
 * @author white230
 */
 
//#include "legacy/Common/Common.h"
#include <cassert>

namespace StructuredGrid
{
  /*!
   *  Given n, compute n^d, where d is the spatial dimension.
   */
  
  template <int dim>
  int dimpower(int n);

  template <> inline int dimpower<1>(int n) { return n; }
  template <> inline int dimpower<2>(int n) { return n*n; }
  template <> inline int dimpower<3>(int n) { return n*n*n; }

  /*!
   * Given a lexographical index N, map back to the original
   * i,j,k indices of the point. The first variation here assumes
   * a uniform number of points nnx in all coordinate directions.
   */

  template <int dim>
  void map_index(const unsigned index,
                 const unsigned nnx,
                 std::vector<unsigned> &indices);

  template <>
  inline
  void map_index<1>(const unsigned index,
                    const unsigned nnx,
                    std::vector<unsigned> &indices)
  {
    assert(index < nnx);
    indices[0] = index;
  }

  template <>
  inline
  void map_index<2>(const unsigned index,
                    const unsigned nnx,
                    std::vector<unsigned> &indices)
  {
    assert(index < nnx*nnx);
    indices[0] = index % nnx;
    indices[1] = index / nnx;
  }

  template <>
  inline
  void map_index<3>(const unsigned index,
                    const unsigned nnx,
                    std::vector<unsigned> &indices)
  {
    assert(index < nnx*nnx*nnx);
    indices[0] = index % nnx;
    indices[1] = (index / nnx) % nnx;
    indices[2] = index / (nnx*nnx);
  }

} // end namespace

#endif

