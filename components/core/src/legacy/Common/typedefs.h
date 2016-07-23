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
/**
 * @file typedefs.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include "intrinsic_typedefs.h"
#include "TensorT/TensorT.h"
#include "ArrayT/Array1dT.h"
#include "ArrayT/Array2dT.h"
#include <set>

/// set number of degrees of freedom
const unsigned int nsdof = 3;

/// define the templated R2TensorT<> with nsdof as its template argument
typedef R2TensorT<nsdof>    R2Tensor;

/// define the templated R2SymTensorT<> with nsdof as its template argument
typedef R2SymTensorT<nsdof> R2SymTensor;

/// define the templated R1TensorT<> with nsdof as its template argument
typedef R1TensorT<nsdof>    R1Tensor;


typedef Array1dT<int> iArray1d;
typedef Array1dT<realT> rArray1d;
typedef Array1dT<std::string> sArray1d;
typedef Array1dT<iArray1d*> iArrayPtrs;
typedef Array1dT<rArray1d*> rArrayPtrs;

typedef Array2dT<int> iArray2d;
typedef Array2dT<localIndex> lArray2d;
typedef Array2dT<globalIndex> gArray2d;
typedef Array2dT<realT> rArray2d;
typedef Array2dT<std::pair<int,localIndex> > pArray2d;

typedef Array1dT<localIndex> lArray1d;
typedef Array1dT<globalIndex> gArray1d;
typedef Array1dT<std::pair<int,localIndex> > pArray1d;


typedef std::set<int> iSet;
typedef std::set<localIndex> lSet;
typedef std::set<globalIndex> gSet;
typedef std::set<std::pair<int,localIndex> > pSet;
typedef std::set<std::string> sSet;

template<typename T> struct type_name
{
    static const char* name();// { static_assert(false, "You are missing a DECL_TYPE_NAME"); }
};

template<> struct type_name<int>                { static const char* name() {return "int";} };
template<> struct type_name<unsigned int>       { static const char* name() {return "unsigned int";} };
template<> struct type_name<long long>          { static const char* name() {return "long long";} };
template<> struct type_name<unsigned long long> { static const char* name() {return "unsigned long long";} };
template<> struct type_name<float>              { static const char* name() {return "float";} };
template<> struct type_name<double>             { static const char* name() {return "double";} };

/// stores pointer to elemenRegion and localIndex of element
class ElementRegionT;
typedef std::pair< ElementRegionT*, localIndex > ElementIdPair;

template <class Type>
struct locallyIndexedType
{
  localIndex m_index;
  Type m_x;

public:
  locallyIndexedType(localIndex indx, Type x):m_index(indx),m_x(x){ /** empty **/}
  virtual localIndex GetIndex(){return m_index;}
  virtual Type GetValue(){return m_x;}
};

typedef locallyIndexedType<realT> locallyIndexedReal;


#endif /* TYPEDEFS_H_ */
