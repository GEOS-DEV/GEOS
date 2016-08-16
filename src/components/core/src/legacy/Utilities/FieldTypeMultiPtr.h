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
 * @file FieldTypeMultiPtr.h
 * @author walsh24
 * @date Jul 7, 2011
 */

#ifndef FIELDTYPEMULTIPTR_H_
#define FIELDTYPEMULTIPTR_H_

#include "Common/Common.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

/// FieldTypeMultiPtr
/// @author walsh24
/// Single class that can contain a pointer to any of the field array types.
class FieldTypeMultiPtr{

  public:

    // Constructors
    FieldTypeMultiPtr():m_fieldType(FieldInfo::realField),m_sArrayPtr(){};
    FieldTypeMultiPtr(Array1dT<int>*   iArrayPtr):m_fieldType(FieldInfo::integerField),m_iArrayPtr(iArrayPtr){};
    FieldTypeMultiPtr(lArray1d*        lArrayPtr):m_fieldType(FieldInfo::localIndexField),m_lArrayPtr(lArrayPtr){};
    FieldTypeMultiPtr(gArray1d*        gArrayPtr):m_fieldType(FieldInfo::globalIndexField),m_gArrayPtr(gArrayPtr){};
    FieldTypeMultiPtr(Array1dT<realT>*        sArrayPtr):m_fieldType(FieldInfo::realField),m_sArrayPtr(sArrayPtr){};
    FieldTypeMultiPtr(Array1dT<R1Tensor>*     vArrayPtr):m_fieldType(FieldInfo::R1TensorField),m_vArrayPtr(vArrayPtr){};
    FieldTypeMultiPtr(Array1dT<R2Tensor>*     tArrayPtr):m_fieldType(FieldInfo::R2TensorField),m_tArrayPtr(tArrayPtr){};
    FieldTypeMultiPtr(Array1dT<R2SymTensor>* stArrayPtr):m_fieldType(FieldInfo::R2SymTensorField),m_stArrayPtr(stArrayPtr){};
    FieldTypeMultiPtr(int* integerPtr):m_fieldType(FieldInfo::integerParameter), m_integerPtr(integerPtr){};
    FieldTypeMultiPtr(realT* scalarPtr): m_fieldType(FieldInfo::realParameter), m_scalarPtr(scalarPtr){};

    // Assignment operators
    FieldTypeMultiPtr& operator=(const FieldTypeMultiPtr& ptr){
      if(this != &ptr){  m_fieldType = ptr.m_fieldType; m_iArrayPtr = ptr.m_iArrayPtr; /*copies union*/};
      return *this;
    };

    FieldTypeMultiPtr& operator=(Array1dT<int>*         const ptr){ m_fieldType = FieldInfo::integerField;  m_iArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(lArray1d*              const ptr){ m_fieldType = FieldInfo::localIndexField;  m_lArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(gArray1d*              const ptr){ m_fieldType = FieldInfo::globalIndexField;  m_gArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(Array1dT<realT>*       const ptr){ m_fieldType = FieldInfo::realField;     m_sArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(Array1dT<R1Tensor>*    const ptr){ m_fieldType = FieldInfo::R1TensorField; m_vArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(Array1dT<R2Tensor>*    const ptr){ m_fieldType = FieldInfo::R2TensorField; m_tArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(Array1dT<R2SymTensor>* const ptr){ m_fieldType = FieldInfo::R2SymTensorField; m_stArrayPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(int* ptr){ m_fieldType = FieldInfo::integerParameter; m_integerPtr = ptr; return *this;};
    FieldTypeMultiPtr& operator=(realT* ptr){ m_fieldType = FieldInfo::realParameter; m_scalarPtr = ptr; return *this;};

    /// GetPtr(ptr) set ptr to the field pointer if pointer type corresponds to m_fieldType, otherwise set ptr to 0;
    void GetPtr(Array1dT<int>*&          iArrayPtr){iArrayPtr  = (m_fieldType == FieldInfo::integerField)? m_iArrayPtr: 0;};
    void GetPtr(lArray1d*&               lArrayPtr){lArrayPtr  = (m_fieldType == FieldInfo::localIndexField)? m_lArrayPtr: 0;};
    void GetPtr(gArray1d*&               gArrayPtr){gArrayPtr  = (m_fieldType == FieldInfo::globalIndexField)? m_gArrayPtr: 0;};
    void GetPtr(Array1dT<realT>*&        sArrayPtr){sArrayPtr  = (m_fieldType == FieldInfo::realField)? m_sArrayPtr: 0;};
    void GetPtr(Array1dT<R1Tensor>*&     vArrayPtr){vArrayPtr  = (m_fieldType == FieldInfo::R1TensorField)? m_vArrayPtr: 0;};
    void GetPtr(Array1dT<R2Tensor>*&     tArrayPtr){tArrayPtr  = (m_fieldType == FieldInfo::R2TensorField)? m_tArrayPtr: 0;};
    void GetPtr(Array1dT<R2SymTensor>*& stArrayPtr){stArrayPtr = (m_fieldType == FieldInfo::R2SymTensorField)? m_stArrayPtr: 0;};
    void GetPtr(int*& integerPtr){integerPtr = (m_fieldType == FieldInfo::integerParameter)? m_integerPtr: 0;};
    void GetPtr(realT*& scalarPtr){scalarPtr = (m_fieldType == FieldInfo::realParameter)? m_scalarPtr: 0;};


    /// Returns the field value at index i, with the specified component.
    /// integer and unsigned fields are converted to real
    realT GetValue(localIndex i,int component){
      switch(m_fieldType){
          case FieldInfo::integerField:     return realT( (*m_iArrayPtr)[i]);
          case FieldInfo::localIndexField:  return realT( (*m_lArrayPtr)[i]);
          case FieldInfo::globalIndexField: return realT( (*m_gArrayPtr)[i]);
          case FieldInfo::realField:        return (*m_sArrayPtr)[i];
          case FieldInfo::R1TensorField:    return (*m_vArrayPtr)[i](component);
          case FieldInfo::R2TensorField:    return (*m_tArrayPtr)[i].Data()[component];
          case FieldInfo::R2SymTensorField: return (*m_stArrayPtr)[i].Data()[component];
          case FieldInfo::integerParameter:     return realT( *m_integerPtr);
          case FieldInfo::realParameter:        return (*m_scalarPtr);
          case FieldInfo::numFieldTypes:
          default:
            throw GPException("Error FieldTypeMultiPtr::GetValue: Unrecognized field type.");

        }
        throw GPException("Error FieldTypeMultiPtr::GetValue: Unrecognized field type.");
        return 0; /** should not get here */
    };

    /// Sets the field value at indx i, with the specified component.
    void SetValue(localIndex i,int component, realT value){
      switch(m_fieldType){
          case FieldInfo::realField:         (*m_sArrayPtr)[i]= value; break;
          case FieldInfo::R1TensorField:     (*m_vArrayPtr)[i](component) = value; break;
          case FieldInfo::R2TensorField:     (*m_tArrayPtr)[i].Data()[component] = value; break;
          case FieldInfo::R2SymTensorField:  (*m_stArrayPtr)[i].Data()[component] = value; break;
          case FieldInfo::realParameter:  *m_scalarPtr = value;break;
          case FieldInfo::integerField:
          case FieldInfo::integerParameter:
          case FieldInfo::numFieldTypes:
            throw GPException("Error FieldTypeMultiPtr::SetValue: Attempting to set integer field with real value.");
            break;
          case FieldInfo::localIndexField:
              throw GPException("Error FieldTypeMultiPtr::SetValue: Attempting to set local index field with real value.");
              break;
          case FieldInfo::globalIndexField:
              throw GPException("Error FieldTypeMultiPtr::SetValue: Attempting to set global index field with real value.");
              break;
          default:
            throw GPException("Error FieldTypeMultiPtr::SetValue: Unrecognized field type.");
        }
    };


    /// Sets the pointer to a new field with type "fieldType" and name "fieldName".
    void SetFieldPtr(ObjectDataStructureBaseT& object, FieldType fieldType, std::string fieldName){
        switch(fieldType){
          case FieldInfo::integerField:     *this = &(object.GetFieldData<int>( fieldName) );  break;
          case FieldInfo::localIndexField:  *this = &(object.GetFieldData<localIndex>( fieldName) );  break;
          case FieldInfo::globalIndexField: *this = &(object.GetFieldData<globalIndex>( fieldName) );  break;
          case FieldInfo::realField:        *this = &(object.GetFieldData<realT>( fieldName) );       break;
          case FieldInfo::R1TensorField:    *this = &(object.GetFieldData<R1Tensor>(fieldName) );     break;
          case FieldInfo::R2TensorField:    *this = &(object.GetFieldData<R2Tensor>( fieldName) );    break;
          case FieldInfo::R2SymTensorField: *this = &(object.GetFieldData<R2SymTensor>( fieldName) ); break;
          case FieldInfo::integerParameter:
          case FieldInfo::realParameter:
          case FieldInfo::numFieldTypes:
            throw GPException("Error FieldTypeMultiPtr::SetFieldPtr: Attempt to set "+ FieldTypeName( fieldType ) + " multipointer to field."); break;
          default:
            throw GPException("Error FieldTypeMultiPtr::SetFieldPtr: Unrecognized field type.");
        }
    };
    void SetFieldPtr(realT* scalarPtr){
      *this = scalarPtr;
    }
    void SetFieldPtr(int* integerPtr){
      *this = integerPtr;
    }

    /// Copy the the contents at entry i in the field into the iterator's container.
    /// Returns an iterator to the next element in the range.
    template<class OutputIterator>
    OutputIterator CopyValues(const localIndex i, OutputIterator itr){
      switch(m_fieldType){
        case FieldInfo::integerField:  *itr++ = (*m_iArrayPtr)[i]; break;
        case FieldInfo::localIndexField:  *itr++ = (*m_lArrayPtr)[i]; break;
        case FieldInfo::globalIndexField:  *itr++ = (*m_gArrayPtr)[i]; break;
        case FieldInfo::realField:     *itr++ = (*m_sArrayPtr)[i]; break;
        case FieldInfo::R1TensorField:  itr = std::copy( (*m_vArrayPtr)[i].begin(),  (*m_vArrayPtr)[i].end(), itr); break;
        case FieldInfo::R2TensorField:  itr = std::copy( (*m_tArrayPtr)[i].begin(),  (*m_tArrayPtr)[i].end(), itr); break;
        case FieldInfo::R2SymTensorField:  itr = std::copy( (*m_stArrayPtr)[i].begin(), (*m_stArrayPtr)[i].end(), itr); break;
        case FieldInfo::integerParameter: *itr++ = *m_integerPtr; break;
        case FieldInfo::realParameter:    *itr++ = *m_scalarPtr; break;
        case FieldInfo::numFieldTypes:
        default:
            throw GPException("Error FieldTypeMultiPtr::CopyValues: Unrecognized field type.");

      }
      return itr;
    }

    // data
    FieldType m_fieldType;

  private:

    union {
      Array1dT<int>*         m_iArrayPtr;
      lArray1d*              m_lArrayPtr;
      gArray1d*              m_gArrayPtr;
      Array1dT<realT>*       m_sArrayPtr;
      Array1dT<R1Tensor>*    m_vArrayPtr;
      Array1dT<R2Tensor>*    m_tArrayPtr;
      Array1dT<R2SymTensor>* m_stArrayPtr;
      realT*                 m_scalarPtr;
      int*                   m_integerPtr;
    };
};

#endif /* FIELDTYPEMULTIPTR_H_ */
