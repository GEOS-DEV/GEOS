/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
/**
 * @file Object.h
 * @author settgast1
 * @date Nov 5, 2010
 */

#ifndef OBJECT_H_
#define OBJECT_H_

#include "Common/Common.h"
#include "FieldRegistry.h"

class Object
{
public:
  Object();
  ~Object();

  Object( const Object& init );

  Object& operator=( const Object& rhs );
  Object& operator++();



  FieldRegistry* m_fieldRegistry;


  template< class T >
  T& GetFieldFromOffset( const size_t offset )
  {
    return (T&)(*((realT*)this + offset));
  }

  template< class T >
  const T& GetFieldFromOffset( const size_t offset ) const
  {
    return (const T&)(*((realT*)this + offset));
  }


  template< FieldKey FIELD >
  inline typename Field<FIELD>::Type& GetField()
  { return (typename Field<FIELD>::Type&)(*((realT*)this + m_fieldRegistry->m_fieldEnumToOffset[FIELD])); }

  template< FieldKey FIELD >
  inline const typename Field<FIELD>::Type& GetField() const
  { return (typename Field<FIELD>::Type&)(*((realT*)this + m_fieldRegistry->m_fieldEnumToOffset[FIELD])); }



  inline Field<FieldInfo::referencePosition>::Type& ReferencePosition()
  { return GetField<FieldInfo::referencePosition>(); }

  inline Field<FieldInfo::displacement>::Type& Displacement()
  { return GetField<FieldInfo::displacement>(); }

  inline Field<FieldInfo::incrementalDisplacement>::Type& IncrementalDisplacement()
  { return GetField<FieldInfo::incrementalDisplacement>(); }

  inline Field<FieldInfo::velocity>::Type& Velocity()
  { return GetField<FieldInfo::velocity>(); }

  inline Field<FieldInfo::acceleration>::Type& Acceleration()
  { return GetField<FieldInfo::acceleration>(); }

  inline Field<FieldInfo::force>::Type& Force()
  { return GetField<FieldInfo::force>(); }

  inline Field<FieldInfo::mass>::Type&   Mass()
  { return GetField<FieldInfo::mass>(); }


  inline const Field<FieldInfo::referencePosition>::Type& ReferencePosition() const
  { return GetField<FieldInfo::referencePosition>(); }

  inline const Field<FieldInfo::displacement>::Type& Displacement() const
  { return GetField<FieldInfo::displacement>(); }

  inline const Field<FieldInfo::incrementalDisplacement>::Type& IncrementalDisplacement() const
  { return GetField<FieldInfo::incrementalDisplacement>(); }

  inline const Field<FieldInfo::velocity>::Type& Velocity() const
  { return GetField<FieldInfo::velocity>(); }

  inline const Field<FieldInfo::acceleration>::Type& Acceleration() const
  { return GetField<FieldInfo::acceleration>(); }

  inline const Field<FieldInfo::force>::Type& Force() const
  { return GetField<FieldInfo::force>(); }

  inline const Field<FieldInfo::mass>::Type&   Mass() const
  { return GetField<FieldInfo::mass>(); }

};



inline Object::Object():
  m_fieldRegistry(NULL)
{}

inline Object::Object(const Object& init):
  m_fieldRegistry(NULL)
{
  operator=(init);
}

inline Object::~Object()
{}


inline Object& Object::operator=( const Object& rhs )
{
  memcpy( this, &rhs, rhs.m_fieldRegistry->m_objectSize * sizeof(realT) );
  return *this;
}

inline Object& Object::operator++()
{
  return( (Object&)(*((realT*)this + m_fieldRegistry->m_objectSize)) );
}



#endif /* OBJECT_H_ */
