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
