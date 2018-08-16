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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ViewWrapper.hpp
 *
 * @date Created on: Jun 8, 2016
 * @authors settgast
 */

#ifndef GEOSX_DATAREPOSITORY_WRAPPERVIEW_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPERVIEW_HPP_

#include "ViewWrapperBase.hpp"

#include "KeyNames.hpp"
#include "common/integer_conversion.hpp"
#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "SFINAE_Macros.hpp"
#include <type_traits>

#include "Macros.hpp"
#include "Buffer.hpp"
#include "BufferOps.hpp"
#include "RestartFlags.hpp"

#include "codingUtilities/GeosxTraits.hpp"
#include "common/GeosxConfig.hpp"


#ifdef USE_ATK
#include "sidre/sidre.hpp"
#include "sidre/SidreTypes.hpp"
#include "SidreWrapper.hpp"
#endif

#include <cstdlib>


namespace geosx
{

namespace dataRepository
{

/**
 * Templated class to serve as a wrapper to arbitrary objects.
 * @tparam T is any object that is to be wrapped by ViewWrapper
 */
template< typename T >
class ViewWrapper : public ViewWrapperBase
{

public:
  /**
   * @param name name of the object
   * @param parent parent group which owns the ViewWrapper
   */
  explicit ViewWrapper( std::string const & name,
                        ManagedGroup * const parent ):
    ViewWrapperBase(name, parent),
    m_ownsData( true ),
    m_data( new T() )
  {}

  /**
   * @param name name of the object
   * @param parent parent group that owns the ViewWrapper
   * @param object object that is being wrapped by the ViewWrapper
   */
  explicit ViewWrapper( std::string const & name,
                        ManagedGroup * const parent,
                        std::unique_ptr<T> object ):
    ViewWrapperBase(name, parent),
  m_ownsData( true ),
  m_data( object.release() )
  {}

  /**
   * @param name name of the object
   * @param parent parent group that owns the ViewWrapper
   * @param object object that is being wrapped by the ViewWrapper
   */
  explicit ViewWrapper( std::string const & name,
                        ManagedGroup * const parent,
                        T * object,
                        bool takeOwnership):
    ViewWrapperBase(name,parent),
    m_ownsData( takeOwnership ),
    m_data( object )
  {}

  /**
   * default destructor
   */
  virtual ~ViewWrapper() noexcept override final
  {
    if( m_ownsData )
    {
      delete m_data;
    }
  }

  /**
   * Copy Constructor
   * @param source source for the copy
   */
  ViewWrapper( ViewWrapper const & source ):
    ViewWrapperBase("copy_constructor_test", nullptr),
    m_data(source.m_data)
  {}

  /**
   * Move Constructor
   * @param source source to be moved
   */
  ViewWrapper( ViewWrapper&& source ):
    ViewWrapperBase(source),
    m_data( std::move(source.m_data) )
  {}

  /**
   * Copy Assignment Operator
   * @param source rhs
   * @return *this
   */
  ViewWrapper& operator=( ViewWrapper const & source )
  {
    m_data = source.m_data;
    return *this;
  }

  /**
   * Move Assignment Operator
   * @param source
   * @return *this
   */
  ViewWrapper& operator=( ViewWrapper && source )
  {
    m_data = std::move(source.m_data);
    return *this;
  }


  /**
   * Factory Method to make a new ViewWrapper<T>, allocating a new T. Only is
   * going to work if T has a default constructor.
   * Perhaps this is worthless in the general case.
   * @param name name of the object
   * @param parent group that owns the ViewWrapper
   * @return A std::unique_ptr<ViewWrapperBase> that holds the newly allocated
   * ViewWrapper.
   */
  template<typename TNEW>
  static std::unique_ptr<ViewWrapperBase> Factory( std::string const & name,
                                                   ManagedGroup * const parent)
  {
    std::unique_ptr<TNEW> newObject = std::make_unique<TNEW>();
    return std::make_unique<ViewWrapper<T> >( name, parent, std::move(newObject));
  }

  virtual std::unique_ptr<ViewWrapperBase> clone( string const & name,
                                                  ManagedGroup * const parent ) override
  {
    std::unique_ptr<ViewWrapperBase>
    clonedWrapper = std::make_unique<ViewWrapper<T> >( name, parent, this->m_data, false );

    return clonedWrapper;
  }



  /**
   * Virtual function to return the typeid of T. Not so sure this does what we
   * want?? TODO
   * @return typeid(T)
   */
  virtual const std::type_info& get_typeid() const noexcept override final
  {
    return typeid(T);
  }

  /**
   * static function to cast a ViewWrapper base to a derived ViewWrapper<T>
   * @param base
   * @return casted ViewWrapper<T>
   */
  static ViewWrapper<T>& cast( ViewWrapperBase& base )
  {
    if( base.get_typeid() != typeid(T) )
    {
      GEOS_ERROR("invalid cast attempt");
    }
    return static_cast< ViewWrapper<T>& >(base);
  }

  static ViewWrapper<T> const & cast( ViewWrapperBase const & base )
  {
    if( base.get_typeid() != typeid(T) )
    {
      GEOS_ERROR("invalid cast attempt");
    }
    return static_cast< ViewWrapper<T> const & >(base);
  }

//  template< bool DO_PACKING >
//  struct pack_wrapper
//  {
//
////    template<typename TT >
////      struct has_square_bracket_operator
////      {
////    private:
////        template<typename U> static constexpr auto test(int)->decltype( std::is_member_function_pointer<decltype(&U::operator[])>::value, bool() )
////        {
////          return std::is_member_function_pointer<decltype(&U::operator[])>::value;
////        }
////        template<typename U> static constexpr auto test(...)->bool
////        {
////          return false;
////        }
////    public:
////        static constexpr bool value = test<TT>(0);
////      };
//
//
//
//    HAS_ALIAS(pointer)
//    HAS_ALIAS(iterator)
//    HAS_MEMBER_FUNCTION(begin,typename U::iterator,,,)
//    HAS_MEMBER_FUNCTION(end,typename U::iterator,,,)
//
//
//    template<class U = T>
//    static typename std::enable_if< CommBufferOps::is_packable<U>::value, localIndex >::type
//    PackT( char *& buffer, U const & m_data )
//    {
//      return CommBufferOps::Pack<DO_PACKING>( buffer, m_data );
//    }
//
//    template<class U = T>
//    static typename std::enable_if< !CommBufferOps::is_packable<U>::value, localIndex >::type
//    PackT( char *& buffer, U const & m_data )
//    {
//      return 0;
//    }
//
//    template<class U = T>
//    static typename std::enable_if< CommBufferOps::is_packable<U>::value, localIndex >::type
//    UnpackT( char const *& buffer, U & m_data )
//    {
//      return CommBufferOps::Unpack( buffer, m_data );
//    }
//
//    template<class U = T>
//    static typename std::enable_if< !CommBufferOps::is_packable<U>::value, localIndex >::type
//    UnpackT( char const *& buffer, U & m_data )
//    {
//      return 0;
//    }
//
//
//
//
//
//    template<class U = T>
//    static typename std::enable_if< CommBufferOps::is_packable<U>::value &&
//                                    !std::is_same<U,string>::value, localIndex >::type
//    PackT( char *& buffer, U const & m_data, localIndex_array const & packList )
//    {
//      return CommBufferOps::Pack<DO_PACKING>( buffer, m_data, packList );
//    }
//
//    template<class U = T>
//    static typename std::enable_if< CommBufferOps::is_packable<U>::value &&
//                                    std::is_same<U,string>::value, localIndex >::type
//    PackT( char *& buffer, U const & m_data, localIndex_array const & packList )
//    {
//      return CommBufferOps::Pack<DO_PACKING>( buffer, m_data );
//    }
//
//    template<class U = T>
//    static typename std::enable_if< !CommBufferOps::is_packable<U>::value , localIndex >::type
//    PackT( char *& buffer, U const & m_data, localIndex_array const & packList )
//    {
//      return 0;
//    }
//
//
//
//
//
//    template<class U = T>
//    static typename std::enable_if< CommBufferOps::is_packable<U>::value &&
//                                    !std::is_same<U,string>::value, localIndex >::type
//    UnpackT( char const *& buffer, U & m_data, localIndex_array const & packList )
//    {
//      return CommBufferOps::Unpack( buffer, m_data, packList );
//    }
//
//    template<class U = T>
//    static typename std::enable_if< CommBufferOps::is_packable<U>::value &&
//                                    std::is_same<U,string>::value, localIndex >::type
//    UnpackT( char const *& buffer, U & m_data, localIndex_array const & packList )
//    {
//      return CommBufferOps::Unpack( buffer, m_data );
//    }
//
//    template<class U = T>
//    static typename std::enable_if< !CommBufferOps::is_packable<U>::value , localIndex >::type
//    UnpackT( char const *& buffer, U & m_data, localIndex_array const & packList )
//    {
//      return 0;
//    }
//
//
//  };

  virtual localIndex Pack( char *& buffer ) const override final
  {
    localIndex packedSize = 0;

    packedSize += bufferOps::Pack<true>( buffer, this->getName() );
    packedSize += bufferOps::Pack<true>( buffer, *m_data);

    return packedSize;
  }

  virtual localIndex Pack( char *& buffer, localIndex_array const & packList ) const override final
  {
    localIndex packedSize = 0;

    static_if( bufferOps::is_packable_by_index<T>::value )
    {
      packedSize += bufferOps::Pack<true>( buffer, this->getName() );
      packedSize += bufferOps::Pack<true>( buffer, *m_data, packList);
    });
    return packedSize;
  }

  virtual localIndex PackSize( ) const override final
  {
    char * buffer = nullptr;
    localIndex packedSize = 0;

    packedSize += bufferOps::Pack<false>( buffer, this->getName() );
    packedSize += bufferOps::Pack<false>( buffer, *m_data);

    return packedSize;
  }

  virtual localIndex PackSize( localIndex_array const & packList ) const override final
  {

    char * buffer = nullptr;
    localIndex packedSize = 0;

    static_if( bufferOps::is_packable_by_index<T>::value )
    {
      packedSize += bufferOps::Pack<false>( buffer, this->getName() );
      packedSize += bufferOps::Pack<false>( buffer, *m_data, packList);
    });

    return packedSize;
  }

  virtual localIndex Unpack( char const *& buffer ) override final
  {
    localIndex unpackedSize = 0;
    string name;
    unpackedSize += bufferOps::Unpack( buffer, name );
    GEOS_ASSERT( name == this->getName(),"buffer unpack leads to viewWrapper names that don't match" )
    unpackedSize += bufferOps::Unpack( buffer, *m_data );
    return unpackedSize;
  }
  virtual localIndex Unpack( char const *& buffer, localIndex_array const & unpackIndices ) override final
  {
    localIndex unpackedSize = 0;
    static_if( bufferOps::is_packable_by_index<T>::value )
    {
      string name;
      unpackedSize += bufferOps::Unpack( buffer, name );
      GEOS_ASSERT( name == this->getName(),"buffer unpack leads to viewWrapper names that don't match" )
      unpackedSize += bufferOps::Unpack( buffer, *m_data, unpackIndices );
    });
    return unpackedSize;
  }




  struct empty_wrapper
  {
    HAS_MEMBER_FUNCTION(empty,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_empty<U>::value, bool>::type
    empty(ViewWrapper<T> const * parent)
    {
      return parent->m_data->empty();
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_empty<U>::value, bool>::type
    empty(ViewWrapper<T> const * parent)
    {
      return parent;
    }
  };
  virtual bool empty() const override final
  {
    return empty_wrapper::empty(this);
  }

  struct size_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT(size, 0, int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT(size, 1, unsigned int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT(size, 2, long int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT(size, 3, unsigned long int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT(size, 4, long long, const, , )
    HAS_MEMBER_FUNCTION_VARIANT(size, 5, unsigned long long, const, , )

    template<class U = T>
    static typename std::enable_if< has_memberfunction_v0_size<U>::value ||
                                    has_memberfunction_v1_size<U>::value ||
                                    has_memberfunction_v2_size<U>::value ||
                                    has_memberfunction_v3_size<U>::value ||
                                    has_memberfunction_v4_size<U>::value ||
                                    has_memberfunction_v5_size<U>::value, localIndex>::type
    size(ViewWrapper<T> const * parent)
    { return static_cast<localIndex>(parent->m_data->size()); }
    
    template<class U = T>
    static typename std::enable_if< !(has_memberfunction_v0_size<U>::value ||
                                      has_memberfunction_v1_size<U>::value ||
                                      has_memberfunction_v2_size<U>::value ||
                                      has_memberfunction_v3_size<U>::value ||
                                      has_memberfunction_v4_size<U>::value ||
                                      has_memberfunction_v5_size<U>::value), localIndex>::type
    size(ViewWrapper<T> const * parent)
    { return 1; }
  };
  virtual localIndex size() const override final
  {
    return size_wrapper::size(this);
  }


  struct num_dimensions_wrapper
  {
    HAS_MEMBER_FUNCTION(numDimensions,int,const,,)

    template<class U = T>
    static typename std::enable_if<has_memberfunction_numDimensions<U>::value, int>::type
    numDimensions(ViewWrapper<T> const * parent)
    { return static_cast<int>(parent->m_data->numDimensions()); }
    
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_numDimensions<U>::value, int>::type
    numDimensions(ViewWrapper<T> const * parent)
    { return 1; }
  };
  virtual int numDimensions() const override final
  {
    return num_dimensions_wrapper::numDimensions(this);
  }


  struct dimension_size_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT(size, 0, int, const, VA_LIST(int), VA_LIST(int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT(size, 1, unsigned int, const, VA_LIST(int), VA_LIST(int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT(size, 2, long int, const, VA_LIST(int), VA_LIST(int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT(size, 3, unsigned long int, , VA_LIST(int), VA_LIST(int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT(size, 4, long long, const, VA_LIST(int), VA_LIST(int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT(size, 5, unsigned long long, const, VA_LIST(int), VA_LIST(int(1)) )

    template<class U = T>
    static typename std::enable_if< has_memberfunction_v0_size<U>::value ||
                                    has_memberfunction_v1_size<U>::value ||
                                    has_memberfunction_v2_size<U>::value ||
                                    has_memberfunction_v3_size<U>::value ||
                                    has_memberfunction_v4_size<U>::value ||
                                    has_memberfunction_v5_size<U>::value, localIndex>::type
    size(ViewWrapper<T> const * const parent, int const i)
    { return integer_conversion<localIndex>(parent->m_data->size(i)); }
    
    template<class U = T>
    static typename std::enable_if< !(has_memberfunction_v0_size<U>::value ||
                                      has_memberfunction_v1_size<U>::value ||
                                      has_memberfunction_v2_size<U>::value ||
                                      has_memberfunction_v3_size<U>::value ||
                                      has_memberfunction_v4_size<U>::value ||
                                      has_memberfunction_v5_size<U>::value), localIndex>::type
    size(ViewWrapper<T> const * const parent, int const i)
    { 
      if (i != 0) 
      {
        GEOS_ERROR("Data is only 1D");
        return 0;
      }
      return parent->size(); 
    }
  };
  virtual localIndex size( int const i) const override final
  {
    return dimension_size_wrapper::size(this, i);
  }


  struct resize_dimension_wrapper
  {
    HAS_MEMBER_FUNCTION(resize, void, , VA_LIST(int, long long const * const), VA_LIST(int(1), nullptr))

    template<class U=T>
    static typename std::enable_if<has_memberfunction_resize<U>::value, void>::type
    resize(ViewWrapper<T> * parent, int num_dims, long long const * const dims)
    { parent->m_data->resize(num_dims, dims); }

    template<class U=T>
    static typename std::enable_if<!has_memberfunction_resize<U>::value, void>::type
    resize(ViewWrapper<T> * parent, int num_dims, long long const * const dims)
    {
      if (num_dims != 1)
      {
        GEOS_ERROR("Data is only 1D");
        return;
      }
      parent->resize( integer_conversion<localIndex>(dims[0]));
    }
  };
  virtual void resize(int num_dims, long long const *  const dims) override final
  { resize_dimension_wrapper::resize(this, num_dims, dims); }


  struct reserve_wrapper
  {
    HAS_MEMBER_FUNCTION(reserve, void, ,VA_LIST(std::size_t),VA_LIST(std::size_t(1)) )
    template<class U = T>
    static typename std::enable_if<has_memberfunction_reserve<U>::value, void>::type reserve(ViewWrapper<T> * const parent, std::size_t new_cap)
    {
      return parent->m_data->reserve(new_cap);
    }
    template<class U = T>
    static typename std::enable_if<!has_memberfunction_reserve<U>::value, void>::type reserve(ViewWrapper<T> * const, std::size_t )
    {
      return; //parent->m_data;
    }
  };
  virtual void reserve( std::size_t new_cap ) override final
  {
    reserve_wrapper::reserve(this, new_cap);
  }
//  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,reserve , void,,
// VA_LIST(std::size_t a), VA_LIST(a) )


  HAS_MEMBER_FUNCTION(capacity,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,capacity,std::size_t,const)

  HAS_MEMBER_FUNCTION(max_size,std::size_t,const,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,max_size,std::size_t,const)

  HAS_MEMBER_FUNCTION(clear,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,clear,void,)

  HAS_MEMBER_FUNCTION(insert,void,,,)
  CONDITIONAL_VIRTUAL_FUNCTION0(ViewWrapper<T>,insert,void,)


  struct resize_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT(resize,0,void,,VA_LIST(int), VA_LIST(static_cast<int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,1,void,,VA_LIST(unsigned int), VA_LIST( static_cast<unsigned int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,2,void,,VA_LIST(long), VA_LIST( static_cast<long int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,3,void,,VA_LIST(unsigned long), VA_LIST(static_cast<unsigned long int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,4,void,,VA_LIST(long long int), VA_LIST(static_cast<long long int>(1)))
    HAS_MEMBER_FUNCTION_VARIANT(resize,5,void,,VA_LIST(unsigned long long), VA_LIST(static_cast<unsigned long long>(1)))


    template<class U = T>
    static typename std::enable_if< has_memberfunction_v0_resize<U>::value ||
                                    has_memberfunction_v1_resize<U>::value ||
                                    has_memberfunction_v2_resize<U>::value ||
                                    has_memberfunction_v3_resize<U>::value ||
                                    has_memberfunction_v4_resize<U>::value ||
                                    has_memberfunction_v5_resize<U>::value, void>::type
    resize(ViewWrapper * const parent, localIndex const new_size)
    {
      return parent->m_data->resize(new_size);
    }


    template<class U = T>
    static typename std::enable_if< !(has_memberfunction_v0_resize<U>::value ||
                                      has_memberfunction_v1_resize<U>::value ||
                                      has_memberfunction_v2_resize<U>::value ||
                                      has_memberfunction_v3_resize<U>::value ||
                                      has_memberfunction_v4_resize<U>::value ||
                                      has_memberfunction_v5_resize<U>::value), void>::type
    resize(ViewWrapper * const, localIndex )
    {
      return;
    }
  };
  virtual void resize( localIndex new_size ) override final
  {
    resize_wrapper::resize(this, new_size);
  }


  struct should_resize_wrapper
  {
    HAS_MEMBER_FUNCTION(isSorted,bool,const,,)
    template<class U = T>
    static typename std::enable_if<has_memberfunction_isSorted<U>::value, bool>::type shouldResize()
    { return false;  }

    template<class U = T>
    static typename std::enable_if<!has_memberfunction_isSorted<U>::value, bool>::type shouldResize()
    { return true; }
  };
  virtual bool shouldResize() const override final
  {
    return should_resize_wrapper::shouldResize();
  }

  /**
   * @name Structure to determine return types for data access functions
   */
  ///@{

  /// Invoke macro to generate test to see if type has an alias named "pointer".
  /// This will be used to determine if the
  /// type is to be treated as an "array" or a single object.
  HAS_ALIAS(pointer)

  /**
   * SFINAE specialized structure to control return type based on properties of
   * T.
   * The default template returns a pointer for all calls to data().
   */
  template< class U=T,
            bool HASPOINTERTYPE = has_alias_pointer<U>::value,
            bool ISSTRING = std::is_same<U,std::string>::value >
  struct Get_Type
  {
    typedef U*       type;
    typedef U const * const_type;

    typedef U *       pointer;
    typedef U const * const_pointer;
  };

  /**
   *  Specialization for case when T has a pointer alias, and it is NOT a
   * string.
   *  In this case, we assume that we are storing an array type. The return type
   * is then a reference, unless the
   *  compilation flag is set such that we require a pointer back (good for
   * speed, but no array class convenience).
   *  The resulting types can both be dereferenced with operator[], so no code
   * changes required
   *  unless array member functions have been called.
   */
  template<class U>
  struct Get_Type<U, true, false>
  {

#ifdef USE_CONTAINERARRAY_RETURN_PTR
    typedef typename U::pointer       type;
    typedef typename U::const_pointer const_type;
#else
    typedef U &       type;
    typedef U const & const_type;
#endif
    typedef typename U::pointer       pointer;
    typedef typename U::const_pointer const_pointer;
  };


  /// Specialization for string. Always return a reference.
  template<class U>
  struct Get_Type<U, true, true>
  {
    typedef U &       type;
    typedef U const & const_type;

    typedef U *       pointer;
    typedef U const * const_pointer;
  };

  using rtype       = typename Get_Type<T>::type;
  using rtype_const = typename Get_Type<T>::const_type;

  using pointer       = typename Get_Type<T>::pointer;
  using const_pointer = typename Get_Type<T>::const_pointer;
  ///@}


  HAS_MEMBER_FUNCTION(data,pointer,,,)
  HAS_MEMBER_FUNCTION_VARIANT(data,_const, pointer,const,,)

  /// Case for if m_data has a member function called "data()", and is not a
  // string
  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && !std::is_same<U,std::string>::value, rtype>::type
  data()
  {
#ifdef USE_CONTAINERARRAY_RETURN_PTR
    return m_data->data();
#else
    return *m_data;
#endif
  }


  template<class U = T>
  typename std::enable_if<has_memberfunction_data<U>::value && !std::is_same<U,string>::value, rtype_const>::type
  data() const
  {
#ifdef USE_CONTAINERARRAY_RETURN_PTR
    return m_data->data();
#else
    return *m_data;
#endif
  }


  /// Case for if m_data is a string
  template<class U = T>
  typename std::enable_if<std::is_same<U,std::string>::value, rtype>::type
  data()
  {
    /// return the object...or a reference to the object
    return *m_data;
  }


  template<class U = T>
  typename std::enable_if<std::is_same<U,std::string>::value, rtype_const>::type
  data() const
  {
    return *m_data;
  }


  /// case for if m_data does NOT have a member function "data()", and is not a
  // string
  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value && !std::is_same<U,std::string>::value, rtype>::type
  data()
  {
    /// return a c-pointer to the object
    return m_data;
  }


  template<class U = T>
  typename std::enable_if<!has_memberfunction_data<U>::value && !std::is_same<U,std::string>::value, rtype_const>::type
  data() const
  {
    return m_data;
  }


  T& reference()
  { return *m_data; }

  T const & reference() const
  { return *m_data; }



  /// Case for if m_data has a member function called "data()"
  template<class U = T>
  typename std::enable_if< ( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value ) &&
                           has_alias_pointer<U>::value && !std::is_same<U,string>::value,typename U::pointer >::type
  dataPtr()
  {
    return m_data->data();
  }

  template<class U = T>
  typename std::enable_if< ( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value ) &&
                           has_alias_pointer<U>::value && !std::is_same<U,string>::value,typename U::const_pointer >::type
  dataPtr() const
  {
    return m_data->data();
  }


  /// Case for if m_data is a string"
  template<class U = T>
  typename std::enable_if< std::is_same<U,string>::value, char const * >::type
  dataPtr()
  {
    return m_data->data();
  }

  template<class U = T>
  typename std::enable_if< std::is_same<U,string>::value, char const * >::type
  dataPtr() const
  {
    return m_data->data();
  }


  /// case for if m_data does NOT have a member function "data()"
  template<class U = T>
  typename std::enable_if<!( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value )&&
                          !std::is_same<U,string>::value, U * >::type
  dataPtr()
  {
    return m_data;
  }

  template<class U = T>
  typename std::enable_if<!( has_memberfunction_data<U>::value || has_memberfunction_v_const_data<U>::value )&&
                          !std::is_same<U,string>::value, U const *>::type
  dataPtr() const
  {
    return m_data;
  }

  HAS_ALIAS(value_type)


  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, size_t>::type
  sizeOfValueType() const
  {
    return sizeof(typename T::value_type);
  }

  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, size_t>::type
  sizeOfValueType() const
  {
    return sizeof(T);
  }

  virtual size_t sizeOfType() const override final
  {
    return sizeOfValueType();
  }


  /// case for if U::value_type exists. Returns the size of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  byteSize() const
  {
    return size() * sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the size of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  byteSize() const
  {
    return size() * sizeof(T);
  }


  /// case for if U::value_type exists. Returns the size of an element of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  elementSize() const
  {
    return sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the size of an element of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  elementSize() const
  {
    return sizeof(T);
  }


  /// case for if U::value_type exists. Returns the typeid of an element of dataPtr
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, const std::type_info&>::type
  elementTypeID() const
  {
    return typeid(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the typeid of an element of dataPtr
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, const std::type_info&>::type
  elementTypeID() const
  {
    return typeid(T);
  }


  /// case for if U::value_type exists. Returns the number of elements given a byte size
  template<class U = T>
  typename std::enable_if<has_alias_value_type<U>::value, localIndex>::type
  numElementsFromByteSize(localIndex d_size) const
  {
    return d_size / sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the number of elements
  // given a byte size
  template<class U = T>
  typename std::enable_if<!has_alias_value_type<U>::value, localIndex>::type
  numElementsFromByteSize(localIndex d_size) const
  {
    return d_size / sizeof(T);
  }


  virtual bool shouldRegisterDataPtr() const override
  {
    std::type_index type_index = std::type_index(elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
    return sidre_type_id != axom::sidre::TypeID::NO_TYPE_ID;
  }


  
  void registerDataPtr(axom::sidre::View * view) const override
  {
#ifdef USE_ATK
    view = (view != nullptr) ? view : getSidreView();

    localIndex num_elements = size();
    if (num_elements > 0)
    {
      std::type_index type_index = std::type_index(elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
      if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
      { return; }

      localIndex sidre_size = rtTypes::getSidreSize(type_index);
      localIndex byte_size = byteSize();
      localIndex element_size = elementSize();

      int ndims = numDimensions();
//      axom::sidre::SidreLength dims[ndims + 1];
      axom::sidre::SidreLength dims[10];
      for (int dim = 0; dim < ndims; ++dim)
      {
        dims[dim] = size(dim);
      }

      if ( byte_size > num_elements * sidre_size )
      {
        dims[ndims++] = element_size / sidre_size;
      }
      
      void * ptr = const_cast<void*>( static_cast<void const *>( dataPtr() ) );
      view->setExternalDataPtr(sidre_type_id, ndims, dims, ptr);
    }
    else
    {
      unregisterDataPtr(view);
    }
#endif
  }  

  /* Register the pointer to data with the associated sidre::View. */
  void registerToWrite(axom::sidre::View * view=nullptr) const override
  {
#ifdef USE_ATK
    if (getRestartFlags() == RestartFlags::NO_WRITE)
    {
      view = (view != nullptr) ? view : getSidreView();
      unregisterDataPtr(view);
      return;
    }

    view = (view != nullptr) ? view : getSidreView();
    storeSizedFromParent(view);

    localIndex num_elements = size();
    if (num_elements > 0) 
    {
      std::type_index type_index = std::type_index(elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
      if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
      {
        localIndex byte_size;
        void * ptr = Buffer::pack(reference(), byte_size);
        view->setExternalDataPtr(axom::sidre::TypeID::INT8_ID, byte_size, ptr);
        return;
      }

      localIndex sidre_size = rtTypes::getSidreSize(type_index);
      localIndex byte_size = byteSize();
      localIndex element_size = elementSize();

      int ndims = numDimensions();
//      axom::sidre::SidreLength dims[ndims + 1];
      axom::sidre::SidreLength dims[10];
      for (int dim = 0; dim < ndims; ++dim)
      {
        dims[dim] = size(dim);
      }

      if ( byte_size > num_elements * sidre_size )
      {
        dims[ndims++] = element_size / sidre_size;
      }
      
      void * ptr = const_cast<void*>(static_cast<void const *>( dataPtr() ) );
      view->setExternalDataPtr(sidre_type_id, ndims, dims, ptr);
    }
    else
    {
      unregisterDataPtr(view);
    }
#endif
  }

  /* Register the pointer to data with the associated sidre::View. */
  void finishWriting(axom::sidre::View * view=nullptr) const override
  {
#ifdef USE_ATK
    if (getRestartFlags() == RestartFlags::NO_WRITE)
    {
      view = (view != nullptr) ? view : getSidreView();
      unregisterDataPtr(view);
      return;
    }

    view = (view != nullptr) ? view : getSidreView();
    view->setAttributeToDefault("__sizedFromParent__");

    if (!view->isExternal() || view->getTotalBytes() == 0)
    {
      return;
    }
    
    std::type_index type_index = std::type_index(elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
    if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
    {
      std::free(view->getVoidPtr());
    }

    unregisterDataPtr(view);
#endif
  }

  void registerToRead(axom::sidre::View * view=nullptr) override
  {
#ifdef USE_ATK
    if (getRestartFlags() != RestartFlags::WRITE_AND_READ)
    {
      unregisterDataPtr(view);
      return;
    }

    view = (view != nullptr) ? view : getSidreView();
    loadSizedFromParent(view);
    if (!view->isExternal() || view->getTotalBytes() == 0)
    {
      return;
    }
    
    std::type_index type_index = std::type_index(elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
    if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
    {
      localIndex byte_size = view->getTotalBytes();
      void * ptr = std::malloc(byte_size);
      view->setExternalDataPtr(axom::sidre::TypeID::INT8_ID, byte_size, ptr);
      return;
    }

    resizeFromSidre(view);
    void * ptr = const_cast<void*>( static_cast<void const *>( dataPtr() ) );
    localIndex sidre_size = rtTypes::getSidreSize(type_index);
    view->setExternalDataPtr(sidre_type_id, byteSize() / sidre_size, ptr);
#endif
  }


  void finishReading(axom::sidre::View * view) override
  {
#ifdef USE_ATK
    if (getRestartFlags() != RestartFlags::WRITE_AND_READ)
    {
      view = (view != nullptr) ? view : getSidreView();
      unregisterDataPtr(view);
      return;
    }
    view = (view != nullptr) ? view : getSidreView();
    if (!view->isExternal() || view->getTotalBytes() == 0)
    {
      return;
    }
    
    std::type_index type_index = std::type_index(elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType(type_index);
    if (sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID)
    {
      localIndex byte_size = view->getTotalBytes();
      void * ptr = view->getVoidPtr();
      Buffer::unpack(reference(), ptr, byte_size);
      std::free(ptr);
    }

    unregisterDataPtr(view);
#endif
  }

  void unregisterDataPtr(axom::sidre::View* view = nullptr) const
  {
#ifdef USE_ATK
    view = (view != nullptr) ? view : getSidreView();
    view->setExternalDataPtr(AXOM_NULLPTR);
#endif
  }

  void storeSizedFromParent(axom::sidre::View* view = nullptr) const
  {
#ifdef USE_ATK
    if (SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
    {
      view = (view != nullptr) ? view : getSidreView();
      view->setAttributeScalar("__sizedFromParent__", sizedFromParent());
    }
#endif
  }

  void loadSizedFromParent(axom::sidre::View* view = nullptr)
  {
#ifdef USE_ATK
    if (SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
    {
      view = (view != nullptr) ? view : getSidreView();
      setSizedFromParent(view->getAttributeScalar("__sizedFromParent__"));
      view->setAttributeToDefault("__sizedFromParent__");
    }
#endif
  }

  void resizeFromSidre(axom::sidre::View* view = nullptr)
  {
#ifdef USE_ATK
    view = (view != nullptr) ? view : getSidreView();
    if (view->isExternal()) 
    { 
      std::type_index type_index = std::type_index(elementTypeID());
      localIndex sidre_size = rtTypes::getSidreSize(type_index);

      localIndex byte_size = view->getTotalBytes();
      localIndex num_elements = numElementsFromByteSize(byte_size);

      int ndims = view->getNumDimensions();
//      axom::sidre::SidreLength dims[ndims];
      axom::sidre::SidreLength dims[10];
      view->getShape(ndims, dims);

      if ( byte_size > num_elements * sidre_size )
      {
        ndims--;
      }

      localIndex num_elems_recorded = 1;
      for (localIndex i = 0; i < ndims; ++i)
      {
        num_elems_recorded *= dims[i];
      }

      if (num_elems_recorded != num_elements)
      {
        GEOS_ERROR("Number of elements recorded not equal to the calculated number: " << 
                   num_elems_recorded << " " << num_elements);
      }

//      long long l_dims[ndims];
      long long l_dims[10];
      for (localIndex i = 0; i < ndims; ++i)
      {
        l_dims[i] = dims[i];
      }

      resize(ndims, l_dims);
    }
#endif
  }


  bool m_ownsData;
  T * m_data;

  ViewWrapper() = delete;
};

template< typename T >
using view_rtype = typename ViewWrapper<T>::rtype;

template< typename T >
using view_rtype_const = typename ViewWrapper<T>::rtype_const;

}
} /* namespace geosx */

#endif /* CORE_SRC_DATAREPOSITORY_DATAOBJECT_HPP_ */
