/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Wrapper.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_WRAPPER_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPER_HPP_

#include <type_traits>

#include "KeyNames.hpp"
#include "IntegerConversion.hpp"
#include "common/DataTypes.hpp"
#include "SFINAE_Macros.hpp"

#include "Macros.hpp"
#include "BufferOps.hpp"
#include "RestartFlags.hpp"
#include "codingUtilities/GeosxTraits.hpp"
#include "common/GeosxConfig.hpp"
#include "DefaultValue.hpp"
#include "cxx-utilities/src/src/StringUtilities.hpp"
#include "WrapperBase.hpp"


#ifdef GEOSX_USE_ATK
#include "axom/sidre/core/sidre.hpp"
#include "SidreWrapper.hpp"
#endif

#include <cstdlib>


namespace geosx
{

namespace dataRepository
{

//template< typename U >
//static void totalViewType( char * const dataType );

/**
 * Templated class to serve as a wrapper to arbitrary objects.
 * @tparam T is any object that is to be wrapped by Wrapper
 */
template< typename T >
class Wrapper : public WrapperBase
{

public:

  using TYPE = T;
  /**
   * @param name name of the object
   * @param parent parent group which owns the Wrapper
   */
  explicit Wrapper( std::string const & name,
                    Group * const parent ):
    WrapperBase( name, parent ),
    m_ownsData( true ),
    m_data( new T() ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value || traits::is_string< T > )
    {
      this->setSizedFromParent( 0 );
    }

    setUserCallBack();
  }

  /**
   * @param name name of the object
   * @param parent parent group that owns the Wrapper
   * @param object object that is being wrapped by the Wrapper
   */
  explicit Wrapper( std::string const & name,
                    Group * const parent,
                    std::unique_ptr< T > object ):
    WrapperBase( name, parent ),
    m_ownsData( true ),
    m_data( object.release() ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value || traits::is_string< T > )
    {
      this->setSizedFromParent( 0 );
    }

    setUserCallBack();
  }

  /**
   * @param name name of the object
   * @param parent parent group that owns the Wrapper
   * @param object object that is being wrapped by the Wrapper
   * @param takeOwnership to indicate whether or not to take ownership of \p object
   */
  explicit Wrapper( std::string const & name,
                    Group * const parent,
                    T * object,
                    bool takeOwnership ):
    WrapperBase( name, parent ),
    m_ownsData( takeOwnership ),
    m_data( object ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value || traits::is_string< T > )
    {
      this->setSizedFromParent( 0 );
    }

    setUserCallBack();
  }

  /**
   * default destructor
   */
  virtual ~Wrapper() noexcept override final
  {
    if( m_ownsData )
    {
      delete m_data;
    }
    //tvTemplateInstantiation();
  }

  /**
   * Copy Constructor
   * @param source source for the copy
   */
  Wrapper( Wrapper const & source ):
    WrapperBase( "copy_constructor_test", nullptr ),
    m_ownsData( source.m_ownsData ),
    m_data( source.m_data ),
    m_default( source.m_default )
  {}

  /**
   * Move Constructor
   * @param source source to be moved
   */
  Wrapper( Wrapper && source ):
    WrapperBase( source ),
    m_ownsData( source.m_ownsData ),
    m_data( std::move( source.m_data ) ),
    m_default( source.m_default )
  {}

  /**
   * Copy Assignment Operator
   * @param source rhs
   * @return *this
   */
  Wrapper & operator=( Wrapper const & source )
  {
    m_data = source.m_data;
    return *this;
  }

  /**
   * Move Assignment Operator
   * @param source
   * @return *this
   */
  Wrapper & operator=( Wrapper && source )
  {
    m_data = std::move( source.m_data );
    return *this;
  }


  /**
   * Factory Method to make a new Wrapper<T>, allocating a new T. Only is
   * going to work if T has a default constructor.
   * Perhaps this is worthless in the general case.
   * @param name name of the object
   * @param parent group that owns the Wrapper
   * @return A std::unique_ptr<WrapperBase> that holds the newly allocated
   * Wrapper.
   */
  template< typename TNEW >
  static std::unique_ptr< WrapperBase > Factory( std::string const & name,
                                                 Group * const parent )
  {
    std::unique_ptr< TNEW > newObject = std::make_unique< TNEW >();
    return std::make_unique< Wrapper< T > >( name, parent, std::move( newObject ));
  }

  virtual std::unique_ptr< WrapperBase > clone( string const & name,
                                                Group * const parent ) override
  {
    std::unique_ptr< WrapperBase >
    clonedWrapper = std::make_unique< Wrapper< T > >( name, parent, this->m_data, false );
    clonedWrapper->CopyWrapperAttributes( *this );

    return clonedWrapper;
  }

  virtual void CopyWrapperAttributes( WrapperBase const & source ) override
  {
    WrapperBase::CopyWrapperAttributes( source );
    Wrapper< T > const & castedSource = *cast( &source );
    m_ownsData = castedSource.m_ownsData;
    m_default = castedSource.m_default;
  }


  virtual const std::type_info & get_typeid() const noexcept override final
  {
    return typeid(T);
  }


  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted Wrapper<T>
   */
  static Wrapper< T > * cast( WrapperBase * const base )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< Wrapper< T > * >(base);
#else
    return static_cast< Wrapper< T > * >(base);
#endif
  }

  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted reference to const Wrapper<T>
   */
  static Wrapper< T > const * cast( WrapperBase const * const base )
  {
#ifdef USE_DYNAMIC_CASTING
    return dynamic_cast< Wrapper< T > const * >(base);
#else
    return static_cast< Wrapper< T > const * >(base);
#endif
  }

  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted Wrapper<T>
   */
  static Wrapper< T > & cast( WrapperBase & base )
  {
    if( base.get_typeid() != typeid(T) )
    {
      GEOS_ERROR( "invalid cast attempt" );
    }
    return static_cast< Wrapper< T > & >(base);
  }

  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted reference to const Wrapper<T>
   */
  static Wrapper< T > const & cast( WrapperBase const & base )
  {
    if( base.get_typeid() != typeid(T) )
    {
      GEOS_ERROR( "invalid cast attempt" );
    }
    return static_cast< Wrapper< T > const & >(base);
  }

  /**
   * @brief function to determine if T is packable by the buffer packing functions
   * @return true if T is packable. false if not.
   */
  virtual bool isPackable() const override final
  {
    return bufferOps::is_packable< T >::value;
  }

  /**
   * @brief function to pack T
   * @param buffer the buffer in which to pack T
   * @return number of packed bytes.
   */
  virtual localIndex Pack( char * & buffer ) const override final
  {
    localIndex packedSize = 0;

    packedSize += bufferOps::Pack< true >( buffer, this->getName() );
    packedSize += bufferOps::Pack< true >( buffer, *m_data );

    return packedSize;
  }

  /**
   * @brief function to pack T
   * @param buffer the buffer in which to pack T
   * @param packList indices of T to pack
   * @return number of packed bytes.
   */
  virtual localIndex Pack( char * & buffer, arrayView1d< localIndex const > const & packList ) const override final
  {
    localIndex packedSize = 0;

    static_if( bufferOps::is_packable_by_index< T >::value )
    {
      if( sizedFromParent()==1 )
      {
        packedSize += bufferOps::Pack< true >( buffer, this->getName() );
        packedSize += bufferOps::Pack< true >( buffer, *m_data, packList );
      }
    }
    end_static_if
    return packedSize;
  }

  /**
   * @brief function to pack return the length of packing...without doing the packing.
   * @return size of packed bytes
   */
  virtual localIndex PackSize( ) const override final
  {
    char * buffer = nullptr;
    localIndex packedSize = 0;

    packedSize += bufferOps::Pack< false >( buffer, this->getName() );
    packedSize += bufferOps::Pack< false >( buffer, *m_data );

    return packedSize;
  }

  /**
   * @brief function to get the the packing size
   * @param packList indices of T to pack
   * @return number of packed bytes.
   */
  virtual localIndex PackSize( arrayView1d< localIndex const > const & packList ) const override final
  {

    char * buffer = nullptr;
    localIndex packedSize = 0;

    static_if( bufferOps::is_packable_by_index< T >::value )
    {
      if( sizedFromParent()==1 )
      {
        packedSize += bufferOps::Pack< false >( buffer, this->getName() );
        packedSize += bufferOps::Pack< false >( buffer, *m_data, packList );
      }
    }
    end_static_if

    return packedSize;
  }

  /**
   * @brief function to unpack a buffer into the object referred to by m_data
   * @param buffer
   * @return
   */
  virtual localIndex Unpack( char const * & buffer ) override final
  {
    localIndex unpackedSize = 0;
    string name;
    unpackedSize += bufferOps::Unpack( buffer, name );
    GEOS_ERROR_IF( name != this->getName(), "buffer unpack leads to wrapper names that don't match" );
    unpackedSize += bufferOps::Unpack( buffer, *m_data );
    return unpackedSize;
  }
  virtual localIndex Unpack( char const * & buffer, arrayView1d< localIndex const > const & unpackIndices ) override final
  {
    localIndex unpackedSize = 0;
    static_if( bufferOps::is_packable_by_index< T >::value )
    {
      if( sizedFromParent()==1 )
      {
        string name;
        unpackedSize += bufferOps::Unpack( buffer, name );
        GEOS_ERROR_IF( name != this->getName(), "buffer unpack leads to wrapper names that don't match" );
        unpackedSize += bufferOps::Unpack( buffer, *m_data, unpackIndices );
      }
    }
    end_static_if

    return unpackedSize;
  }



  /// @cond DO_NOT_DOCUMENT
  struct empty_wrapper
  {
    HAS_MEMBER_FUNCTION( empty, bool, const, , )
    template< class U = T >
    static typename std::enable_if< has_memberfunction_empty< U >::value, bool >::type
    empty( Wrapper< T > const * parent )
    {
      return parent->m_data->empty();
    }
    template< class U = T >
    static typename std::enable_if< !has_memberfunction_empty< U >::value, bool >::type
    empty( Wrapper< T > const * parent )
    {
      return parent;
    }
  };/// @endcond DO_NOT_DOCUMENT
  virtual bool empty() const override final
  {
    return empty_wrapper::empty( this );
  }

  /// @cond DO_NOT_DOCUMENT
  struct size_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT( size, 0, int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT( size, 1, unsigned int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT( size, 2, long int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT( size, 3, unsigned long int, const, , )
    HAS_MEMBER_FUNCTION_VARIANT( size, 4, long long, const, , )
    HAS_MEMBER_FUNCTION_VARIANT( size, 5, unsigned long long, const, , )

    template< class U = T >
    static typename std::enable_if< has_memberfunction_v0_size< U >::value ||
                                    has_memberfunction_v1_size< U >::value ||
                                    has_memberfunction_v2_size< U >::value ||
                                    has_memberfunction_v3_size< U >::value ||
                                    has_memberfunction_v4_size< U >::value ||
                                    has_memberfunction_v5_size< U >::value, localIndex >::type
    size( Wrapper< T > const * parent )
    { return static_cast< localIndex >(parent->m_data->size()); }

    template< class U = T >
    static typename std::enable_if< !(has_memberfunction_v0_size< U >::value ||
                                      has_memberfunction_v1_size< U >::value ||
                                      has_memberfunction_v2_size< U >::value ||
                                      has_memberfunction_v3_size< U >::value ||
                                      has_memberfunction_v4_size< U >::value ||
                                      has_memberfunction_v5_size< U >::value), localIndex >::type
    size( Wrapper< T > const * GEOSX_UNUSED_ARG( parent ) )
    { return 1; }
  };/// @endcond DO_NOT_DOCUMENT
  virtual localIndex size() const override final
  {
    return size_wrapper::size( this );
  }


  /// @cond DO_NOT_DOCUMENT
  struct num_dimensions_wrapper
  {
    HAS_MEMBER_FUNCTION( numDimensions, int, const, , )

    template< class U = T >
    static typename std::enable_if< has_memberfunction_numDimensions< U >::value, int >::type
    numDimensions( Wrapper< T > const * parent )
    { return static_cast< int >(parent->m_data->numDimensions()); }

    template< class U = T >
    static typename std::enable_if< !has_memberfunction_numDimensions< U >::value, int >::type
    numDimensions( Wrapper< T > const * GEOSX_UNUSED_ARG( parent ) )
    { return 1; }
  };/// @endcond DO_NOT_DOCUMENT
  virtual int numDimensions() const override final
  {
    return num_dimensions_wrapper::numDimensions( this );
  }

  /// @cond DO_NOT_DOCUMENT
  struct dimension_size_wrapper
  {
    HAS_MEMBER_FUNCTION_VARIANT( size, 0, int, const, VA_LIST( int ), VA_LIST( int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT( size, 1, unsigned int, const, VA_LIST( int ), VA_LIST( int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT( size, 2, long int, const, VA_LIST( int ), VA_LIST( int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT( size, 3, unsigned long int, , VA_LIST( int ), VA_LIST( int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT( size, 4, long long, const, VA_LIST( int ), VA_LIST( int(1)) )
    HAS_MEMBER_FUNCTION_VARIANT( size, 5, unsigned long long, const, VA_LIST( int ), VA_LIST( int(1)) )

    template< class U = T >
    static typename std::enable_if< has_memberfunction_v0_size< U >::value ||
                                    has_memberfunction_v1_size< U >::value ||
                                    has_memberfunction_v2_size< U >::value ||
                                    has_memberfunction_v3_size< U >::value ||
                                    has_memberfunction_v4_size< U >::value ||
                                    has_memberfunction_v5_size< U >::value, localIndex >::type
    size( Wrapper< T > const * const parent, int const i )
    { return integer_conversion< localIndex >( parent->m_data->size( i )); }

    template< class U = T >
    static typename std::enable_if< !(has_memberfunction_v0_size< U >::value ||
                                      has_memberfunction_v1_size< U >::value ||
                                      has_memberfunction_v2_size< U >::value ||
                                      has_memberfunction_v3_size< U >::value ||
                                      has_memberfunction_v4_size< U >::value ||
                                      has_memberfunction_v5_size< U >::value), localIndex >::type
    size( Wrapper< T > const * const parent, int const i )
    {
      if( i != 0 )
      {
        GEOS_ERROR( "Data is only 1D" );
        return 0;
      }
      return parent->size();
    }
  };/// @endcond DO_NOT_DOCUMENT
  virtual localIndex size( int const i ) const override final
  {
    return dimension_size_wrapper::size( this, i );
  }


  /// @cond DO_NOT_DOCUMENT
  struct resize_dimension_wrapper
  {
    HAS_MEMBER_FUNCTION( resize, void, , VA_LIST( int, localIndex const * ),
                         VA_LIST( int(1), static_cast< localIndex const * >(nullptr)))

    template< class U=T >
    static typename std::enable_if< has_memberfunction_resize< U >::value, void >::type
    resize( Wrapper< T > * parent, int num_dims, localIndex const * const dims )
    { parent->m_data->resize( num_dims, dims ); }

    template< class U=T >
    static typename std::enable_if< !has_memberfunction_resize< U >::value, void >::type
    resize( Wrapper< T > * parent, int num_dims, localIndex const * const dims )
    {
      if( num_dims != 1 )
      {
        GEOS_ERROR( "Data is only 1D" );
        return;
      }
      parent->resize( integer_conversion< localIndex >( dims[0] ));
    }
  };/// @endcond DO_NOT_DOCUMENT
  virtual void resize( int num_dims, localIndex const * const dims ) override final
  { resize_dimension_wrapper::resize( this, num_dims, dims ); }

  /// @cond DO_NOT_DOCUMENT
  struct reserve_wrapper
  {
    HAS_MEMBER_FUNCTION( reserve, void, , VA_LIST( std::size_t ), VA_LIST( std::size_t( 1 )) )
    template< class U = T >
    static typename std::enable_if< has_memberfunction_reserve< U >::value, void >::type reserve( Wrapper< T > * const parent, std::size_t new_cap )
    {
      return parent->m_data->reserve( new_cap );
    }
    template< class U = T >
    static typename std::enable_if< !has_memberfunction_reserve< U >::value, void >::type reserve( Wrapper< T > * const, std::size_t )
    {
      return; //parent->m_data;
    }
  };/// @endcond DO_NOT_DOCUMENT
  virtual void reserve( std::size_t new_cap ) override final
  {
    reserve_wrapper::reserve( this, new_cap );
  }
//  CONDITIONAL_VIRTUAL_FUNCTION( Wrapper<T>,reserve , void,,
// VA_LIST(std::size_t a), VA_LIST(a) )


  HAS_MEMBER_FUNCTION( capacity, std::size_t, const, , )
  CONDITIONAL_VIRTUAL_FUNCTION0( Wrapper< T >, capacity, std::size_t, const )

  HAS_MEMBER_FUNCTION( max_size, std::size_t, const, , )
  CONDITIONAL_VIRTUAL_FUNCTION0( Wrapper< T >, max_size, std::size_t, const )

  HAS_MEMBER_FUNCTION( clear, void, , , )
  CONDITIONAL_VIRTUAL_FUNCTION0( Wrapper< T >, clear, void, )

  HAS_MEMBER_FUNCTION( insert, void, , , )
  CONDITIONAL_VIRTUAL_FUNCTION0( Wrapper< T >, insert, void, )


  /// @cond DO_NOT_DOCUMENT
  struct resize_wrapper
  {
    template< typename UU >
    struct has_memberfunction_resize
    {
      HAS_MEMBER_FUNCTION_VARIANT( resize, 0, void, , VA_LIST( int ), VA_LIST( static_cast< int >(1)))
      HAS_MEMBER_FUNCTION_VARIANT( resize, 1, void, , VA_LIST( unsigned int ), VA_LIST( static_cast< unsigned int >(1)))
      HAS_MEMBER_FUNCTION_VARIANT( resize, 2, void, , VA_LIST( long ), VA_LIST( static_cast< long int >(1)))
      HAS_MEMBER_FUNCTION_VARIANT( resize, 3, void, , VA_LIST( unsigned long ), VA_LIST( static_cast< unsigned long int >(1)))
      HAS_MEMBER_FUNCTION_VARIANT( resize, 4, void, , VA_LIST( long long int ), VA_LIST( static_cast< long long int >(1)))
      HAS_MEMBER_FUNCTION_VARIANT( resize, 5, void, , VA_LIST( unsigned long long ), VA_LIST( static_cast< unsigned long long >(1)))

      static constexpr bool value = has_memberfunction_v0_resize< UU >::value ||
                                    has_memberfunction_v1_resize< UU >::value ||
                                    has_memberfunction_v2_resize< UU >::value ||
                                    has_memberfunction_v3_resize< UU >::value ||
                                    has_memberfunction_v4_resize< UU >::value ||
                                    has_memberfunction_v5_resize< UU >::value;
    };

    template< typename UU, typename ENABLE=void >
    struct has_memberfunction_resize2
    {
      static constexpr bool value = false;
    };

    template< typename UU >
    struct has_memberfunction_resize2< UU, typename std::enable_if< DefaultValue< UU >::has_default_value >::type >
    {
      typedef typename DefaultValue< UU >::value_type DVT;
      HAS_MEMBER_FUNCTION_VARIANT( resizeDefault, 0, void, , VA_LIST( int, DVT const & ), VA_LIST( static_cast< int >(1), DVT()))
      HAS_MEMBER_FUNCTION_VARIANT( resizeDefault, 1, void, , VA_LIST( unsigned int, DVT const & ), VA_LIST( static_cast< unsigned int >(1), DVT()))
      HAS_MEMBER_FUNCTION_VARIANT( resizeDefault, 2, void, , VA_LIST( long, DVT const & ), VA_LIST( static_cast< long int >(1), DVT()))
      HAS_MEMBER_FUNCTION_VARIANT( resizeDefault, 3, void, , VA_LIST( unsigned long, DVT const & ), VA_LIST( static_cast< unsigned long int >(1), DVT()))
      HAS_MEMBER_FUNCTION_VARIANT( resizeDefault, 4, void, , VA_LIST( long long int, DVT const & ), VA_LIST( static_cast< long long int >(1), DVT()))
      HAS_MEMBER_FUNCTION_VARIANT( resizeDefault, 5, void, , VA_LIST( unsigned long long, DVT const & ), VA_LIST( static_cast< unsigned long long >(1), DVT()))

      static constexpr bool value = has_memberfunction_v0_resizeDefault< UU >::value ||
                                    has_memberfunction_v1_resizeDefault< UU >::value ||
                                    has_memberfunction_v2_resizeDefault< UU >::value ||
                                    has_memberfunction_v3_resizeDefault< UU >::value ||
                                    has_memberfunction_v4_resizeDefault< UU >::value ||
                                    has_memberfunction_v5_resizeDefault< UU >::value;
    };

    template< class U = T >
    static typename std::enable_if< has_memberfunction_resize2< U >::value &&
                                    DefaultValue< U >::has_default_value, void >::type
    resize( Wrapper * const parent, localIndex const new_size )
    {
      return parent->m_data->resizeDefault( new_size, parent->m_default.value );
    }

    template< class U = T >
    static typename std::enable_if< !(has_memberfunction_resize2< U >::value &&
                                      DefaultValue< U >::has_default_value) &&
                                    has_memberfunction_resize< U >::value, void >::type
    resize( Wrapper * const parent, localIndex const new_size )
    {
      return parent->m_data->resize( new_size );
    }


    template< class U = T >
    static typename std::enable_if< !(has_memberfunction_resize2< U >::value &&
                                      DefaultValue< U >::has_default_value) &&
                                    !has_memberfunction_resize< U >::value, void >::type
    resize( Wrapper * const, localIndex )
    {
      return;
    }
  };/// @endcond DO_NOT_DOCUMENT
  virtual void resize( localIndex new_size ) override final
  {
    resize_wrapper::resize( this, new_size );
  }

  /// @cond DO_NOT_DOCUMENT
  struct should_resize_wrapper
  {
    HAS_MEMBER_FUNCTION( isSorted, bool, const, , )
    template< class U = T >
    static typename std::enable_if< traits::is_set< U >, bool >::type shouldResize()
    { return false;  }

    template< class U = T >
    static typename std::enable_if< !traits::is_set< U >, bool >::type shouldResize()
    { return true; }
  };/// @endcond DO_NOT_DOCUMENT
  virtual bool shouldResize() const override final
  {
    return should_resize_wrapper::shouldResize();
  }

  /// @cond DO_NOT_DOCUMENT
  struct copy_wrapper
  {
    HAS_ALIAS( isArray )

    template< class U=T >
    static typename std::enable_if< has_alias_isArray< U >::value, void >::type
    copy( T * const data, localIndex const sourceIndex, localIndex const destIndex )
    {
      data->copy( destIndex, sourceIndex );
//      (*data)[destIndex] = (*data)[sourceIndex];
    }

    template< class U=T >
    static typename std::enable_if< !has_alias_isArray< U >::value, void >::type
    copy( T * const GEOSX_UNUSED_ARG( data ), localIndex const GEOSX_UNUSED_ARG( sourceIndex ), localIndex const GEOSX_UNUSED_ARG( destIndex ) )
    {}

  };/// @endcond DO_NOT_DOCUMENT
  virtual void copy( localIndex const sourceIndex, localIndex const destIndex ) override final
  {
    if( this->sizedFromParent() )
    {
      copy_wrapper::copy( this->m_data, sourceIndex, destIndex );
    }
  }


  /**
   * @name Structure to determine return types for data access functions
   */
  ///@{

  /// Invoke macro to generate test to see if type has an alias named "pointer".
  /// This will be used to determine if the
  /// type is to be treated as an "array" or a single object.
  HAS_ALIAS( pointer )

  /**
   * SFINAE specialized structure to control return type based on properties of
   * T.
   * The default template returns a pointer for all calls to data().
   */
  template< class U=T,
            bool HASPOINTERTYPE = has_alias_pointer< U >::value >
  struct Get_Type
  {
    /// pointer type
    typedef U * pointer;

    /// pointer to const type
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
  template< class U >
  struct Get_Type< U, true >
  {
    /// pointer type
    typedef typename U::pointer       pointer;

    /// pointer to const type
    typedef typename U::const_pointer const_pointer;
  };

  /// the valid pointer type for T
  using pointer       = typename Get_Type< T >::pointer;

  /// the valid pointer to const type for T
  using const_pointer = typename Get_Type< T >::const_pointer;
  ///@}


  HAS_MEMBER_FUNCTION( data, pointer, , , )
  HAS_MEMBER_FUNCTION_VARIANT( data, _const, pointer, const, , )


  /**
   * @brief accessor for m_data
   * @return reference to T
   */
  T & reference()
  { return *m_data; }

  /**
   * @brief accessor for m_data
   * @return reference to const T
   */
  T const & reference() const
  { return *m_data; }

  /**
   * @brief accessor for m_data
   * @return pointer to T
   */
  T * getPointer()
  { return m_data; }

  /**
   * @brief accessor for m_data
   * @return pointer to const T
   */
  T const * getPointer() const
  { return m_data; }

  HAS_ALIAS( ViewType )

  template< class U=T,
            bool HASPOINTERTYPE = has_alias_ViewType< U >::value >
  struct Get_View_Type
  {
    using ViewType = T &;
    using ViewTypeConst = T const &;
  };

  template< class U >
  struct Get_View_Type< U, true >
  {
    using ViewType = typename T::ViewType const &;
    using ViewTypeConst = typename T::ViewTypeConst const &;
  };

  using ViewType      = typename Get_View_Type< T >::ViewType;

  using ViewTypeConst = typename Get_View_Type< T >::ViewTypeConst;

  template< class U=T >
  typename std::enable_if< has_alias_ViewType< U >::value, ViewType >::type
  referenceAsView()
  { return m_data->toView(); }

  template< class U=T >
  typename std::enable_if< !has_alias_ViewType< U >::value, ViewType >::type
  referenceAsView()
  { return *m_data; }

  template< class U=T >
  typename std::enable_if< has_alias_ViewType< U >::value, ViewTypeConst >::type
  referenceAsView() const
  { return m_data->toViewConst(); }

  template< class U=T >
  typename std::enable_if< !has_alias_ViewType< U >::value, ViewType >::type
  referenceAsView() const
  { return *m_data; }

  /**
   * @brief accessor for m_default
   * @return reference to const m_default member
   */
  template< typename U=T >
  DefaultValue< T > const &
  getDefaultValueStruct() const
  {
    return m_default;
  }


  /**
   * @brief accessor for default value
   * @return reference to const T
   */
  template< typename U=T >
  typename std::enable_if< DefaultValue< U >::has_default_value, T const & >::type
  getDefaultValue() const
  {
    return m_default.value;
  }

  /**
   * @brief setter for default value
   * @return pointer to Wrapper<T>
   */
  template< typename U=T >
  typename std::enable_if< DefaultValue< U >::has_default_value, Wrapper< T > * >::type
  setDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    return this;
  }

  /**
   * @brief set and apply for default value
   * @return pointer to Wrapper<T>
   */
  template< typename U=T >
  typename std::enable_if< DefaultValue< U >::has_default_value, Wrapper< T > * >::type
  setApplyDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    *m_data = m_default.value;
    return this;
  }

  /**
   * @brief Case for if m_data has a member function called "data()"
   * @return pointer type specified by T::pointer
   */
  template< class U = T >
  typename std::enable_if< ( has_memberfunction_data< U >::value || has_memberfunction_v_const_data< U >::value ) &&
                           has_alias_pointer< U >::value && !std::is_same< U, string >::value, typename U::pointer >::type
  dataPtr()
  {
    return m_data->data();
  }

  /**
   * @brief Case for if m_data has a member function called "data()"
   * @return pointer type specified by T::const_pointer
   */
  template< class U = T >
  typename std::enable_if< ( has_memberfunction_data< U >::value || has_memberfunction_v_const_data< U >::value ) &&
                           has_alias_pointer< U >::value && !std::is_same< U, string >::value, typename U::const_pointer >::type
  dataPtr() const
  {
    return m_data->data();
  }


  /**
   * @brief  Case for if m_data is a string"
   * @return pointer type specified by char *
   */
  template< class U = T >
  typename std::enable_if< std::is_same< U, string >::value, char * >::type
  dataPtr()
  {
    return const_cast< char * >(m_data->data());
  }

  /**
   * @brief  Case for if m_data is a string"
   * @return pointer type specified by char const *
   */
  template< class U = T >
  typename std::enable_if< std::is_same< U, string >::value, char const * >::type
  dataPtr() const
  {
    return m_data->data();
  }


  /**
   * @brief  case for if m_data does NOT have a member function "data()"
   * @return pointer type specified by T *
   */
  template< class U = T >
  typename std::enable_if< !( has_memberfunction_data< U >::value || has_memberfunction_v_const_data< U >::value )&&
                           !std::is_same< U, string >::value, U * >::type
  dataPtr()
  {
    return m_data;
  }

  /**
   * @brief  case for if m_data does NOT have a member function "data()"
   * @return pointer type specified by T *
   */
  template< class U = T >
  typename std::enable_if< !( has_memberfunction_data< U >::value || has_memberfunction_v_const_data< U >::value )&&
                           !std::is_same< U, string >::value, U const * >::type
  dataPtr() const
  {
    return m_data;
  }

  HAS_MEMBER_FUNCTION( setUserCallBack,
                       void,
                       ,
                       std::string const &,
                       "" )

  template< class U = T >
  typename std::enable_if< has_memberfunction_setUserCallBack< U >::value, void >::type
  setUserCallBack()
  {
    std::string const path = getSidreView()->getPathName();
    m_data->setUserCallBack( path );
  }

  template< class U = T >
  typename std::enable_if< !has_memberfunction_setUserCallBack< U >::value, void >::type
  setUserCallBack()
  {}



  struct move_wrapper
  {
    HAS_MEMBER_FUNCTION( move,
                         void,
                         ,
                         VA_LIST( chai::ExecutionSpace, bool ),
                         VA_LIST( chai::CPU, true ) )

    template< class U = T >
    static typename std::enable_if< has_memberfunction_move< U >::value, void >::type
    move( U & data, chai::ExecutionSpace space, bool touch )
    {
      data.move( space, touch );
    }

    template< class U = T >
    static typename std::enable_if< !has_memberfunction_move< U >::value, void >::type
    move( U &, chai::ExecutionSpace, bool )
    {}
  };

  virtual void move( chai::ExecutionSpace space, bool touch ) override
  { return move_wrapper::move( *m_data, space, touch ); }

  HAS_ALIAS( value_type )


  /**
   * @brief function to get the size of T
   * @return size of T
   */
  template< class U = T >
  typename std::enable_if< has_alias_value_type< U >::value, size_t >::type
  sizeOfValueType() const
  {
    return sizeof(typename T::value_type);
  }

  /**
   * @brief function to get the size of T
   * @return size of T
   */
  template< class U = T >
  typename std::enable_if< !has_alias_value_type< U >::value, size_t >::type
  sizeOfValueType() const
  {
    return sizeof(T);
  }

  virtual size_t sizeOfType() const override final
  {
    return sizeOfValueType();
  }


  /**
   * @brief case for if U::value_type exists. Returns the size of dataPtr
   * @return size of T::value_type
   */
  template< class U = T >
  typename std::enable_if< has_alias_value_type< U >::value, localIndex >::type
  byteSize() const
  {
    return size() * sizeof(typename T::value_type);
  }


  /**
   * @brief case for if U::value_type doesn't exists. Returns the size of dataPtr
   * @return size of T::value_type
   */
  template< class U = T >
  typename std::enable_if< !has_alias_value_type< U >::value, localIndex >::type
  byteSize() const
  {
    return size() * sizeof(T);
  }


  /**
   * @brief case for if U::value_type exists. Returns the size of an element of dataPtr
   * @return size of T::value_type
   */
  template< class U = T >
  typename std::enable_if< has_alias_value_type< U >::value, localIndex >::type
  elementSize() const
  {
    return sizeof(typename T::value_type);
  }

  /**
   * @brief case for if U::value_type doesn't exists. Returns the size of an element of dataPtr
   * @return size of T::value_type
   */
  template< class U = T >
  typename std::enable_if< !has_alias_value_type< U >::value, localIndex >::type
  elementSize() const
  {
    return sizeof(T);
  }


  /// case for if U::value_type exists. Returns the typeid of an element of dataPtr
  template< class U = T >
  typename std::enable_if< has_alias_value_type< U >::value, const std::type_info & >::type
  elementTypeID() const
  {
    return typeid(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the typeid of an element of dataPtr
  template< class U = T >
  typename std::enable_if< !has_alias_value_type< U >::value, const std::type_info & >::type
  elementTypeID() const
  {
    return typeid(T);
  }


  /// case for if U::value_type exists. Returns the number of elements given a byte size
  template< class U = T >
  typename std::enable_if< has_alias_value_type< U >::value, localIndex >::type
  numElementsFromByteSize( localIndex d_size ) const
  {
    return d_size / sizeof(typename T::value_type);
  }


  /// case for if U::value_type doesn't exists. Returns the number of elements
  /// given a byte size
  template< class U = T >
  typename std::enable_if< !has_alias_value_type< U >::value, localIndex >::type
  numElementsFromByteSize( localIndex d_size ) const
  {
    return d_size / sizeof(T);
  }

  /// @cond DO_NOT_DOCUMENT

  virtual bool shouldRegisterDataPtr() const override
  {
    std::type_index type_index = std::type_index( elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType( type_index );
    return sidre_type_id != axom::sidre::TypeID::NO_TYPE_ID;
  }



  void registerDataPtr( axom::sidre::View * view ) const override
  {
#ifdef GEOSX_USE_ATK
    view = (view != nullptr) ? view : getSidreView();

    localIndex num_elements = size();
    if( num_elements > 0 )
    {
      std::type_index type_index = std::type_index( elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType( type_index );
      if( sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID )
      { return; }

      localIndex sidre_size = rtTypes::getSidreSize( type_index );
      localIndex byte_size = byteSize();
      localIndex element_size = elementSize();

      int ndims = numDimensions();
      axom::sidre::IndexType dims[10];
      for( int dim = 0 ; dim < ndims ; ++dim )
      {
        dims[dim] = size( dim );
      }

      if( byte_size > num_elements * sidre_size )
      {
        dims[ndims++] = element_size / sidre_size;
      }

      void * ptr = const_cast< void * >( static_cast< void const * >( dataPtr() ) );
      view->setExternalDataPtr( sidre_type_id, ndims, dims, ptr );
    }
    else
    {
      unregisterDataPtr( view );
    }
#endif
  }

  /* Register the pointer to data with the associated sidre::View. */
  void registerToWrite( axom::sidre::View * view=nullptr ) override
  {
#ifdef GEOSX_USE_ATK
    if( getRestartFlags() == RestartFlags::NO_WRITE )
    {
      view = (view != nullptr) ? view : getSidreView();
      unregisterDataPtr( view );
      return;
    }

    move( chai::CPU, false );

    view = (view != nullptr) ? view : getSidreView();
    storeSizedFromParent( view );

    localIndex num_elements = size();
    if( num_elements > 0 )
    {
      std::type_index type_index = std::type_index( elementTypeID());
      axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType( type_index );
      if( sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID )
      {
        localIndex byte_size = bufferOps::PackSize( reference());
        char * const buffer = new char[byte_size];
        char * buffer_cpy = buffer;
        bufferOps::Pack< true >( buffer_cpy, reference());
        view->setExternalDataPtr( axom::sidre::TypeID::INT8_ID, byte_size, buffer );
        return;
      }

      localIndex sidre_size = rtTypes::getSidreSize( type_index );
      localIndex byte_size = byteSize();
      localIndex element_size = elementSize();

      int ndims = numDimensions();
      axom::sidre::IndexType dims[10];
      for( int dim = 0 ; dim < ndims ; ++dim )
      {
        dims[dim] = size( dim );
      }

      if( byte_size > num_elements * sidre_size )
      {
        dims[ndims++] = element_size / sidre_size;
      }

      void * ptr = const_cast< void * >(static_cast< void const * >( dataPtr() ) );
      view->setExternalDataPtr( sidre_type_id, ndims, dims, ptr );
    }
    else
    {
      unregisterDataPtr( view );
    }
#endif
  }

  /* Register the pointer to data with the associated sidre::View. */
  void finishWriting( axom::sidre::View * view=nullptr ) const override
  {
#ifdef GEOSX_USE_ATK
    if( getRestartFlags() == RestartFlags::NO_WRITE )
    {
      view = (view != nullptr) ? view : getSidreView();
      unregisterDataPtr( view );
      return;
    }

    view = (view != nullptr) ? view : getSidreView();
    view->setAttributeToDefault( "__sizedFromParent__" );

    if( !view->isExternal() || view->getTotalBytes() == 0 )
    {
      return;
    }

    std::type_index type_index = std::type_index( elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType( type_index );
    if( sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID )
    {
      delete[] static_cast< char * >(view->getVoidPtr());
    }

    unregisterDataPtr( view );
#endif
  }

  void registerToRead( axom::sidre::View * view=nullptr ) override
  {
#ifdef GEOSX_USE_ATK
    if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
    {
      unregisterDataPtr( view );
      return;
    }

    view = (view != nullptr) ? view : getSidreView();
    loadSizedFromParent( view );
    if( !view->isExternal() || view->getTotalBytes() == 0 )
    {
      return;
    }

    std::type_index type_index = std::type_index( elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType( type_index );
    if( sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID )
    {
      localIndex byte_size = integer_conversion< localIndex >( view->getTotalBytes());
      char * ptr = new char[byte_size];
      view->setExternalDataPtr( axom::sidre::TypeID::INT8_ID, byte_size, ptr );
      return;
    }

    resizeFromSidre( view );
    void * ptr = const_cast< void * >( static_cast< void const * >( dataPtr() ) );
    localIndex sidre_size = rtTypes::getSidreSize( type_index );
    view->setExternalDataPtr( sidre_type_id, byteSize() / sidre_size, ptr );
#endif
  }


  void finishReading( axom::sidre::View * view ) override
  {
#ifdef GEOSX_USE_ATK
    if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
    {
      view = (view != nullptr) ? view : getSidreView();
      unregisterDataPtr( view );
      return;
    }
    view = (view != nullptr) ? view : getSidreView();
    if( !view->isExternal() || view->getTotalBytes() == 0 )
    {
      return;
    }

    std::type_index type_index = std::type_index( elementTypeID());
    axom::sidre::TypeID sidre_type_id = rtTypes::toSidreType( type_index );
    if( sidre_type_id == axom::sidre::TypeID::NO_TYPE_ID )
    {
      localIndex const byte_size = integer_conversion< localIndex >( view->getTotalBytes());
      const char * const buffer = static_cast< char * >(view->getVoidPtr());
      const char * buffer_cpy = buffer;
      localIndex const bytes_read = bufferOps::Unpack( buffer_cpy, reference());
      GEOS_ERROR_IF( bytes_read != byte_size, bytes_read << " != " << byte_size );
      delete[] buffer;
    }

    unregisterDataPtr( view );
#endif
  }

  void unregisterDataPtr( axom::sidre::View * view = nullptr ) const
  {
#ifdef GEOSX_USE_ATK
    view = (view != nullptr) ? view : getSidreView();
    view->setExternalDataPtr( nullptr );
#endif
  }

  void storeSizedFromParent( axom::sidre::View * view = nullptr ) const
  {
#ifdef GEOSX_USE_ATK
    if( SidreWrapper::dataStore().hasAttribute( "__sizedFromParent__" ))
    {
      view = (view != nullptr) ? view : getSidreView();
      view->setAttributeScalar( "__sizedFromParent__", sizedFromParent());
    }
#endif
  }

  void loadSizedFromParent( axom::sidre::View * view = nullptr )
  {
#ifdef GEOSX_USE_ATK
    if( SidreWrapper::dataStore().hasAttribute( "__sizedFromParent__" ))
    {
      view = (view != nullptr) ? view : getSidreView();
      setSizedFromParent( view->getAttributeScalar( "__sizedFromParent__" ));
      view->setAttributeToDefault( "__sizedFromParent__" );
    }
#endif
  }

  /**
   *
   * @param view
   */
  void resizeFromSidre( axom::sidre::View * view = nullptr )
  {
#ifdef GEOSX_USE_ATK
    view = (view != nullptr) ? view : getSidreView();
    if( view->isExternal())
    {
      std::type_index type_index = std::type_index( elementTypeID());
      localIndex sidre_size = rtTypes::getSidreSize( type_index );

      localIndex byte_size = integer_conversion< localIndex >( view->getTotalBytes());
      localIndex num_elements = numElementsFromByteSize( byte_size );

      int ndims = view->getNumDimensions();
      axom::sidre::IndexType dims[10];
      view->getShape( ndims, dims );

      if( byte_size > num_elements * sidre_size )
      {
        ndims--;
      }

      localIndex num_elems_recorded = 1;
      for( localIndex i = 0 ; i < ndims ; ++i )
      {
        num_elems_recorded *= dims[i];
      }

      if( num_elems_recorded != num_elements )
      {
        GEOS_ERROR( "Number of elements recorded not equal to the calculated number: " <<
                    num_elems_recorded << " " << num_elements );
      }

//      long long l_dims[ndims];
      localIndex l_dims[10];
      for( localIndex i = 0 ; i < ndims ; ++i )
      {
        l_dims[i] = dims[i];
      }

      resize( ndims, l_dims );
    }
#endif
  }
  /// @endcond DO_NOT_DOCUMENT


  /** @name overridden setters
   *  Group of setters that override non-virtual functions in WrapperBase
   */
  ///@{

  /**
   * @brief set whether this wrapper is resized when its parent is resized
   * @param val an int that is converted into a bool
   * @return a pointer to this wrapper
   */
  Wrapper< T > * setSizedFromParent( int val )
  {
    WrapperBase::setSizedFromParent( val );
    return this;
  }

  /**
   * @brief set the RestartFlags of the wrapper
   * @param flags the new RestartFlags value
   * @return a pointer to this wrapper
   */
  Wrapper< T > * setRestartFlags( RestartFlags flags )
  {
    WrapperBase::setRestartFlags( flags );
    return this;
  }

  /**
   * @brief set the PlotLevel of the wrapper
   * @param flag the new PlotLevel value
   * @return a pointer to this wrapper
   */
  Wrapper< T > * setPlotLevel( PlotLevel const flag )
  {
    WrapperBase::setPlotLevel( flag );
    return this;
  }

  /**
   * @brief set the plotLevel of the wrapper
   * @param flag an integer that specifies the new plotLevel value
   * @return a pointer to this wrapper
   */
  Wrapper< T > * setPlotLevel( int const flag )
  {
    WrapperBase::setPlotLevel( flag );
    return this;
  }

  /**
   * @brief set the InputFlag of the wrapper
   * @param input the new InputFlags value
   * @return a pointer to this wrapper
   */
  Wrapper< T > * setInputFlag( InputFlags const input )
  {
    WrapperBase::setInputFlag( input );
    return this;
  }

  /**
   * @brief set the description string of the wrapper
   * @param description the description
   * @return a pointer to this wrapper
   */
  Wrapper< T > * setDescription( string const & description )
  {
    WrapperBase::setDescription( description );
    return this;
  }


  ///@}

#if defined(USE_TOTALVIEW_OUTPUT)
  virtual string totalviewTypeName() const override final
  {
    return cxx_utilities::demangle( typeid( Wrapper< T > ).name() );
  }

  virtual int setTotalviewDisplay() const override final
  {
    //std::cout<<"executing Wrapper::setTotalviewDisplay()"<<std::endl;
    WrapperBase::setTotalviewDisplay();
    TV_ttf_add_row( "m_ownsData", "bool", &m_ownsData );
    TV_ttf_add_row( "m_data", totalview::typeName< T >().c_str(), m_data );
    TV_ttf_add_row( "m_default", totalview::typeName< DefaultValue< T > >().c_str(), &m_default );
    return 0;
  }
//  void tvTemplateInstantiation();
#endif

private:
  /// flag to indicate whether or not this wrapper is responsible for allocation/deallocation of the object at the
  /// address of m_data
  bool m_ownsData;

  /// the object being wrapped by this wrapper
  T * m_data;

  /// the default value of the object being wrapped
  DefaultValue< T > m_default;


  Wrapper() = delete;
};

}
} /* namespace geosx */

//template< typename T >
//int TV_ttf_display_type( geosx::dataRepository::Wrapper<T> const * wrapper)
//{
//  std::cout<<"Executing "<<wrapper->totalviewTypeName()<<"::TV_ttf_display_type()"<<std::endl;
//  return TV_ttf_format_raw;
//}
//
//template int TV_ttf_display_type( geosx::dataRepository::Wrapper<int> const * wrapper );
//
//template< typename T >
//void geosx::dataRepository::Wrapper<T>::tvTemplateInstantiation()
//{
//  TV_ttf_display_type<T>(this);
//}


#endif /* GEOSX_DATAREPOSITORY_WRAPPER_HPP_ */
