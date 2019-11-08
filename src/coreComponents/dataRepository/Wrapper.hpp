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

// Source inclues
#include "WrapperHelpers.hpp"
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

// System includes
#include <type_traits>
#include <cstdlib>
#include <type_traits>


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
    return dynamicCast< Wrapper< T > * >( base );
  }

  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted reference to const Wrapper<T>
   */
  static Wrapper< T > const * cast( WrapperBase const * const base )
  {
    return dynamicCast< Wrapper< T > const * >( base );
  }

  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted Wrapper<T>
   */
  static Wrapper< T > & cast( WrapperBase & base )
  {
    return dynamicCast< Wrapper< T > & >( base );
  }

  /**
   * static function to cast a Wrapper base to a derived Wrapper<T>
   * @param base
   * @return casted reference to const Wrapper<T>
   */
  static Wrapper< T > const & cast( WrapperBase const & base )
  {
    return dynamicCast< Wrapper< T > const & >( base );
  }

  /**
   * @brief function to determine if T is packable by the buffer packing functions
   * @return true if T is packable. false if not.
   */
  virtual bool isPackable() const override final
  {
    return bufferOps::is_packable< T >;
  }

  /**
   * @brief function to pack T
   * @param buffer the buffer in which to pack T
   * @return number of packed bytes.
   */
  virtual localIndex Pack( buffer_unit_type * & buffer ) const override final
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
  virtual localIndex Pack( buffer_unit_type * & buffer, arrayView1d< localIndex const > const & packList ) const override final
  {
    localIndex packedSize = 0;

    static_if( bufferOps::is_packable_by_index< T > )
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
    buffer_unit_type * buffer = nullptr;
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

    buffer_unit_type * buffer = nullptr;
    localIndex packedSize = 0;

    static_if( bufferOps::is_packable_by_index< T > )
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
  virtual localIndex Unpack( buffer_unit_type const * & buffer ) override final
  {
    localIndex unpackedSize = 0;
    string name;
    unpackedSize += bufferOps::Unpack( buffer, name );
    GEOS_ERROR_IF( name != this->getName(), "buffer unpack leads to wrapper names that don't match" );
    unpackedSize += bufferOps::Unpack( buffer, *m_data );
    return unpackedSize;
  }

  virtual localIndex Unpack( buffer_unit_type const * & buffer, arrayView1d< localIndex const > const & unpackIndices ) override final
  {
    localIndex unpackedSize = 0;
    static_if( bufferOps::is_packable_by_index< T > )
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

  virtual localIndex size() const override final
  {
    return wrapperHelpers::size( *m_data );
  }

  virtual void resize( int ndims, localIndex const * const dims ) override final
  { wrapperHelpers::resizeDimensions( *m_data, ndims, dims ); }

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

  HAS_MEMBER_FUNCTION( capacity, std::size_t, const, , )
  CONDITIONAL_VIRTUAL_FUNCTION0( Wrapper< T >, capacity, std::size_t, const )

  virtual void resize( localIndex const newSize ) override final
  {
    wrapperHelpers::resizeDefault( *m_data, newSize, m_default );
  }

  /// @cond DO_NOT_DOCUMENT
  struct copy_wrapper
  {
    template< class U=T >
    static typename std::enable_if< traits::is_array< U >, void >::type
    copy( T * const data, localIndex const sourceIndex, localIndex const destIndex )
    {
      data->copy( destIndex, sourceIndex );
    }

    template< class U=T >
    static typename std::enable_if< !traits::is_array< U >, void >::type
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
  typename std::enable_if< DefaultValue< U >::has_default_value, typename DefaultValue< U >::value_type const & >::type
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


  traits::Pointer< T > dataPtr()
  {
    return wrapperHelpers::dataPtr( *m_data );
  }

  traits::ConstPointer< T > dataPtr() const
  {
    return wrapperHelpers::dataPtr( *m_data );
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
    std::string const path = getConduitNode().path();
    m_data->setUserCallBack( path );
  }

  template< class U = T >
  typename std::enable_if< !has_memberfunction_setUserCallBack< U >::value, void >::type
  setUserCallBack()
  {}



  struct move_wrapper
  {
    template< class U = T >
    static typename std::enable_if< traits::has_chai_move_method< U >, void >::type
    move( U & data, chai::ExecutionSpace space, bool touch )
    { data.move( space, touch ); }

    template< class U = T >
    static typename std::enable_if< !traits::has_chai_move_method< U >, void >::type
    move( U &, chai::ExecutionSpace, bool )
    {}
  };

  virtual void move( chai::ExecutionSpace space, bool touch ) override
  { return move_wrapper::move( *m_data, space, touch ); }

  /// @cond DO_NOT_DOCUMENT

  /* Register the pointer to data with the associated conduit::Node. */
  void registerToWrite() override
  {
    getConduitNode().reset();

    if( getRestartFlags() == RestartFlags::NO_WRITE )
    {
      return;
    }

    move( chai::CPU, false );

    getConduitNode()[ "__sizedFromParent__" ].set( sizedFromParent() );

    wrapperHelpers::pushDataToConduitNode( *m_data, getConduitNode() );
  }

  /* Register the pointer to data with the associated conduit::Node. */
  void finishWriting() override
  {
    getConduitNode().reset();
  }

  void loadFromConduit() override
  {
    if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
    {
      getConduitNode().reset();
      return;
    }

    setSizedFromParent( getConduitNode()["__sizedFromParent__"].value() );

    wrapperHelpers::pullDataFromConduitNode( *m_data, getConduitNode() );

    getConduitNode().reset();
  }


  /// @endcond DO_NOT_DOCUMENT


  /**
   *  @name overridden setters
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
