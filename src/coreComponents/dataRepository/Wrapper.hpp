/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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
#include "wrapperHelpers.hpp"
#include "KeyNames.hpp"
#include "LvArray/src/limits.hpp"
#include "common/DataTypes.hpp"
#include "codingUtilities/SFINAE_Macros.hpp"
#include "LvArray/src/Macros.hpp"
#include "BufferOps.hpp"
#include "BufferOpsDevice.hpp"
#include "RestartFlags.hpp"
#include "codingUtilities/traits.hpp"
#include "common/GeosxConfig.hpp"
#include "DefaultValue.hpp"
#include "LvArray/src/system.hpp"
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
 * @tparam T is any type that is to be wrapped by Wrapper
 */
template< typename T >
class Wrapper final : public WrapperBase
{
public:

  /**
   * @brief Alias for the wrapped type @p T
   */
  using TYPE = T;

  /**
   * @name Constructors, destructor, copy/move assignment
   */
  ///@{

  /**
   * @brief Constructor that creates a new instance of wrapped type
   * @param name name of the object
   * @param parent parent group which owns the Wrapper
   */
  explicit Wrapper( string const & name,
                    Group & parent ):
    WrapperBase( name, parent ),
    m_ownsData( true ),
    m_data( new T() ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value || traits::is_string< T > )
    {
      setSizedFromParent( 0 );
    }

    setName();
  }

  /**
   * @brief Constructor that takes ownership of an existing object.
   * @param name name of the object
   * @param parent parent group that owns the Wrapper
   * @param object object that is being wrapped by the Wrapper
   */
  explicit Wrapper( string const & name,
                    Group & parent,
                    std::unique_ptr< T > object ):
    WrapperBase( name, parent ),
    m_ownsData( true ),
    m_data( object.release() ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value || traits::is_string< T > )
    {
      setSizedFromParent( 0 );
    }

    setName();
  }

  /**
   * @brief Constructor that does not take ownership of an existing object.
   * @param name name of the object
   * @param parent parent group that owns the Wrapper
   * @param object object that is being wrapped by the Wrapper
   */
  explicit Wrapper( string const & name,
                    Group & parent,
                    T * object ):
    WrapperBase( name, parent ),
    m_ownsData( false ),
    m_data( object ),
    m_default()
  {
    if( traits::is_tensorT< T > || std::is_arithmetic< T >::value || traits::is_string< T > )
    {
      setSizedFromParent( 0 );
    }

    setName();
  }

  /**
   * @brief Default destructor
   *
   * Deletes wrapped object if the wrapper is owning
   */
  virtual ~Wrapper() noexcept override
  {
    if( m_ownsData )
    {
      delete m_data;
    }
    //tvTemplateInstantiation();
  }

  /**
   * @brief Copy Assignment Operator
   *
   * @param source rhs
   * @return @p *this
   */
  Wrapper & operator=( Wrapper const & source )
  {
    m_data = source.m_data;
    return *this;
  }

  /**
   * @brief Move Assignment Operator
   * @param source
   * @return *this
   */
  Wrapper & operator=( Wrapper && source )
  {
    m_data = std::move( source.m_data );
    return *this;
  }

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual std::unique_ptr< WrapperBase > clone( string const & name,
                                                Group & parent ) override
  {
    std::unique_ptr< WrapperBase >
    clonedWrapper = std::make_unique< Wrapper< T > >( name, parent, m_data );
    clonedWrapper->copyWrapperAttributes( *this );

    return clonedWrapper;
  }

  virtual void copyWrapper( WrapperBase const & source ) override
  {
    GEOSX_ERROR_IF( source.getName() != m_name, "Tried to clone wrapper of with different name" );
    WrapperBase::copyWrapperAttributes( source );
    Wrapper< T > const & castedSource = dynamicCast< Wrapper< T > const & >( source );
    m_ownsData = castedSource.m_ownsData;
    m_default = castedSource.m_default;
    copyData( source );

  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void copyWrapperAttributes( WrapperBase const & source ) override
  {
    WrapperBase::copyWrapperAttributes( source );
    Wrapper< T > const & castedSource = dynamicCast< Wrapper< T > const & >( source );
    m_ownsData = castedSource.m_ownsData;
    m_default = castedSource.m_default;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual const std::type_info & getTypeId() const noexcept override
  {
    return typeid(T);
  }

  ///@}

  /// @copydoc WrapperBase::getHistoryMetadata
  virtual
  HistoryMetadata getHistoryMetadata( localIndex const packCount = -1 ) const override final
  {
    return geosx::getHistoryMetadata( getName(), referenceAsView( ), packCount );
  }

  /**
   * @name Packing/unpacking methods.
   */
  ///@{

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::isPackable
  virtual
  bool isPackable( bool onDevice ) const override
  {
    if( onDevice )
    {
      // this isn't accurate if array/arraview return false for this, which I think they do
      return bufferOps::can_memcpy< T >;
    }
    else
    {
      return bufferOps::is_packable< T >;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::pack
  virtual
  localIndex pack( buffer_unit_type * & buffer, bool withMetadata, bool onDevice, parallelDeviceEvents & events ) const override final
  {
    localIndex packedSize = 0;
    if( withMetadata ) packedSize += bufferOps::Pack< true >( buffer, getName() );
    if( onDevice )
    {
      if( withMetadata )
      {
        packedSize += wrapperHelpers::PackDevice< true >( buffer, reference(), events );
      }
      else
      {
        packedSize += wrapperHelpers::PackDataDevice< true >( buffer, reference(), events );
      }
    }
    else
    {
      packedSize += bufferOps::Pack< true >( buffer, *m_data );
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::packByIndex
  virtual
  localIndex packByIndex( buffer_unit_type * & buffer, arrayView1d< localIndex const > const & packList, bool withMetadata, bool onDevice, parallelDeviceEvents & events ) const override final
  {
    localIndex packedSize = 0;
    if( sizedFromParent() == 1 )
    {
      if( withMetadata ) packedSize += bufferOps::Pack< true >( buffer, getName() );
      if( onDevice )
      {
        if( withMetadata )
        {
          packedSize += wrapperHelpers::PackByIndexDevice< true >( buffer, reference(), packList, events );
        }
        else
        {
          packedSize += wrapperHelpers::PackDataByIndexDevice< true >( buffer, reference(), packList, events );
        }
      }
      else
      {
        packedSize += wrapperHelpers::PackByIndex< true >( buffer, *m_data, packList );
      }
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::packSize
  virtual
  localIndex packSize( bool withMetadata, bool onDevice, parallelDeviceEvents & events ) const override final
  {
    buffer_unit_type * buffer = nullptr;
    localIndex packedSize = 0;
    if( withMetadata ) packedSize += bufferOps::Pack< false >( buffer, getName() );
    if( onDevice )
    {
      if( withMetadata )
      {
        packedSize += wrapperHelpers::PackDevice< false >( buffer, reference(), events );
      }
      else
      {
        packedSize += wrapperHelpers::PackDataDevice< false >( buffer, reference(), events );
      }
    }
    else
    {
      packedSize += bufferOps::Pack< false >( buffer, *m_data );
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::packByIndexSize
  virtual
  localIndex packByIndexSize( arrayView1d< localIndex const > const & packList, bool withMetadata, bool onDevice, parallelDeviceEvents & events ) const override final
  {
    localIndex packedSize = 0;
    buffer_unit_type * buffer = nullptr;
    if( sizedFromParent() == 1 )
    {
      if( withMetadata ) packedSize += bufferOps::Pack< false >( buffer, getName() );
      if( onDevice )
      {
        if( withMetadata )
        {
          packedSize += wrapperHelpers::PackByIndexDevice< false >( buffer, reference(), packList, events );
        }
        else
        {
          packedSize += wrapperHelpers::PackDataByIndexDevice< false >( buffer, reference(), packList, events );
        }
      }
      else
      {
        packedSize += wrapperHelpers::PackByIndex< false >( buffer, *m_data, packList );
      }
    }
    return packedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::unpack
  virtual
  localIndex unpack( buffer_unit_type const * & buffer, bool withMetadata, bool onDevice, parallelDeviceEvents & events ) override final
  {
    localIndex unpackedSize = 0;
    if( withMetadata )
    {
      string name;
      unpackedSize += bufferOps::Unpack( buffer, name );
      GEOSX_ERROR_IF( name != getName(), "buffer unpack leads to wrapper names that don't match" );
    }
    if( onDevice )
    {
      if( withMetadata )
      {
        unpackedSize += wrapperHelpers::UnpackDevice( buffer, referenceAsView(), events );
      }
      else
      {
        unpackedSize += wrapperHelpers::UnpackDataDevice( buffer, referenceAsView(), events );
      }
    }
    else
    {
      unpackedSize += bufferOps::Unpack( buffer, *m_data );
    }
    return unpackedSize;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// @copydoc WrapperBase::unpackByIndex
  virtual
  localIndex unpackByIndex( buffer_unit_type const * & buffer, arrayView1d< localIndex const > const & unpackIndices, bool withMetadata, bool onDevice, parallelDeviceEvents & events ) override final
  {
    localIndex unpackedSize = 0;
    if( sizedFromParent()==1 )
    {
      if( withMetadata )
      {
        string name;
        unpackedSize += bufferOps::Unpack( buffer, name );
        GEOSX_ERROR_IF( name != getName(), "buffer unpack leads to wrapper names that don't match" );
      }
      if( onDevice )
      {
        if( withMetadata )
        {
          unpackedSize += wrapperHelpers::UnpackByIndexDevice( buffer, referenceAsView(), unpackIndices, events );
        }
        else
        {
          unpackedSize += wrapperHelpers::UnpackDataByIndexDevice( buffer, referenceAsView(), unpackIndices, events );
        }
      }
      else
      {
        unpackedSize += wrapperHelpers::UnpackByIndex( buffer, *m_data, unpackIndices );
      }
    }
    return unpackedSize;
  }

  ///@}

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void const * voidPointer() const override
  { return wrapperHelpers::dataPtr( *m_data ); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex elementByteSize() const override
  { return wrapperHelpers::byteSizeOfElement< T >(); }

  /**
   * @name Methods that delegate to the wrapped type
   *
   * These functions will call the corresponding method on the wrapped
   * object, if such method is declared in wrapped type.
   */
  ///@{

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex size() const override
  { return wrapperHelpers::size( *m_data ); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void resize( int ndims, localIndex const * const dims ) override
  {
    wrapperHelpers::move( *m_data, LvArray::MemorySpace::CPU, true );
    wrapperHelpers::resizeDimensions( *m_data, ndims, dims );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void reserve( localIndex const newCapacity ) override
  {
    wrapperHelpers::move( *m_data, LvArray::MemorySpace::CPU, true );
    wrapperHelpers::reserve( reference(), newCapacity );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual localIndex capacity() const override
  {
    // We don't use reference() here because that would return an ArrayView which has no capacity method.
    return wrapperHelpers::capacity( *m_data );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void resize( localIndex const newSize ) override
  {
    wrapperHelpers::move( *m_data, LvArray::MemorySpace::CPU, true );
    wrapperHelpers::resizeDefault( reference(), newSize, m_default );
  }

  /// @cond DO_NOT_DOCUMENT
  struct copy_wrapper
  {
    template< typename U, int NDIM, typename PERMUTATION >
    static void copy( Array< U, NDIM, PERMUTATION > const & array, localIndex const sourceIndex, localIndex const destIndex )
    {
      LvArray::forValuesInSliceWithIndices( array[ sourceIndex ],
                                            [destIndex, &array]( U const & sourceVal, auto const ... indices )
      {
        array( destIndex, indices ... ) = sourceVal;
      } );
    }

    template< typename U >
    static void copy( U const &, localIndex const, localIndex const )
    {}

    template< typename U=T >
    static std::enable_if_t< traits::hasCopyAssignmentOp< U > >
    copyData( U & destinationData, U const & sourceData )
    {
      destinationData = sourceData;
    }

    template< typename U=T >
    static std::enable_if_t< !traits::hasCopyAssignmentOp< U > >
    copyData( U &, U const & )
    {}
  };
  /// @endcond

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void copy( localIndex const sourceIndex, localIndex const destIndex ) override
  {
    if( sizedFromParent() )
    {
      copy_wrapper::copy( reference(), sourceIndex, destIndex );
    }
  }



  virtual void copyData( WrapperBase const & source ) override
  {
    Wrapper< T > const & castedSource = dynamicCast< Wrapper< T > const & >( source );
    copy_wrapper::copyData( *m_data, *castedSource.m_data );
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual void move( LvArray::MemorySpace const space, bool const touch ) const override
  { return wrapperHelpers::move( *m_data, space, touch ); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  virtual string typeRegex() const override
  { return TypeRegex< T >::get(); }

  ///@}

  /**
   * @name Methods for wrapped value data access
   */
  ///@{

  /**
   * @brief Accessor for m_data
   * @return reference to T
   */
  T & reference()
  { return *m_data; }

  /**
   * @brief const Accessor for m_data
   * @return reference to T, or in the case of an Array, a reference to an
   *         ArrayView<T const> const.
   */
  GEOSX_DECLTYPE_AUTO_RETURN reference() const
  { return referenceAsView(); }

  /**
   * @brief Provide access to wrapped object converted to a view, if possible.
   * @tparam U dummy template parameter to enable SFINAE (do not change)
   * @return the view for wrapped object
   *
   * This is used mainly for @p LvArray classes (arrays, etc.) that can convert
   * themselves into views. For other types, a regular reference is returned.
   */
  template< typename _T=T, typename=std::enable_if_t< traits::HasMemberFunction_toView< _T > > >
  GEOSX_DECLTYPE_AUTO_RETURN referenceAsView()
  { return m_data->toView(); }

  /**
   * @copydoc referenceAsView()
   */
  template< typename _T=T, typename=std::enable_if_t< !traits::HasMemberFunction_toView< _T > > >
  T & referenceAsView()
  { return *m_data; }

  /**
   * @copydoc referenceAsView()
   */
  template< typename _T=T, typename=std::enable_if_t< traits::HasMemberFunction_toView< _T > > >
  GEOSX_DECLTYPE_AUTO_RETURN referenceAsView() const
  { return m_data->toViewConst(); }

  /**
   * @copydoc referenceAsView()
   */
  template< typename _T=T, typename=std::enable_if_t< !traits::HasMemberFunction_toView< _T > > >
  T const & referenceAsView() const
  { return *m_data; }

  ///@}

  /**
   * @name Methods to manipulate default value for wrapped object
   */
  ///@{

  /**
   * @copydoc WrapperBase::hasDefaultValue()
   */
  virtual bool hasDefaultValue() const final override
  {
    return m_default.has_default_value;
  }

  /**
   * @brief Accessor for m_default.
   * @return reference to const m_default member
   */
  template< typename U=T >
  DefaultValue< T > const &
  getDefaultValueStruct() const
  {
    return m_default;
  }


  /**
   * @brief Accessor for default value.
   * @return reference to const T
   */
  template< typename U=T >
  std::enable_if_t< DefaultValue< U >::has_default_value, typename DefaultValue< U >::value_type const & >
  getDefaultValue() const
  {
    return m_default.value;
  }

  /**
   * @brief Setter for default value.
   * @param defaultVal the new default value
   * @return pointer to Wrapper<T>
   */
  template< typename U=T >
  std::enable_if_t< DefaultValue< U >::has_default_value, Wrapper< T > & >
  setDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    return *this;
  }

  /**
   * @brief Set and apply for default value.
   * @param defaultVal the new default value
   * @return pointer to Wrapper<T>
   */
  template< typename U=T >
  std::enable_if_t< !traits::is_array< U > && DefaultValue< U >::has_default_value, Wrapper< T > & >
  setApplyDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    *m_data = m_default.value;
    return *this;
  }

  /**
   * @brief Set and apply for default value.
   * @param defaultVal the new default value
   * @return pointer to Wrapper<T>
   */
  template< typename U=T >
  std::enable_if_t< traits::is_array< U > && DefaultValue< U >::has_default_value, Wrapper< T > & >
  setApplyDefaultValue( typename DefaultValue< U >::value_type const & defaultVal )
  {
    m_default.value = defaultVal;
    m_data->template setValues< serialPolicy >( m_default.value );
    return *this;
  }

  /**
   * @copydoc WrapperBase::getDefaultValueString()
   */
  virtual string getDefaultValueString() const override
  {
    // Find the dimensionality of the wrapper value
    string wrapper_type = rtTypes::typeNames( std::type_index( getTypeId()));
    integer value_dim = 0;
    if( wrapper_type.find( "array3d" ) != string::npos )
    {
      value_dim = 3;
    }
    else if( wrapper_type.find( "array2d" ) != string::npos )
    {
      value_dim = 2;
    }
    else if( wrapper_type.find( "array" ) != string::npos )
    {
      value_dim = 1;
    }

    // Compose the default string
    std::stringstream ss;

    for( integer ii=0; ii<value_dim; ++ii )
    {
      ss << "{";
    }

    ss << m_default;

    for( integer ii=0; ii<value_dim; ++ii )
    {
      ss << "}";
    }

    return ss.str();
  }

  virtual bool processInputFile( xmlWrapper::xmlNode const & targetNode ) override
  {
    InputFlags const inputFlag = getInputFlag();
    if( inputFlag >= InputFlags::OPTIONAL )
    {
      if( inputFlag == InputFlags::REQUIRED || !hasDefaultValue() )
      {
        bool const readSuccess = xmlWrapper::readAttributeAsType( reference(),
                                                                  getName(),
                                                                  targetNode,
                                                                  inputFlag == InputFlags::REQUIRED );
        GEOSX_THROW_IF( !readSuccess,
                        "Input variable " << getName() << " is required in " << targetNode.path() <<
                        ". Available options are: \n" << dumpInputOptions( true ) <<
                        "\nFor more details, please refer to documentation at: \n" <<
                        "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html \n",
                        InputError );
      }
      else
      {
        xmlWrapper::readAttributeAsType( reference(), getName(), targetNode, getDefaultValueStruct() );
      }

      return true;
    }

    return false;
  }

  /// @endcond DO_NOT_DOCUMENT

  /**
   * @brief Associate the path to this wrapper with the object being held.
   * @note This calls T::setName( string ) if it exists and is a no-op if not.
   *       LvArray objects implement this method.
   */
  void setName()
  { wrapperHelpers::setName( reference(), m_conduitNode.path() ); }


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void addBlueprintField( conduit::Node & fields,
                          string const & name,
                          string const & topology,
                          std::vector< string > const & componentNames = {} ) const override
  { wrapperHelpers::addBlueprintField( reference(), fields, name, topology, componentNames ); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void populateMCArray( conduit::Node & node, std::vector< string > const & componentNames = {} ) const override
  { wrapperHelpers::populateMCArray( reference(), node, componentNames ); }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  std::unique_ptr< WrapperBase > averageOverSecondDim( string const & name, Group & group ) const override
  {
    auto ptr = wrapperHelpers::averageOverSecondDim( reference() );
    using U = typename decltype( ptr )::element_type;

    GEOSX_ERROR_IF( ptr == nullptr, "Failed to average over the second dimension of." );

    return std::make_unique< Wrapper< U > >( name, group, std::move( ptr ) );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void registerToWrite() const override
  {
    m_conduitNode.reset();

    if( getRestartFlags() == RestartFlags::NO_WRITE )
    {
      return;
    }

    move( LvArray::MemorySpace::CPU, false );

    m_conduitNode[ "__sizedFromParent__" ].set( sizedFromParent() );

    wrapperHelpers::pushDataToConduitNode( *m_data, m_conduitNode );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void finishWriting() const override
  { m_conduitNode.reset(); }


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  bool loadFromConduit() override
  {
    if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
    {
      m_conduitNode.reset();
      return false;
    }

    setSizedFromParent( m_conduitNode[ "__sizedFromParent__" ].value() );

    wrapperHelpers::pullDataFromConduitNode( *m_data, m_conduitNode );

    m_conduitNode.reset();

    return true;
  }

  /**
   *  @name Overridden setters
   *  Setters that replace (hide) corresponding non-virtual functions in WrapperBase
   *  with the goal of returning Wrapper<T> pointers rather than WrapperBase pointers.
   */
  ///@{

  /*
   * @brief Set whether this wrapper is resized when its parent is resized.
   * @param val an int that is converted into a bool
   * @return a pointer to this wrapper
   */

  /**
   * @copydoc WrapperBase::setSizedFromParent(int)
   */
  Wrapper< T > & setSizedFromParent( int val )
  {
    WrapperBase::setSizedFromParent( val );
    return *this;
  }

  /**
   * @copydoc WrapperBase::setRestartFlags(RestartFlags)
   */
  Wrapper< T > & setRestartFlags( RestartFlags flags )
  {
    WrapperBase::setRestartFlags( flags );
    return *this;
  }

  /**
   * @copydoc WrapperBase::setPlotLevel(PlotLevel const)
   */
  Wrapper< T > & setPlotLevel( PlotLevel const flag )
  {
    WrapperBase::setPlotLevel( flag );
    return *this;
  }

  /**
   * @copydoc WrapperBase::setInputFlag(InputFlags const)
   */
  Wrapper< T > & setInputFlag( InputFlags const input )
  {
    WrapperBase::setInputFlag( input );
    return *this;
  }

  /**
   * @copydoc WrapperBase::setDescription(string const &)
   */
  Wrapper< T > & setDescription( string const & description )
  {
    WrapperBase::setDescription( description );
    return *this;
  }

  /**
   * @copydoc WrapperBase::setRegisteringObjects(string const &)
   */
  Wrapper< T > & setRegisteringObjects( string const & objectName )
  {
    WrapperBase::setRegisteringObjects( objectName );
    return *this;
  }

  ///@}

#if defined(USE_TOTALVIEW_OUTPUT)
  virtual string totalviewTypeName() const override
  {
    return LvArray::system::demangle( typeid( Wrapper< T > ).name() );
  }

  virtual int setTotalviewDisplay() const override
  {
    //std::cout<<"executing Wrapper::setTotalviewDisplay()"<<std::endl;
    WrapperBase::setTotalviewDisplay();
    TV_ttf_add_row( "m_ownsData", "bool", &m_ownsData );
    TV_ttf_add_row( "m_data", LvArray::system::demangle< T >().c_str(), m_data );
    TV_ttf_add_row( "m_default", LvArray::system::demangle< DefaultValue< T > >().c_str(), &m_default );
    return 0;
  }
//  void tvTemplateInstantiation();
#endif

#if defined(GEOSX_USE_PYGEOSX)
  virtual PyObject * createPythonObject( ) override
  { return wrapperHelpers::createPythonObject( reference() ); }
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
