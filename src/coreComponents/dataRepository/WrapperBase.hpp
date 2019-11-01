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

/** @file */

#ifndef GEOSX_DATAREPOSITORY_WRAPPERBASE_HPP_
#define GEOSX_DATAREPOSITORY_WRAPPERBASE_HPP_

#include <string>
#include <memory>
#include "common/DataTypes.hpp"
#include "InputFlags.hpp"
#include "RestartFlags.hpp"

namespace conduit
{
class Node;
}


namespace geosx
{
namespace dataRepository
{

class Group;

/**
 * @class WrapperBase
 */
class WrapperBase
{
public:

  /**
   * @brief constructor
   * @param[in] name name of the object
   * @param[in] parent pointer to Group that holds this WrapperBase
   */
  explicit WrapperBase( string const & name,
                        Group * const parent );

  WrapperBase() = delete;
  WrapperBase( WrapperBase const & ) = delete;
  WrapperBase( WrapperBase && ) = delete;
  WrapperBase & operator=( WrapperBase const & ) = delete;
  WrapperBase & operator=( WrapperBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~WrapperBase();

  /**
   * @brief move operator
   * @param[in] source
   */
  // WrapperBase( WrapperBase && source );


  virtual void CopyWrapperAttributes( WrapperBase const & source );

  /**
   * @brief Virtual function to return the typeid of T.
   * @return type_info of the wrapped type "typeid(T)"
   */
  virtual std::type_info const & get_typeid() const = 0;


  /**
   * @brief function call T::empty()
   * @return boolean true if T is empty, false if not.
   */
  virtual bool empty() const = 0;

  /**
   * @brief function to call T::size()
   * @return result of T::size()
   */
  virtual localIndex size() const = 0;

  /**
   * @brief function to call T::numDimensions()
   * @return result of T::numDimensions()
   */
  virtual int numDimensions() const = 0;

  /**
   * @brief function to call T::size(int)
   * @return result of T::size(int)
   */
  virtual localIndex size( int i ) const = 0;

  /**
   * @brief function to call T::resize( num_dims, dims )
   * @param[in] num_dims number of dimensions in T
   * @param[in] dims pointer to the new dims
   */
  virtual void resize( int num_dims, localIndex const * const dims ) = 0;

  /**
   * @brief function to call T::resize( new_cap )
   * @param[in] new_cap the new capacity of the T
   */
  virtual void reserve( std::size_t new_cap ) = 0;

  /**
   * @brief function to call T::capacity()
   * @return result of T::capacity()
   */
  virtual std::size_t capacity() const = 0;

  /**
   * @brief function to call T::max_size()
   * @return result of T::max_size()
   */
  virtual std::size_t max_size() const = 0;

  /**
   * @brief function to call T::clear()
   * @return result of T::clear()
   */
  virtual void clear() = 0;

  /**
   * @brief function to call T::insert()
   * @return result of T::insert()
   */
  virtual void insert() = 0;

  /**
   * @brief function to call T::resize(newsize)
   * @param[in] newsize parameter to pass to T::resize(newsize)
   * @return result of T::resize(newsize)
   */
  virtual void resize( localIndex newsize ) = 0;


  /**
   * @brief     function to create a clone of *this WrapperBase
   * @param[in] name name of the clone
   * @param[in] parent parent Group that will hold this clone
   * @return
   *
   * The overridden function will create a copy of the derived Wrapper<T> the using the provided
   * values of name and parent to differentiate itself from the source.
   */
  virtual std::unique_ptr< WrapperBase > clone( string const & name,
                                                Group * const parent ) = 0;

  virtual void move( chai::ExecutionSpace space, bool touch ) = 0;

  /**
   *
   * @return
   */
  virtual size_t sizeOfType() const = 0;

  /**
   *
   * @param view
   */
  virtual void registerToWrite() = 0;

  /**
   *
   * @param view
   */
  virtual void finishWriting() = 0;

  /**
   *
   * @param view
   */
  virtual void loadFromConduit() = 0;

  /**
   * @brief function to call resize( newsize ) where newsize is taken from the parent Group
   */
  void resize();

  /**
   *
   * @param sourceIndex
   * @param destIndex
   */
  virtual void copy( localIndex const sourceIndex, localIndex const destIndex ) = 0;

  /**
   *
   * @return
   */
  virtual bool isPackable() const = 0;

  /**
   *
   * @param buffer
   * @return
   */
  virtual localIndex Pack( buffer_unit_type * & buffer ) const = 0;

  /**
   *
   * @param buffer
   * @param packList
   * @return
   */
  virtual localIndex Pack( buffer_unit_type * & buffer, arrayView1d< localIndex const > const & packList ) const = 0;

  /**
   *
   * @return
   */
  virtual localIndex PackSize( ) const = 0;

  /**
   *
   * @param packList
   * @return
   */
  virtual localIndex PackSize( arrayView1d< localIndex const > const & packList ) const = 0;

  /**
   *
   * @param buffer
   * @return
   */
  virtual localIndex Unpack( buffer_unit_type const * & buffer ) = 0;

  /**
   *
   * @param buffer
   * @param unpackIndices
   * @return
   */
  virtual localIndex Unpack( buffer_unit_type const * & buffer, arrayView1d< localIndex const > const & unpackIndices ) = 0;

  /**
   *
   * @return
   */
  int sizedFromParent() const
  {
    return m_sizedFromParent;
  }

  /**
   *
   * @param val
   * @return
   */
  WrapperBase * setSizedFromParent( int val )
  {
    m_sizedFromParent = val;
    return this;
  }

  /**
   *
   * @return
   */
  RestartFlags getRestartFlags() const { return m_restart_flags; }

  /**
   *
   * @param flags
   * @return
   */
  WrapperBase * setRestartFlags( RestartFlags flags )
  {
    m_restart_flags = flags;
    return this;
  }


  /**
   *
   * @return
   */
  PlotLevel getPlotLevel() const { return m_plotLevel; }

  /**
   *
   * @param flag
   * @return
   */
  WrapperBase * setPlotLevel( PlotLevel const flag )
  {
    m_plotLevel = flag;
    return this;
  }

  /**
   *
   * @param flag
   * @return
   */
  WrapperBase * setPlotLevel( int const flag )
  {
    m_plotLevel = IntToPlotLevel( flag );
    return this;
  }

  /**
   *
   * @return
   */
  string const & getName() const
  {
    return m_name;
  }

  /**
   *
   * @param input
   * @return
   */
  WrapperBase * setInputFlag( InputFlags const input )
  {
    if( input == InputFlags::OPTIONAL || input == InputFlags::REQUIRED )
    {
      this->setSizedFromParent( 0 );
    }
    m_inputFlag = input;
    return this;
  }

  /**
   *
   * @return
   */
  InputFlags getInputFlag() const
  {
    return m_inputFlag;
  }

  /**
   *
   * @param description
   * @return
   */
  WrapperBase * setDescription( string const & description )
  {
    m_description = description;
    return this;
  }

  std::vector< string > const & getRegisteringObjects() const
  {
    return m_registeringObjects;
  }

  WrapperBase * setRegisteringObjects( string const & objectName )
  {
    m_registeringObjects.push_back( objectName );
    return this;
  }

  /**
   *
   * @return
   */
  string const & getDescription() const
  {
    return m_description;
  }

#if defined(USE_TOTALVIEW_OUTPUT)
  /**
   * @brief Virtual function to return the the typename for a Wrapper derived type that is
   *                represented by a WrapperBase *.
   * @return A string that contains the typename for use in totalview.
   */
  virtual string totalviewTypeName() const = 0;

  /**
   * @brief Function to execute the TV_tff_add_row() calls that represent each data member that
   *        will be displayed.
   * @return 0
   */
  virtual int setTotalviewDisplay() const;
//  static int TV_ttf_display_type( const WrapperBase * wrapper);
#endif

protected:

  conduit::Node & getConduitNode()
  {
    return m_conduitNode;
  }

private:

  /// name of the object that is being wrapped
  string m_name;

  /// pointer to Group that holds this WrapperBase
  Group * m_parent;

  /// integer to indicate whether or not this wrapped object should be resized when m_parent is resized
  int m_sizedFromParent;

  /// flag to determine the restart behavior for this wrapped object
  RestartFlags m_restart_flags;

  /// flag to store the plotLevel
  PlotLevel m_plotLevel;

  /// flag to store if this wrapped object should be read from input
  InputFlags m_inputFlag;

  /// a string description of the wrapped object
  string m_description;

  std::vector< string > m_registeringObjects;

  /// a reference to the corresponding conduit::Node
  conduit::Node & m_conduitNode;
};

}
} /* namespace geosx */

#endif /* GEOSX_DATAREPOSITORY_WRAPPERBASE_HPP_ */
