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

/** @file */

#ifndef GEOSX_DATAREPOSITORY_VIEWWRAPPERBASE_HPP_
#define GEOSX_DATAREPOSITORY_VIEWWRAPPERBASE_HPP_

#include <string>
#include <memory>
#include "common/DataTypes.hpp"
#include "InputFlags.hpp"
#include "RestartFlags.hpp"

namespace axom
{
namespace sidre
{
class View;
}
}


namespace geosx
{
namespace dataRepository
{

class ManagedGroup;

/**
 * @class ViewWrapperBase
 */
class ViewWrapperBase
{
public:

  /**
   * @brief default destuctor
   */
  virtual ~ViewWrapperBase();

  /**
   * @brief constructor
   * @param[in] name name of the object
   * @param[in] parent pointer to ManagedGroup that holds this ViewWrapperBase
   */
  explicit ViewWrapperBase( string const & name,
                            ManagedGroup * const parent);


  /**
   * @brief move operator
   * @param[in] source
   */
  ViewWrapperBase( ViewWrapperBase&& source );


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
  virtual localIndex size(int i) const = 0;

  /**
   * @brief function to call T::resize( num_dims, dims )
   * @param[in] num_dims number of dimensions in T
   * @param[in] dims pointer to the new dims
   */
  virtual void resize(int num_dims, localIndex const * const dims) = 0;

  /**
   * @brief function to call T::resize( new_cap )
   * @param[in] new_cap the new capacity of the T
   */
  virtual void reserve(std::size_t new_cap) = 0;

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
  virtual void resize(localIndex newsize) = 0;


  /**
   * @brief     function to create a clone of *this ViewWrapperBase
   * @param[in] name name of the clone
   * @param[in] parent parent ManagedGroup that will hold this clone
   * @return
   *
   * The overridden function will create a copy of the derived ViewWrapper<T> the using the provided
   * values of name and parent to differentiate itself from the source.
   */
  virtual std::unique_ptr<ViewWrapperBase> clone( string const & name,
                                                  ManagedGroup * const parent ) = 0;

  /**
   *
   * @return
   */
  virtual bool shouldResize() const = 0;

  /**
   *
   * @return
   */
  virtual size_t sizeOfType() const = 0;

  /**
   *
   * @return
   */
  virtual bool shouldRegisterDataPtr() const = 0;

  /**
   *
   * @param view
   */
  virtual void registerDataPtr(axom::sidre::View * view=nullptr) const = 0; 

  /**
   *
   * @param view
   */
  virtual void registerToWrite(axom::sidre::View * view=nullptr) const = 0;

  /**
   *
   * @param view
   */
  virtual void finishWriting(axom::sidre::View * view=nullptr) const = 0;

  /**
   *
   * @param view
   */
  virtual void registerToRead(axom::sidre::View * view=nullptr) = 0;

  /**
   *
   * @param view
   */
  virtual void finishReading(axom::sidre::View * view=nullptr) = 0;

  /**
   * @brief function to call resize( newsize ) where newsize is taken from the parent ManagedGroup
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
  virtual localIndex Pack( char *& buffer ) const = 0;

  /**
   *
   * @param buffer
   * @param packList
   * @return
   */
  virtual localIndex Pack( char *& buffer, arrayView1d<localIndex const> const & packList ) const = 0;

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
  virtual localIndex PackSize( arrayView1d<localIndex const> const & packList ) const = 0;

  /**
   *
   * @param buffer
   * @return
   */
  virtual localIndex Unpack( char const *& buffer ) = 0;

  /**
   *
   * @param buffer
   * @param unpackIndices
   * @return
   */
  virtual localIndex Unpack( char const *& buffer, arrayView1d<localIndex const> const & unpackIndices ) = 0;

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
  ViewWrapperBase * setSizedFromParent( int val )
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
  ViewWrapperBase * setRestartFlags( RestartFlags flags)
  {
    m_restart_flags = flags;
    return this;
  }

#ifdef GEOSX_USE_ATK
  /**
   *
   * @return
   */
  axom::sidre::View * getSidreView() const
  {
    return m_sidreView;
  }
#endif

  /**
   *
   * @return
   */
  PlotLevel getPlotLevel() const {return m_plotLevel;}

  /**
   *
   * @param flag
   * @return
   */
  ViewWrapperBase * setPlotLevel( PlotLevel const flag )
  {
    m_plotLevel = flag;
    return this;
  }

  /**
   *
   * @param flag
   * @return
   */
  ViewWrapperBase * setPlotLevel( int const flag )
  {
    m_plotLevel = IntToPlotLevel(flag);
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
  ViewWrapperBase * setInputFlag( InputFlags const input )
  {
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
  ViewWrapperBase * setDescription( string const & description )
  {
    m_description = description;
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


private:

  /// name of the object that is being wrapped
  string m_name;

  /// pointer to ManagedGroup that holds this ViewWrapperBase
  ManagedGroup * m_parent;

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

  #ifdef GEOSX_USE_ATK
  /// a pointer to the corrosponding sidre view
  axom::sidre::View* m_sidreView;
#endif


  ViewWrapperBase() = delete;
  ViewWrapperBase( ViewWrapperBase const & ) = delete;
  ViewWrapperBase& operator=( ViewWrapperBase const & ) = delete;
  ViewWrapperBase& operator=( ViewWrapperBase&& ) = delete;

};

}
} /* namespace geosx */

#endif
