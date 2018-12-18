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

/*
 * DataObjectBase.hpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

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

class ViewWrapperBase
{
public:

  /*!
   * \brief default destuctor
   */
  virtual ~ViewWrapperBase();

  /*!
   *
   * @param name name of the object
   * \brief constructor
   */
  explicit ViewWrapperBase( std::string const & name,
                            ManagedGroup * const parent);


  ViewWrapperBase( ViewWrapperBase&& source );


  /*!
   *
   * @return type_info of the DataObject
   */
  virtual std::type_info const & get_typeid() const = 0;


  virtual bool empty() const = 0;
  virtual localIndex size() const = 0;
  virtual int numDimensions() const = 0;
  virtual localIndex size(int i) const = 0;
  virtual void resize(int num_dims, long long const * const dims) = 0;
  virtual void reserve(std::size_t new_cap) = 0;
  virtual std::size_t capacity() const = 0;
  virtual std::size_t max_size() const = 0;
  virtual void clear() = 0;
  virtual void insert() = 0;
  virtual void resize(localIndex newsize) = 0;


  /**
   * @brief  function to create a clone of *this ViewWrapperBase
   * @param name name of the clone
   * @param parent parent ManagedGroup that will hold this clone
   * @return
   *
   * The overridden function will create a copy of the derived ViewWrapper<T> the using the provided
   * values of name and parent to differentiate itself from the source.
   */
  virtual std::unique_ptr<ViewWrapperBase> clone( string const & name,
                                                  ManagedGroup * const parent ) = 0;

  virtual bool shouldResize() const = 0;

  virtual size_t sizeOfType() const = 0;

  virtual bool shouldRegisterDataPtr() const = 0;
  virtual void registerDataPtr(axom::sidre::View * view=nullptr) const = 0; 
  virtual void registerToWrite(axom::sidre::View * view=nullptr) const = 0;
  virtual void finishWriting(axom::sidre::View * view=nullptr) const = 0;
  virtual void registerToRead(axom::sidre::View * view=nullptr) = 0;
  virtual void finishReading(axom::sidre::View * view=nullptr) = 0;


  void resize();

  virtual void copy( localIndex const sourceIndex, localIndex const destIndex ) = 0;

  virtual bool isPackable() const = 0;
  virtual localIndex Pack( char *& buffer ) const = 0;
  virtual localIndex Pack( char *& buffer, arrayView1d<localIndex const> const & packList ) const = 0;
  virtual localIndex PackSize( ) const = 0;
  virtual localIndex PackSize( arrayView1d<localIndex const> const & packList ) const = 0;

  virtual localIndex Unpack( char const *& buffer ) = 0;
  virtual localIndex Unpack( char const *& buffer, arrayView1d<localIndex const> const & unpackIndices ) = 0;

//  virtual int PackingSize( char *& buffer, localIndex_array const & packList ) = 0;

  int sizedFromParent() const
  {
    return m_sizedFromParent;
  }

  void setSizedFromParent( int val )
  {
    m_sizedFromParent = val;
  }

  RestartFlags getRestartFlags() const { return m_restart_flags; }

  void setRestartFlags( RestartFlags flags) { m_restart_flags = flags; } 

#ifdef GEOSX_USE_ATK
  axom::sidre::View * getSidreView() const
  {
    return m_sidreView;
  }
#endif

  PlotLevel getPlotLevel() const {return m_plotLevel;}

  void setPlotLevel( PlotLevel const flag )
  {
    m_plotLevel = flag;
  }

  ViewWrapperBase * setPlotLevel( int const flag )
  {
    m_plotLevel = IntToPlotLevel(flag);
    return this;
  }

  string const & getName() const
  {
    return m_name;
  }

  ViewWrapperBase * setInputFlag( InputFlags const input )
  {
    m_inputFlag = input;
    return this;
  }

  InputFlags getInputFlag() const
  {
    return m_inputFlag;
  }

  ViewWrapperBase * setDescription( string const & description )
  {
    m_description = description;
    return this;
  }

  string const & getDescription() const
  {
    return m_description;
  }


private:
  std::string m_name;
  ManagedGroup* m_parent;
  int m_sizedFromParent;
  RestartFlags m_restart_flags;
  PlotLevel m_plotLevel;
  InputFlags m_inputFlag;
  string m_description;
#ifdef GEOSX_USE_ATK
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
