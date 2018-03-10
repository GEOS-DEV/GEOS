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

#ifdef USE_ATK
namespace axom
{
namespace sidre
{
class View;
}
}
#endif

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
  virtual bool shouldResize() const = 0;

  virtual void registerDataPtr(axom::sidre::View * view=nullptr) const = 0; 
  virtual void registerToWrite(axom::sidre::View * view=nullptr) const = 0;
  virtual void finishWriting(axom::sidre::View * view=nullptr) const = 0;
  virtual void registerToRead(axom::sidre::View * view=nullptr) = 0;
  virtual void finishReading(axom::sidre::View * view=nullptr) = 0;


  void resize();

  int sizedFromParent() const
  {
    return m_sizedFromParent;
  }

  void setSizedFromParent( int val )
  {
    m_sizedFromParent = val;
  }

  bool getWriteToRestart() const
  { 
    return m_write_to_restart;
  }

  void setWriteToRestart( bool write_to_restart )
  {
    m_write_to_restart = write_to_restart;
  }

  bool getReadFromRestart() const
  { 
    return m_read_from_restart;
  }

  void setReadFromRestart( bool read_from_restart )
  {
    m_read_from_restart = read_from_restart;
  }

#ifdef USE_ATK
  axom::sidre::View * getSidreView() const
  {
    return m_sidreView;
  }
#endif

  string const & getName() const
  {
    return m_name;
  }


private:
  std::string m_name;
  ManagedGroup* m_parent;
  int m_sizedFromParent;
  bool m_write_to_restart;
  bool m_read_from_restart;

#ifdef USE_ATK
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
