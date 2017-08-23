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
                            ManagedGroup * const parent );


  ViewWrapperBase( ViewWrapperBase&& source );


  /*!
   *
   * @return type_info of the DataObject
   */
  virtual std::type_info const & get_typeid() const = 0;


  virtual bool empty() const = 0;
  virtual localIndex size( ) const = 0;
  virtual void reserve( std::size_t new_cap ) = 0;
  virtual std::size_t capacity() const = 0;
  virtual std::size_t max_size() const = 0;
  virtual void clear() = 0;
  virtual void insert() = 0;
  virtual void registerDataPtr() = 0;
  virtual void unregisterDataPtr() = 0;
  virtual void resizeFromSidre() = 0;
  virtual void storeSizedFromParent() = 0;
  virtual void loadSizedFromParent() = 0;

  virtual void resize( localIndex newsize ) = 0;

  int sizedFromParent() const
  {
    return m_sizedFromParent;
  }

  void setSizedFromParent( int val )
  {
    m_sizedFromParent = val;
  }

  axom::sidre::View const * getSidreView() const
  {
    return m_sidreView;
  }
  axom::sidre::View * getSidreView()
  {
    return m_sidreView;
  }

  string const & getName() const
  {
    return m_name;
  }

private:
  std::string m_name;
  ManagedGroup* m_parent;
  int m_sizedFromParent;
  axom::sidre::View* m_sidreView;

  ViewWrapperBase() = delete;
  ViewWrapperBase( ViewWrapperBase const & ) = delete;
  ViewWrapperBase& operator=( ViewWrapperBase const & ) = delete;
  ViewWrapperBase& operator=( ViewWrapperBase&& ) = delete;

};

}
} /* namespace geosx */

#endif
