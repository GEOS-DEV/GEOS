/*
 * DataObjectBase.hpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAOBJECTBASE_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAOBJECTBASE_HPP_

#include <string>
#include <memory>

namespace asctoolkit
{
namespace sidre
{
class DataView;
}
}

namespace geosx {


template< typename T > class DataObject;

class DataObjectBase
{
public:

  /*!
   * \brief default destuctor
   */
  virtual ~DataObjectBase();

  /*!
   *
   * @param name name of the object
   * \brief constructor
   */
  DataObjectBase( std::string const & name );


  DataObjectBase( DataObjectBase&& source );


  /*!
   *
   * @return type_info of the DataObject
   */
  virtual std::type_info const & get_typeid() const = 0 ;


//  virtual bool empty() const = 0;


  virtual std::size_t size( ) = 0;
  /*
  size_type max_size() const;


  virtual void clear() = 0 ;
  virtual void insert() = 0;
  iterator erase( iterator pos );
  iterator erase( const_iterator pos );
  iterator erase( const_iterator first, const_iterator last );
  size_type erase( const key_type& key );

  iterator erase( const_iterator pos );
  iterator erase( iterator first, iterator last );
  iterator erase( const_iterator first, const_iterator last );

  void swap( unordered_map& other );
  void swap( vector& other );
*/

  virtual void resize( std::size_t newsize ) = 0;



  template< typename T >
  static std::unique_ptr<DataObjectBase> Factory( std::string const & name )
  {
    return std::move(std::make_unique<DataObject<T> >( name ) );
  }

  template< typename T >
  DataObject<T>& getObject()
  {
    return dynamic_cast<DataObject<T>&>(*this);
  }

  template< typename T >
  DataObject<T> const & getObject() const
  {
    return dynamic_cast<DataObject<T> const &>(*this);
  }



  template< typename T >
  typename DataObject<T>::const_rtype getObjectData() const
  { return (dynamic_cast<DataObject<T> const &>(*this)).getObjectData(); }

  template< typename T >
  typename DataObject<T>::rtype getObjectData()
  { return (dynamic_cast<DataObject<T> &>(*this)).getObjectData(); }
//  { return const_cast<typename DataObject<T>::rtype>( const_cast<DataObjectBase const *>(this)->getObjectData<T>() ); }


private:
  std::string m_name;

  asctoolkit::sidre::DataView* m_sidreView;

  DataObjectBase() = delete;
  DataObjectBase( DataObjectBase const & ) = delete;
  DataObjectBase& operator=( DataObjectBase const & ) = delete;
  DataObjectBase& operator=( DataObjectBase&& ) = delete;

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAOBJECTBASE_HPP_ */
