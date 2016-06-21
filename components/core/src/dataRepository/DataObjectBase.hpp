/*
 * DataObjectBase.hpp
 *
 *  Created on: Jun 17, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAOBJECTBASE_HPP_
#define COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAOBJECTBASE_HPP_

#include <string>

namespace geosx {


template< typename T > class DataObject;

class DataObjectBase
{
public:
  /*!
   *
   * @param name name of the object
   * \brief constructor
   */
  DataObjectBase( const std::string& name );

  /*!
   * \brief default destuctor
   */
  virtual ~DataObjectBase();

  /*!
   *
   * @return type_info of the DataObject
   */
  virtual const std::type_info& get_typeid() const = 0 ;

  virtual void resize( const std::size_t newsize ) = 0;

  virtual std::size_t size( ) = 0;


  template< typename T >
  static std::unique_ptr<DataObjectBase> Factory( const std::string& name )
  {
    return std::move(std::unique_ptr< DataObject<T> >( new DataObject<T>( name ) ) );
//    return std::move(std::make_unique<DataObject<T> >( name ) );
  }

  template< typename T >
  DataObject<T>& getObject()
  {
    return dynamic_cast<DataObject<T>&>(*this);
  }



//#ifndef OBJECTDATA_PTR_RETURN

  template< typename T >
  typename DataObject<T>::const_rtype getObjectData() const
  { return (dynamic_cast<const DataObject<T>&>(*this)).getObjectData(); }

  template< typename T >
  typename DataObject<T>::rtype getObjectData()
  { return const_cast<typename DataObject<T>::rtype>( const_cast<const DataObjectBase*>(this)->getObjectData<T>() ); }
/*
#else

  template< typename T >
  typename std::enable_if<!(std::is_fundamental<T>::value), T >::type::type* getObjectData()
  {
    return (dynamic_cast<DataObject<T>&>(*this).getObjectData()).data();
  }

  template< typename T >
  const typename std::enable_if<std::is_fundamental<T>::value, T >::type& getObjectData() const
  {
    return &((dynamic_cast<DataObject<T>&>(*this)).getObjectData());
  }


#endif
*/

private:
  std::string m_fieldName;
//  attributeMap m_attributes;

  DataObjectBase();
  DataObjectBase(const DataObjectBase&);

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_DATAOBJECTBASE_HPP_ */
