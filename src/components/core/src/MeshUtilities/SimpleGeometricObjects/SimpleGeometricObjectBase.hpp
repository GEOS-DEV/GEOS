/*
 * SimpleGeometricObjects.h
 *
 *  Created on: Dec 4, 2012
 *      Author: settgast1
 */

#ifndef SIMPLEGEOMETRICOBJECTS_H_
#define SIMPLEGEOMETRICOBJECTS_H_

//#include "common/Common.h"
#include "codingUtilities/StringUtilities.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "ObjectCatalog.hpp"

class Function;

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const geometricObjects("GeometricObjects");
}
}


class SimpleGeometricObjectBase// : public dataRepository::ManagedGroup
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  using CatalogInterface = cxx_utilities::CatalogInterface< SimpleGeometricObjectBase >;

  static typename CatalogInterface::CatalogType& GetCatalog()
  {
    static SimpleGeometricObjectBase::CatalogInterface::CatalogType catalog;
    return catalog;
  }
  ///@}

  SimpleGeometricObjectBase(  )
  {}

  virtual ~SimpleGeometricObjectBase()
  {}

//  virtual void ReadXML( pugi::xml_node& hdn ) = 0;

  virtual bool IsCoordInObject( const R1Tensor& coord ) const = 0;

  virtual void ReadXML( xmlWrapper::xmlNode const & xmlNode ) = 0;
};


}
#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
