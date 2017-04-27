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
class SimpleGeometricObjectBase : public dataRepository::ManagedGroup
{
public:

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  using CatalogInterface = cxx_utilities::CatalogInterface< SimpleGeometricObjectBase, string const & >;

  static typename CatalogInterface::CatalogType& GetCatalog()
  {
    static SimpleGeometricObjectBase::CatalogInterface::CatalogType catalog;
    return catalog;
  }
  ///@}

  SimpleGeometricObjectBase( string const & name, ManagedGroup * const parent ):
    ManagedGroup(name,parent)
  {

  }

  virtual ~SimpleGeometricObjectBase()
  {
  }

//  virtual void ReadXML( pugi::xml_node& hdn ) = 0;

  virtual bool IsCoordInObject( const R1Tensor& coord ) = 0;

};


}
#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
