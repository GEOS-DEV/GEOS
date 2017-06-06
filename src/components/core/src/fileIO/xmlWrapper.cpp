/*
 * xmlWrapper.cpp
 *
 *  Created on: Jun 3, 2017
 *      Author: rrsettgast
 */

#include "common/DataTypes.hpp"
#include "xmlWrapper.hpp"
#include "DocumentationNode.hpp"
#include <slic/slic.hpp>
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/StringUtilities.hpp"

using namespace cxx_utilities;

namespace geosx
{
using namespace dataRepository;

xmlWrapper::xmlWrapper()
{
  // TODO Auto-generated constructor stub

}

xmlWrapper::~xmlWrapper()
{
  // TODO Auto-generated destructor stub
}



void xmlWrapper::ReadAttributeAsType( dataRepository::ManagedGroup & group,
                                      DocumentationNode const & subDocNode,
                                      xmlNode const & targetNode )
{
  std::string childType = subDocNode.getSchemaType();
  rtTypes::TypeIDs const typeID = rtTypes::typeID(childType);
  rtTypes::ApplyIntrinsicTypeLambda2 ( typeID,
                                       [&]( auto a, auto b ) -> void
  {
    string defVal = subDocNode.getDefault();

    pugi::xml_attribute xmlatt = targetNode.attribute(subDocNode.getStringKey().c_str());
    ViewWrapper<decltype(a)>& dataView = group.getWrapper<decltype(a)>(subDocNode.getStringKey());
    std::vector<decltype(b)> xmlVal;

    if( !xmlatt.empty() )
    {
      as_type( xmlVal, xmlatt.value(), defVal );
    }
    else
    {
      if( defVal == "REQUIRED")
      {
        string message = "variable " + subDocNode.getName() + " is required in " + targetNode.path();
        SLIC_ERROR( message );
      }
      else
      {
        stringutilities::StringToType( xmlVal, defVal );
      }


    }
    localIndex const size = xmlVal.size();
    dataView.resize( size );
    typename ViewWrapper<decltype(a)>::rtype data = dataView.data();
//        decltype(a) * data = dataView.pointer();
    cxx_utilities::equateStlVector(data,xmlVal);
  });


}

} /* namespace geosx */
