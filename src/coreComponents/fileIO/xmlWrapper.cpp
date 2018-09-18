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
 * xmlWrapper.cpp
 *
 *  Created on: Jun 3, 2017
 *      Author: rrsettgast
 */

#include "common/DataTypes.hpp"
#include "xmlWrapper.hpp"
#include "DocumentationNode.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/StringUtilities.hpp"

#ifdef GEOSX_USE_ATK
#include <slic/slic.hpp>
#endif


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
      ViewWrapper<decltype(a)>& dataView = *(group.getWrapper<decltype(a)>(subDocNode.getStringKey()));
      std::vector<decltype(b)> xmlVal;

      if( !xmlatt.empty() )
      {
        as_type( xmlVal, xmlatt.value(), defVal );
      }
      else
      {
        if( defVal == "REQUIRED" )
        {
          string message = "variable " + subDocNode.getName() + " is required in " + targetNode.path();
#ifdef GEOSX_USE_ATK
          SLIC_ERROR( message );
#endif
        }
        else
        {
          stringutilities::StringToType( xmlVal, defVal );
        }


      }
      localIndex const size = multidimensionalArray::integer_conversion<localIndex>(xmlVal.size());
      dataView.resize( size );
      typename ViewWrapper<decltype(a)>::rtype data = dataView.data();
//        decltype(a) * data = dataView.pointer();
      cxx_utilities::equateStlVector(data,xmlVal);
    });


}



R1Tensor xmlWrapper::as_type( xmlNode const & node, std::string const name, R1Tensor defValue )
{
  R1Tensor rval = defValue;
  pugi::xml_attribute att = node.attribute( name.c_str() );

  if( !att.empty() )
  {
    string inputValue = att.value();
    if( inputValue!="" )
    {
      std::string csvstr = inputValue;
      std::istringstream ss( csvstr );

      while( ss.peek() == ',' || ss.peek() == ' ' )
      {
        ss.ignore();
      }
      for( int i=0 ; i<3 ; ++i )
      {
        ss>>rval[i];
        while( ss.peek() == ',' || ss.peek() == ' ' )
        {
          ss.ignore();
        }
      }
    }
  }

  return rval;
}

} /* namespace geosx */
