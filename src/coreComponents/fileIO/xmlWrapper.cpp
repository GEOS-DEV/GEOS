/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 * @file xmlWrapper.cpp
 */

#include "ArrayUtilities.hpp"
#include "common/DataTypes.hpp"
#include "xmlWrapper.hpp"
#include "IntegerConversion.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "fileIO/utils/utils.hpp"

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

void xmlWrapper::StringToInputVariable( R1Tensor & target, string inputValue )
{
  string csvstr = inputValue;
  std::istringstream ss( csvstr );

  real64 value;
  int count = 0;
  while( ss.peek() == ',' || ss.peek() == ' ' )
  {
    ss.ignore();
  }
  while( !((ss>>value).fail()) )
  {
    target[count++] = value ;
    while( ss.peek() == ',' || ss.peek() == ' ' )
    {
      ss.ignore();
    }
  }
  GEOS_ERROR_IF(count!=3, "incorrect number of components specified for R1Tensor");
}

void xmlWrapper::addIncludedXML( xmlNode & targetNode )
{

  xmlNode rootNode = targetNode.root();
  string path = rootNode.child( filePathString ).attribute( filePathString ).value();

  xmlNode includedNode = targetNode.child("Included");

  for (xmlWrapper::xmlNode childNode=includedNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    // Get the child tag and name
    string childName = childNode.name();
    GEOS_ERROR_IF( childName!="File", "Child nodes of \"Included\" should be named \"File\" ");

    string filePathName = childNode.attribute("name").value();
    if( filePathName[0] != '/' )
    {
      filePathName = path + filePathName;
    }
    xmlDocument includedXmlDocument;
    xmlResult result;
    result = includedXmlDocument.load_file(filePathName.c_str());
    GEOS_ERROR_IF( !result, "Attempt to include file ("<<filePathName.c_str()<<") failed\n");

    for ( xmlNode importNode=includedXmlDocument.first_child() ; importNode ; importNode=importNode.next_sibling())
    {
      targetNode.append_copy(importNode);
    }
  }
}


} /* namespace geosx */
