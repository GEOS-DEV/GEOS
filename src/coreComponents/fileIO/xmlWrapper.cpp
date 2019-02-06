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

/*
 * xmlWrapper.cpp
 *
 *  Created on: Jun 3, 2017
 *      Author: rrsettgast
 */

#include "ArrayUtilities.hpp"
#include "common/DataTypes.hpp"
#include "xmlWrapper.hpp"
#include "IntegerConversion.hpp"
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

} /* namespace geosx */
