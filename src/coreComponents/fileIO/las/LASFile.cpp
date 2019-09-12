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

#include "LASFile.hpp"

namespace geosx
{

LASInformationSection * LASInformationSection::CreateLASInformationSection( char const & name )
{
  if( name == 'V' )          // Version Information
  {
    return new LASVersionInformationSection();
  }
  else if( name == 'W' )     // Well Information
  {
    return new LASWellInformationSection();
  }
  else if( name == 'C' )     // Curve Information
  {
    return new LASCurveInformationSection();
  }
  else if( name == 'P' )     // Parameter Information
  {
    return new LASParameterInformationSection();
  }
  else if( name == 'O' )     // Other Information
  {
    return new LASOtherInformationSection();
  }
  else
  {
    GEOS_ERROR( name << " is not a valid section for LAS files" );
    return nullptr;
  }
}
}
