// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include "bufvector.h"
#include "Utilities/StringUtilities.h"

const unsigned bufvector::sizeOfPackedFieldString = 20;

unsigned int bufvector::PackFieldname( std::string fieldname ){ // nb not passed
                                                                // by reference
  fieldname.resize(sizeOfPackedFieldString,' ');
  return Pack(fieldname);
}


unsigned int bufvector::PackFieldname(char*& buffer,  std::string fieldname ){ // nb
                                                                               // not
                                                                               // passed
                                                                               // by
                                                                               // reference
  fieldname.resize(sizeOfPackedFieldString,' ');
  return Pack(buffer, fieldname);
}

bool bufvector::FieldnameMatchesIdString(std::string fieldname, std::string idString){// nb
                                                                                      // not
                                                                                      // passed
                                                                                      // by
                                                                                      // reference
  fieldname.resize(sizeOfPackedFieldString,' ');
  //idString.resize(sizeOfPackedFieldString,' ');
  return streq(fieldname,idString);
}
