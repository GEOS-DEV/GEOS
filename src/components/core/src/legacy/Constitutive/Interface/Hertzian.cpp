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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */



/*
 * Hertzian.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */


#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "Hertzian.h"
#include "Constitutive/Interface/InterfaceFactory.h"
#include <typeinfo>
#include <assert.h>
#include "../../IO/ticpp/HierarchicalDataNode.h.old"
#include "../../IO/ticpp/HierarchicalDataNode.h.old"

Hertzian::Hertzian( ):
  HertzianIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

Hertzian::~Hertzian()
{}


/// Register class in the class factory
REGISTER_INTERFACE( Hertzian )
