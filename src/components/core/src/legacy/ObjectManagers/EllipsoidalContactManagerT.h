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

/**
 * @file ContactManagerT.h
 * @author Scott Johnson
 * @date July 13, 2011
 */
#ifndef ELLIPSOIDALCONTACTMANAGERT_H_
#define ELLIPSOIDALCONTACTMANAGERT_H_

#include "../../dataRepository/Group.hpp"
#include "Common/Common.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ContactManagerT.h"

/**
 * @author Scott Johnson
 * @brief Class to manage collections of ellipsoid-ellipsoid contacts
 */
class EllipsoidalContactManagerT : public ContactManagerBaseT
{
public:
  /**
   * @brief Ellipsoidal contact manager constructor
   * @author Scott Johnson
   */
  EllipsoidalContactManagerT();

  /**
   * @brief Contact manager destructor
   * @author Scott Johnson
   */
  virtual ~EllipsoidalContactManagerT(){}
};

#endif /* ELLIPSOIDALCONTACTMANAGERT_H_ */
