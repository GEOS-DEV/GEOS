// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file EllipsoidalContactManagerT.cpp
 * @author Scott Johnson
 * @date created on July 13, 2011
 */

#include "ContactManagerBaseT.h"
#include "EllipsoidalContactManagerT.h"

/**
 * @brief Constructor for the contacts / neighbor pairs
 * @author Scott Johnson
 * @date Jun 15, 2011
 *
 * note: any of the models can use latching spring
 * note: any of the cohesive models can be told to be cementitious or adhesive
 *
 */
EllipsoidalContactManagerT::EllipsoidalContactManagerT():
  ContactManagerBaseT(ObjectDataStructureBaseT::EllipsoidalContactManager)
{}
