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


#ifndef SRC_COMPONENTS_CORE_SRC_TRIBOLCOUPLING_HPP_
#define SRC_COMPONENTS_CORE_SRC_TRIBOLCOUPLING_HPP_

#include "dataRepository/ManagedGroup.hpp"

#define PARALLEL 1
#define HAVE_LLNL_GLOBALID 1
#include "SlideWorldAdapter.h"

namespace geosx
{

/**
 * @class TribolCoupling
 *
 * A class for managing Tribol code coupling.
 */
class TribolCoupling
{
public:
   static void Initialize(dataRepository::ManagedGroup * eventManager, dataRepository::ManagedGroup * domain) ;
   static void SyncTermination(int* terminate) ;
   static void SyncTimestep(real64* newDt) ;
   static void Cleanup() ;
   static vista::View *s_tribolProblem ;
   static vista::View *s_tribolDomain ;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_TRIBOLCOUPLING_HPP_ */
