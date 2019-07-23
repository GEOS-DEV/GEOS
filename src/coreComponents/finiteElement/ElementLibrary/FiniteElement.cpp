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

#include "FiniteElement.h"


namespace geosx
{



/**
 * Explicit instantiations
 */

//template class FiniteElement<1>;
//template class FiniteElement<2>;
//template class FiniteElement<3>;

typedef FiniteElement<3> FiniteElement3d;

REGISTER_CATALOG_ENTRY( FiniteElementBase, FiniteElement3d, BasisBase const &, QuadratureBase const &, const int )


}
