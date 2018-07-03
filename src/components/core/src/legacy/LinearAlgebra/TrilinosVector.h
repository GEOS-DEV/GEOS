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

#ifndef TRILINOS_VECTOR_H
#define TRILINOS_VECTOR_H

/**
 * @file TrilinosVector.h
 * @author white230
 */

#include "Common/Common.h"

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "Epetra_FEVector.h"
#include "Teuchos_RCP.hpp"

/** \namespace Trilinos */

namespace Trilinos
{


/**
 * Vector class which wraps an underlying Trilinos data type.  The
 * vector may be partitioned across parallel processes.
 */

class Vector
{
public:

  /** Default (empty) constructor. */

  Vector();

  /** Copy constructor */

  //Vector(const Vector &vector);

  /** Destructor */

  //~Vector();

  /** Constructor taking a parallel partitioning map */

  //Vector(const &Partitioning);

  /** Reinit, taking a parallel partitioning map */

  //reinit(const &Partitioning);

  /** Compute L2 norm of the vector */

  //double l2_norm();

  /** Scalar assignment */

  //Vector & operator = (const double scalar);

  /** Scalar multiplication */

  //Vector & operator *= (const double scalar);

  /** Vector assignment */

  //Vector & operator = (const Vector &vector)

  /** Vector addition */

  //Vector & operator += (const Vector &vector);

private:

  Teuchos::RCP<Epetra_FEVector> epetra_fevector_ptr;

};


}



#endif
