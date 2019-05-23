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
 * @file BlasLapackFunctions.hpp
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKFUNCTIONS_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKFUNCTIONS_HPP_

double dasum_( int const * N,
               double const * DX,
               int const * INCX );
void daxpy_( int const * N,
             double const * DA,
             double const * DX,
             int const * INCX,
             double * DY,
             int const * INCY );
void dcopy_( int const * N,
             double const * DX,
             int const * INCX,
             double * DY,
             int const * INCY );
double ddot_( int const * N,
              double const * DX,
              int const * INCX,
              double const * DY,
              int const * INCY );
int idamax_( int const * N,
             double const * DX,
             int const * INCX );
void dgemm_( char const * TRANSA,
             char const * TRANSB,
             int const * M,
             int const * N,
             int const * K,
             double const * ALPHA,
             double const * A,
             int const * LDA,
             double const * B,
             int const * LDB,
             double const * BETA,
             double * C,
             int const * LDC );
void dgetrf_( int const * M,
              int const * N,
              double * A,
              int const * LDA,
              int * IPIV,
              int * INFO );
void dgetri_( int const * N,
              double * A,
              int const * LDA,
              int const * IPIV,
              double * WORK,
              int const * LWORK,
              int * INFO );
double dlange_( char const * NORM,
                int const * M,
                int const * N,
                double const * A,
                int const * LDA,
                double * WORK );
double dnrm2_( int const * N,
               double const * X,
               int const * INCX );
void dscal_( int const * N,
             double const * DA,
             double * DX,
             int const * INCX );

#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKFUNCTIONS_HPP_ */
