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
 * @file BlasLapackFunctions.h
 */

#ifdef FORTRAN_MANGLE_NO_UNDERSCORE
#define FORTRAN_MANGLE(name) name
#else
#define FORTRAN_MANGLE(name) name ## _
#endif

#ifdef __cplusplus
extern "C"
{

#define GEOSX_dasum FORTRAN_MANGLE(dasum)
double GEOSX_dasum( int const * N,
                    double const * DX,
                    int const * INCX );

#define GEOSX_daxpy FORTRAN_MANGLE(daxpy)
void GEOSX_daxpy( int const * N,
                  double const * DA,
                  double const * DX,
                  int const * INCX,
                  double * DY,
                  int const * INCY );

#define GEOSX_dcopy FORTRAN_MANGLE(dcopy)
void GEOSX_dcopy( int const * N,
                  double const * DX,
                  int const * INCX,
                  double * DY,
                  int const * INCY );

#define GEOSX_ddot FORTRAN_MANGLE(ddot)
double GEOSX_ddot( int const * N,
                   double const * DX,
                   int const * INCX,
                   double const * DY,
                   int const * INCY );

#define GEOSX_idamax FORTRAN_MANGLE(idamax)
int GEOSX_idamax( int const * N,
                  double const * DX,
                  int const * INCX );

#define GEOSX_dgemm FORTRAN_MANGLE(dgemm)
void GEOSX_dgemm( char const * TRANSA,
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

#define GEOSX_dgetrf FORTRAN_MANGLE(dgetrf)
void GEOSX_dgetrf( int const * M,
                   int const * N,
                   double * A,
                   int const * LDA,
                   int * IPIV,
                   int * INFO );

#define GEOSX_dgetri FORTRAN_MANGLE(dgetri)
void GEOSX_dgetri( int const * N,
                   double * A,
                   int const * LDA,
                   int const * IPIV,
                   double * WORK,
                   int const * LWORK,
                   int * INFO );

#define GEOSX_dlange FORTRAN_MANGLE(dlange)
double GEOSX_dlange( char const * NORM,
                     int const * M,
                     int const * N,
                     double const * A,
                     int const * LDA,
                     double * WORK );

#define GEOSX_dlarnv FORTRAN_MANGLE(dlarnv)
void GEOSX_dlarnv( int const * IDIST,
                   int * ISEED,
                   int const * N,
                   double * X );

#define GEOSX_dnrm2 FORTRAN_MANGLE(dnrm2)
double GEOSX_dnrm2( int const * N,
                    double const * X,
                    int const * INCX );

#define GEOSX_dscal FORTRAN_MANGLE(dscal)
void GEOSX_dscal( int const * N,
                  double const * DA,
                  double * DX,
                  int const * INCX );
}
#endif
