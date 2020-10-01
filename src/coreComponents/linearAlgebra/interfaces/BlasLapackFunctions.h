/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BlasLapackFunctions.h
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_

/// This macro provide a flexible interface for Fortran naming convetion for compiled objects
#ifdef FORTRAN_MANGLE_NO_UNDERSCORE
#define FORTRAN_MANGLE( name ) name
#else
#define FORTRAN_MANGLE( name ) name ## _
#endif

#ifdef __cplusplus
extern "C"
{

#define GEOSX_dasum FORTRAN_MANGLE( dasum )
double GEOSX_dasum( int const * n,
                    double const * dx,
                    int const * incx );

#define GEOSX_daxpy FORTRAN_MANGLE( daxpy )
void GEOSX_daxpy( int const * n,
                  double const * da,
                  double const * dx,
                  int const * incx,
                  double * dy,
                  int const * incy );

#define GEOSX_dcopy FORTRAN_MANGLE( dcopy )
void GEOSX_dcopy( int const * n,
                  double const * dx,
                  int const * incx,
                  double * dy,
                  int const * incy );

#define GEOSX_ddot FORTRAN_MANGLE( ddot )
double GEOSX_ddot( int const * n,
                   double const * dx,
                   int const * incx,
                   double const * dy,
                   int const * incy );

#define GEOSX_idamax FORTRAN_MANGLE( idamax )
int GEOSX_idamax( int const * n,
                  double const * dx,
                  int const * incx );

#define GEOSX_dgemm FORTRAN_MANGLE( dgemm )
void GEOSX_dgemm( char const * transa,
                  char const * transb,
                  int const * m,
                  int const * n,
                  int const * k,
                  double const * alpha,
                  double const * a,
                  int const * lda,
                  double const * b,
                  int const * ldb,
                  double const * beta,
                  double * c,
                  int const * ldc );

#define GEOSX_dgetrf FORTRAN_MANGLE( dgetrf )
void GEOSX_dgetrf( int const * m,
                   int const * n,
                   double * a,
                   int const * lda,
                   int * ipiv,
                   int * info );

#define GEOSX_dgetri FORTRAN_MANGLE( dgetri )
void GEOSX_dgetri( int const * n,
                   double * a,
                   int const * lda,
                   int const * ipiv,
                   double * work,
                   int const * lwork,
                   int * info );

#define GEOSX_dlange FORTRAN_MANGLE( dlange )
double GEOSX_dlange( char const * norm,
                     int const * m,
                     int const * n,
                     double const * a,
                     int const * lda,
                     double * work );

#define GEOSX_dlarnv FORTRAN_MANGLE( dlarnv )
void GEOSX_dlarnv( int const * idist,
                   int * iseed,
                   int const * n,
                   double * x );

#define GEOSX_dnrm2 FORTRAN_MANGLE( dnrm2 )
double GEOSX_dnrm2( int const * n,
                    double const * x,
                    int const * incx );

#define GEOSX_dscal FORTRAN_MANGLE( dscal )
void GEOSX_dscal( int const * n,
                  double const * da,
                  double * dx,
                  int const * incx );

#define GEOSX_dgesvd FORTRAN_MANGLE( dgesvd )
void GEOSX_dgesvd( char const * jobu,
                   char const * jobvt,
                   int const * m,
                   int const * n,
                   double * a,
                   int const * lda,
                   double * s,
                   double * u,
                   int const * ldu,
                   double * vt,
                   int const * ldvt,
                   double * wkopt,
                   int const * lwork,
                   int * info );

}
#endif

#endif //GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_
