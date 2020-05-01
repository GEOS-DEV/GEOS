/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @brief This file contains the definition of the R2TensorT class
 * @file R2TensorT.h
 */

#ifndef R2_TENSOR_T_H_
#define R2_TENSOR_T_H_

#include "TensorBaseT.h"
//#include "TensorT.h"

template< int T_dim > class R2SymTensorT;
template< int T_dim > class R1TensorT;
template< int T_dim > class R4minSymTensorT;
template< int T_dim > class R6minSymTensorT;

/**
 * @brief R2TensorT is a rank-2 tensor object type
 * @tparam T_dim length of tensor index
 *
 * R2TensorT derives from TensorBaseT, and defines basic operations that can be
 * done on a rank-2 tensor, as well as operations that result in a rank-2
 * tensor.
 */
template< int T_dim >
class R2TensorT : public TensorBaseT< T_dim *T_dim >
{
public:
  //**** CONSTRUCTORS AND DESTRUCTORS *****************************************
  /// default constructor
  GEOSX_HOST_DEVICE
  R2TensorT( void );

  /**
   * @param[in] data use for initialization of t_data
   */
  explicit R2TensorT( const realT data ): TensorBaseT< T_dim *T_dim >( data ) {}

  /// copy constructor
  R2TensorT( const R2TensorT< T_dim > & rhs ) = default;

  /// constructor initialized by raw data
  explicit R2TensorT( const realT data[T_dim*T_dim] ): TensorBaseT< T_dim *T_dim >( data ){}

  /// explicit constructors
  R2TensorT( realT Txx, realT Txy, realT Tyx, realT Tyy );  //2D only - throws error
                                                            // otherwise
  R2TensorT( realT Txx, realT Txy, realT Txz, realT Tyx, realT Tyy, realT Tyz, realT Tzx, realT Tzy, realT Tzz ); //3D
                                                                                                                  // only
                                                                                                                  // -
                                                                                                                  // throws
                                                                                                                  // error
                                                                                                                  // otherwise

  /// non-virtual destructor
  ~R2TensorT( void ) = default;

  //***** ASSIGNMENT OPERATORS ************************************************
  /// assignment of all data to an integer
  GEOSX_HOST_DEVICE
  R2TensorT< T_dim > & operator=( const int & rhs );

  /// assignment to all data to a realT
  GEOSX_HOST_DEVICE
  R2TensorT< T_dim > & operator=( const realT & rhs );

  /// assignment to another R2TensorT
  R2TensorT< T_dim > & operator=( const R2TensorT< T_dim > & rhs ) = default;

  /// assignment to another R2SymTensorT
  R2TensorT< T_dim > & operator=( const R2SymTensorT< T_dim > & rhs );

  /// add another R2SymTensorT
  R2TensorT< T_dim > & operator+=( const R2SymTensorT< T_dim > & rhs );

  GEOSX_HOST_DEVICE
  R2TensorT & operator+=( const R2TensorT< T_dim > & rhs );

  //***** ACCESS OPERATORS ****************************************************
  /// const access to data
  GEOSX_HOST_DEVICE
  const realT & operator()( const int i, const int j ) const;

  /// non-const access to data
  GEOSX_HOST_DEVICE
  realT & operator()( const int i, const int j );


  //***** MULTIPLICATION OPERATIONS *******************************************

  /// multiply (inner product) Rank2 tensor with Rank 2 tensor
  GEOSX_HOST_DEVICE
  void AijBjk( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B );

  /// multiply (inner product) Rank2 tensor with Rank 2 tensor
  GEOSX_HOST_DEVICE
  void AijBkj( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B );

  /// multiply (inner product) Rank2 tensor with Rank 2 tensor
  void AjiBjk( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B );

  /// multiply (inner product) Rank2 tensor with Rank 2 tensor
  void AjiBkj( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B );

  /// multiply (inner product) Symmetric Rank2 tensor with Rank 2 tensor
  void AijBjk( const R2SymTensorT< T_dim > & A, const R2TensorT< T_dim > & B );

  /// multiply (inner product) Symmetric Rank2 tensor with Rank 2 tensor
  GEOSX_HOST_DEVICE
  void AijBkj( const R2SymTensorT< T_dim > & A, const R2TensorT< T_dim > & B );

  /// multiply (inner product) Rank2 tensor with Symmetric Rank 2 tensor
  void AijBjk( const R2TensorT< T_dim > & A, const R2SymTensorT< T_dim > & B );

  /// multiply (inner product) Rank2 tensor with Symmetric Rank 2 tensor
  void AjiBjk( const R2TensorT< T_dim > & A, const R2SymTensorT< T_dim > & B );

  /// multiply (dyadic product) Rank1 tensor with Rank 1 tensor
  GEOSX_HOST_DEVICE
  void dyadic_ab( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b );

  /// multiply (dyadic product) Rank1 tensor with itself
  void dyadic_aa( const R1TensorT< T_dim > & a );

  /// multiply (dyadic product) Rank1 tensor with Rank 1 tensor and add
  GEOSX_HOST_DEVICE
  void plus_dyadic_ab( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b );

  void FillColumn( const int col, const R1TensorT< T_dim > & a );
  void AddToColumn( const int col, const R1TensorT< T_dim > & a );

  void FillRow( const int row, const R1TensorT< T_dim > & a );
  void AddToRow( const int row, const R1TensorT< T_dim > & a );

  /*
     /// junk
     void AB_plus_CD( const R2TensorT<T_dim>& A,
                   const R2TensorT<T_dim>& B,
                   const R2TensorT<T_dim>& C,
                   const R2TensorT<T_dim>& D )
     {
     for( int i=0 ; i<(T_dim*T_dim) ; ++i )
      this->t_data[i] += A.t_data[i] * B.t_data[i] + C.t_data[i] * D.t_data[i] ;
     }
   */
  //****** SELF TENSOR OPERATIONS **********************************************
  /// Inner Product
  realT Inner( void ) const;

  /// Trace
  GEOSX_HOST_DEVICE
  realT Trace( void ) const;

  /// Determinant
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  realT Det( void ) const;

  /// One minus Determinant
  realT OneMinusDet( void ) const;

  /// Inverse
  realT Inverse( void ) { return Inverse( *this ); }

  /// Inverse
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  realT Inverse( R2TensorT< T_dim > & tensor );

  /// Inverse of tensor minus identity
  R2TensorT< T_dim > & Inverse_I( void ) { return Inverse_I( *this ); }

  /// Inverse of tensor minus identity
  R2TensorT< T_dim > & Inverse_I( R2TensorT< T_dim > & tensor );

  /// add a realT to the diagonal
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void PlusIdentity( const realT rhs )
  {
    this->t_data[0] += rhs;
    this->t_data[4] += rhs;
    this->t_data[8] += rhs;
  }

  //****** MATRIX MANIPULATIONS ***********************************************
  /// Swap Rows of a matrix
  void RowSwap( const int i1, const int i2 );

  /// Pivot rows to a row-reduced echelon form
  void RREF( void );

  //***** OUTPUT **************************************************************
  /// print function
  void print( std::ostream & os ) const;

  //****** TENSOR TRSNSFORMATION **********************************************

  void Aijkl_to_Bmn( const R4minSymTensorT< 3 > & A );

  //***** FRIEND DECLARATIONS *************************************************
  /// declare R2SymTensorT a friend so that it can access t_data directly
  friend class R2SymTensorT< T_dim >;

  /// declare R2TensorT a friend so that it can access t_data directly
  friend class R1TensorT< T_dim >;
  friend class R4minSymTensorT< T_dim >;
  friend class R4minSymTensorT< 3 >;
  friend class R6minSymTensorT< T_dim >;

private:

};
//*****************************************************************************
//***** END DECLARATIONS ******************************************************
//*****************************************************************************



template< int T_dim >
void R2TensorT< T_dim >::print( std::ostream & os ) const
{
  for( int i=0; i<T_dim; ++i )
    for( int j=0; j<T_dim; ++j )
      os<<(*this)( i, j )<<'\t';
}



#include "R2SymTensorT.h"
#include "R1TensorT.h"
#include "R4minSymTensorT.h"

//*****************************************************************************
//***** R2TensorT Member Function Definition **********************************
//*****************************************************************************

//**** CONSTRUCTORS AND DESTRUCTORS *******************************************
/**
 * @return none
 */
template< int T_dim >
GEOSX_HOST_DEVICE
R2TensorT< T_dim >::R2TensorT( void ):
  TensorBaseT< T_dim *T_dim >()
{}

/// Explicit 2D constructor
///
/// Template specialisation - if templated on another dimension constructor will
// throw a compile time error.
template<>
inline R2TensorT< 2 >::R2TensorT( realT Txx, realT Txy,
                                  realT Tyx, realT Tyy ):
  TensorBaseT< 2*2 >()
{
  this->t_data[0] = Txx;
  this->t_data[1] = Txy;
  this->t_data[2] = Tyx;
  this->t_data[3] = Tyy;
}


/// Explicit 3D constructor
///
/// Template specialisation - if templated on another dimension constructor will
// throw a compile time error.
template<>
inline R2TensorT< 3 >::R2TensorT( realT Txx, realT Txy, realT Txz,
                                  realT Tyx, realT Tyy, realT Tyz,
                                  realT Tzx, realT Tzy, realT Tzz ):
  TensorBaseT< 3*3 >()
{
  this->t_data[0] = Txx;
  this->t_data[1] = Txy;
  this->t_data[2] = Txz;
  this->t_data[3] = Tyx;
  this->t_data[4] = Tyy;
  this->t_data[5] = Tyz;
  this->t_data[6] = Tzx;
  this->t_data[7] = Tzy;
  this->t_data[8] = Tzz;
}

//***** ACCESS OPERATORS ******************************************************
/**
 * @param[in] i first index of the data to be returned
 * @param[in] j second index of the data to be returned
 * @return non-modifiable reference to the tensor data at index (i,j)
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
const realT & R2TensorT< T_dim >::operator()( const int i, const int j ) const
{
  return (this->t_data[ i*T_dim + j ]);
}

/**
 * @param[in] i first index of the data to be returned
 * @param[in] j second index of the data to be returned
 * @return modifiable reference to the tensor data at index (i,j)
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
realT & R2TensorT< T_dim >::operator()( const int i, const int j )
{
  return (this->t_data[ i*T_dim + j ]);
}


//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @param[in] rhs value to set each member of t_data to
 * @return reference to this
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
R2TensorT< T_dim > & R2TensorT< T_dim >::operator=( const int & rhs )
{
  TensorBaseT< T_dim *T_dim >::operator=( rhs );
  return *this;
}

/**
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template< int T_dim >
GEOSX_HOST_DEVICE
inline R2TensorT< T_dim > & R2TensorT< T_dim >::operator=( const realT & rhs )
{
  TensorBaseT< T_dim *T_dim >::operator=( rhs );
  return *this;
}


/**
 * @param[in] rhs symmetic tensor to copy
 * @return reference to *this
 */
template< int T_dim >
inline R2TensorT< T_dim > & R2TensorT< T_dim >::operator=( const R2SymTensorT< T_dim > & rhs )
{
  for( int i=1; i<=T_dim; ++i )
    for( int j=1; j<=T_dim; ++j )
      (*this)( i, j ) = rhs( i, j );

  return *this;
}
template<>
inline R2TensorT< 2 > & R2TensorT< 2 >::operator=( const R2SymTensorT< 2 > & rhs )
{
  t_data[0] = rhs.t_data[0];
  t_data[1] = rhs.t_data[1];
  t_data[2] = rhs.t_data[1];
  t_data[3] = rhs.t_data[2];

  return *this;
}
template<>
inline R2TensorT< 3 > & R2TensorT< 3 >::operator=( const R2SymTensorT< 3 > & rhs )
{
  t_data[0] = rhs.t_data[0];
  t_data[1] = rhs.t_data[1];
  t_data[2] = rhs.t_data[3];

  t_data[3] = rhs.t_data[1];
  t_data[4] = rhs.t_data[2];
  t_data[5] = rhs.t_data[4];

  t_data[6] = rhs.t_data[3];
  t_data[7] = rhs.t_data[4];
  t_data[8] = rhs.t_data[5];

  return *this;
}


template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
R2TensorT< T_dim > & R2TensorT< T_dim >::operator+=( const R2TensorT< T_dim > & rhs )
{
  TensorBaseT< T_dim *T_dim >::operator+=( rhs );
  return *this;
}


/**
 * @param[in] rhs symmetic tensor to add
 * @return reference to *this
 */
template< int T_dim >
inline R2TensorT< T_dim > & R2TensorT< T_dim >::operator+=( const R2SymTensorT< T_dim > & rhs )
{
  for( int i=1; i<=T_dim; ++i )
    for( int j=1; j<=T_dim; ++j )
      (*this)( i, j ) += rhs( i, j );

  return *this;
}

//***** MULTIPLICATION OPERATORS **********************************************
/**
 * @param[in] A rank-2 tensor
 * @param[in] B rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {AB}\f$ -or-
 *\f$A_{ij} B_{jk}\f$
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void R2TensorT< T_dim >::AijBjk( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B )
{
  // int ij;
  // int jk;
  // int ik;

//  if( T_dim == 2 )
//  {
//    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[2];
//    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[1]*B.t_data[3];
//
//    this->t_data[2] = A.t_data[2]*B.t_data[0] + A.t_data[3]*B.t_data[2];
//    this->t_data[3] = A.t_data[2]*B.t_data[1] + A.t_data[3]*B.t_data[3];
//  }
//  else if( T_dim == 3 )
//  {
  this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[3] + A.t_data[2]*B.t_data[6];
  this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[1]*B.t_data[4] + A.t_data[2]*B.t_data[7];
  this->t_data[2] = A.t_data[0]*B.t_data[2] + A.t_data[1]*B.t_data[5] + A.t_data[2]*B.t_data[8];

  this->t_data[3] = A.t_data[3]*B.t_data[0] + A.t_data[4]*B.t_data[3] + A.t_data[5]*B.t_data[6];
  this->t_data[4] = A.t_data[3]*B.t_data[1] + A.t_data[4]*B.t_data[4] + A.t_data[5]*B.t_data[7];
  this->t_data[5] = A.t_data[3]*B.t_data[2] + A.t_data[4]*B.t_data[5] + A.t_data[5]*B.t_data[8];

  this->t_data[6] = A.t_data[6]*B.t_data[0] + A.t_data[7]*B.t_data[3] + A.t_data[8]*B.t_data[6];
  this->t_data[7] = A.t_data[6]*B.t_data[1] + A.t_data[7]*B.t_data[4] + A.t_data[8]*B.t_data[7];
  this->t_data[8] = A.t_data[6]*B.t_data[2] + A.t_data[7]*B.t_data[5] + A.t_data[8]*B.t_data[8];
//  }
//  else
//  {
//    for( int i=1 ; i<=T_dim ; ++i )
//      for( int k=1 ; k<=T_dim ; ++k )
//      {
//        ik =  (i-1)*T_dim + (k-1);
//        this->t_data[ik] = 0.0;
//        for( int j=1 ; j<=T_dim ; ++j )
//        {
//          ij = (i-1)*T_dim + (j-1);
//          jk = (j-1)*T_dim + (k-1);
//          this->t_data[ik] += A.t_data[ij] * B.t_data[jk];
//        }
//      }
//  }
}


/**
 * @param[in] A rank-2 tensor
 * @param[in] B rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {AB^T}\f$ -or-
 *\f$A_{ij} B_{kj}\f$
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void R2TensorT< T_dim >::AijBkj( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[1];
    this->t_data[1] = A.t_data[0]*B.t_data[2] + A.t_data[1]*B.t_data[3];

    this->t_data[2] = A.t_data[2]*B.t_data[0] + A.t_data[3]*B.t_data[1];
    this->t_data[3] = A.t_data[2]*B.t_data[2] + A.t_data[3]*B.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[1] + A.t_data[2]*B.t_data[2];
    this->t_data[1] = A.t_data[0]*B.t_data[3] + A.t_data[1]*B.t_data[4] + A.t_data[2]*B.t_data[5];
    this->t_data[2] = A.t_data[0]*B.t_data[6] + A.t_data[1]*B.t_data[7] + A.t_data[2]*B.t_data[8];

    this->t_data[3] = A.t_data[3]*B.t_data[0] + A.t_data[4]*B.t_data[1] + A.t_data[5]*B.t_data[2];
    this->t_data[4] = A.t_data[3]*B.t_data[3] + A.t_data[4]*B.t_data[4] + A.t_data[5]*B.t_data[5];
    this->t_data[5] = A.t_data[3]*B.t_data[6] + A.t_data[4]*B.t_data[7] + A.t_data[5]*B.t_data[8];

    this->t_data[6] = A.t_data[6]*B.t_data[0] + A.t_data[7]*B.t_data[1] + A.t_data[8]*B.t_data[2];
    this->t_data[7] = A.t_data[6]*B.t_data[3] + A.t_data[7]*B.t_data[4] + A.t_data[8]*B.t_data[5];
    this->t_data[8] = A.t_data[6]*B.t_data[6] + A.t_data[7]*B.t_data[7] + A.t_data[8]*B.t_data[8];
  }
  else
  {
    int ij;
    int kj;
    int ik;
    for( int i=1; i<=T_dim; ++i )
      for( int k=1; k<=T_dim; ++k )
      {
        ik =  (i-1)*T_dim + (k-1);
        this->t_data[ik] = 0.0;
        for( int j=1; j<=T_dim; ++j )
        {
          ij = (i-1)*T_dim + (j-1);
          kj = (k-1)*T_dim + (j-1);
          this->t_data[ik] += A.t_data[ij] * B.t_data[kj];
        }
      }
  }
}

/**
 * @param[in] A rank-2 tensor
 * @param[in] B rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {A^TB}\f$ -or-
 *\f$A_{ji} B_{jk}\f$
 */
template< int T_dim >
inline void R2TensorT< T_dim >::AjiBjk( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B )
{
  int ji;
  int jk;
  int ik;

  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[2]*B.t_data[2];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[2]*B.t_data[3];

    this->t_data[2] = A.t_data[1]*B.t_data[0] + A.t_data[3]*B.t_data[2];
    this->t_data[3] = A.t_data[1]*B.t_data[1] + A.t_data[3]*B.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[3]*B.t_data[3] + A.t_data[6]*B.t_data[6];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[3]*B.t_data[4] + A.t_data[6]*B.t_data[7];
    this->t_data[2] = A.t_data[0]*B.t_data[2] + A.t_data[3]*B.t_data[5] + A.t_data[6]*B.t_data[8];

    this->t_data[3] = A.t_data[1]*B.t_data[0] + A.t_data[4]*B.t_data[3] + A.t_data[7]*B.t_data[6];
    this->t_data[4] = A.t_data[1]*B.t_data[1] + A.t_data[4]*B.t_data[4] + A.t_data[7]*B.t_data[7];
    this->t_data[5] = A.t_data[1]*B.t_data[2] + A.t_data[4]*B.t_data[5] + A.t_data[7]*B.t_data[8];

    this->t_data[6] = A.t_data[2]*B.t_data[0] + A.t_data[5]*B.t_data[3] + A.t_data[8]*B.t_data[6];
    this->t_data[7] = A.t_data[2]*B.t_data[1] + A.t_data[5]*B.t_data[4] + A.t_data[8]*B.t_data[7];
    this->t_data[8] = A.t_data[2]*B.t_data[2] + A.t_data[5]*B.t_data[5] + A.t_data[8]*B.t_data[8];
  }
  else
  {
    for( int i=1; i<=T_dim; ++i )
      for( int k=1; k<=T_dim; ++k )
      {
        ik =  (i-1)*T_dim + (k-1);
        this->t_data[ik] = 0.0;
        for( int j=1; j<=T_dim; ++j )
        {
          ji = (j-1)*T_dim + (i-1);
          jk = (j-1)*T_dim + (k-1);
          this->t_data[ik] += A.t_data[ji] * B.t_data[jk];
        }
      }
  }
}

/**
 * @param[in] A rank-2 tensor
 * @param[in] B rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {A^TB^T}\f$ -or-
 *\f$A_{ji} B_{kj}\f$
 */
template< int T_dim >
inline void R2TensorT< T_dim >::AjiBkj( const R2TensorT< T_dim > & A, const R2TensorT< T_dim > & B )
{
  int ji;
  int kj;
  int ik;

  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[2]*B.t_data[1];
    this->t_data[1] = A.t_data[0]*B.t_data[2] + A.t_data[2]*B.t_data[3];

    this->t_data[2] = A.t_data[1]*B.t_data[0] + A.t_data[3]*B.t_data[1];
    this->t_data[3] = A.t_data[1]*B.t_data[2] + A.t_data[3]*B.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[3]*B.t_data[1] + A.t_data[6]*B.t_data[2];
    this->t_data[1] = A.t_data[0]*B.t_data[3] + A.t_data[3]*B.t_data[4] + A.t_data[6]*B.t_data[5];
    this->t_data[2] = A.t_data[0]*B.t_data[6] + A.t_data[3]*B.t_data[7] + A.t_data[6]*B.t_data[8];

    this->t_data[3] = A.t_data[1]*B.t_data[0] + A.t_data[4]*B.t_data[1] + A.t_data[7]*B.t_data[2];
    this->t_data[4] = A.t_data[1]*B.t_data[3] + A.t_data[4]*B.t_data[4] + A.t_data[7]*B.t_data[5];
    this->t_data[5] = A.t_data[1]*B.t_data[6] + A.t_data[4]*B.t_data[7] + A.t_data[7]*B.t_data[8];

    this->t_data[6] = A.t_data[2]*B.t_data[0] + A.t_data[5]*B.t_data[1] + A.t_data[8]*B.t_data[2];
    this->t_data[7] = A.t_data[2]*B.t_data[3] + A.t_data[5]*B.t_data[4] + A.t_data[8]*B.t_data[5];
    this->t_data[8] = A.t_data[2]*B.t_data[6] + A.t_data[5]*B.t_data[7] + A.t_data[8]*B.t_data[8];
  }
  else
  {
    for( int i=1; i<=T_dim; ++i )
      for( int k=1; k<=T_dim; ++k )
      {
        ik =  (i-1)*T_dim + (k-1);
        this->t_data[ik] = 0.0;
        for( int j=1; j<=T_dim; ++j )
        {
          ji = (j-1)*T_dim + (i-1);
          kj = (k-1)*T_dim + (j-1);
          this->t_data[ik] += A.t_data[ji] * B.t_data[kj];
        }
      }
  }
}


/**
 * @param[in] A symmetric rank-2 tensor
 * @param[in] B rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {AB}\f$ -or-
 *\f$A_{ij} B_{jk}\f$
 */
template< int T_dim >
inline void R2TensorT< T_dim >::AijBjk( const R2SymTensorT< T_dim > & A, const R2TensorT< T_dim > & B )
{

  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[2];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[1]*B.t_data[3];

    this->t_data[2] = A.t_data[1]*B.t_data[0] + A.t_data[2]*B.t_data[2];
    this->t_data[3] = A.t_data[1]*B.t_data[1] + A.t_data[2]*B.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[3] + A.t_data[3]*B.t_data[6];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[1]*B.t_data[4] + A.t_data[3]*B.t_data[7];
    this->t_data[2] = A.t_data[0]*B.t_data[2] + A.t_data[1]*B.t_data[5] + A.t_data[3]*B.t_data[8];

    this->t_data[3] = A.t_data[1]*B.t_data[0] + A.t_data[2]*B.t_data[3] + A.t_data[4]*B.t_data[6];
    this->t_data[4] = A.t_data[1]*B.t_data[1] + A.t_data[2]*B.t_data[4] + A.t_data[4]*B.t_data[7];
    this->t_data[5] = A.t_data[1]*B.t_data[2] + A.t_data[2]*B.t_data[5] + A.t_data[4]*B.t_data[8];

    this->t_data[6] = A.t_data[3]*B.t_data[0] + A.t_data[4]*B.t_data[3] + A.t_data[5]*B.t_data[6];
    this->t_data[7] = A.t_data[3]*B.t_data[1] + A.t_data[4]*B.t_data[4] + A.t_data[5]*B.t_data[7];
    this->t_data[8] = A.t_data[3]*B.t_data[2] + A.t_data[4]*B.t_data[5] + A.t_data[5]*B.t_data[8];
  }
  else
  {
    GEOSX_WARNING( "R2TensorT::AijBjk(R2SymTensorT,R2TensorT) not implemented for dimension > 3 " );
  }
}

/**
 * @param[in] A symmetric rank-2 tensor
 * @param[in] B rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {AB^T}\f$ -or-
 *\f$A_{ij} B_{kj}\f$
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void R2TensorT< T_dim >::AijBkj( const R2SymTensorT< T_dim > & A, const R2TensorT< T_dim > & B )
{
//  int ij;
//  int jk;
//  int ik;

//  if( T_dim == 2 )
//  {
//    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[1];
//    this->t_data[1] = A.t_data[0]*B.t_data[2] + A.t_data[1]*B.t_data[3];
//
//    this->t_data[2] = A.t_data[1]*B.t_data[0] + A.t_data[2]*B.t_data[1];
//    this->t_data[3] = A.t_data[1]*B.t_data[2] + A.t_data[2]*B.t_data[3];
//  }
//  else if( T_dim == 3 )
//  {
  this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[1] + A.t_data[3]*B.t_data[2];
  this->t_data[1] = A.t_data[0]*B.t_data[3] + A.t_data[1]*B.t_data[4] + A.t_data[3]*B.t_data[5];
  this->t_data[2] = A.t_data[0]*B.t_data[6] + A.t_data[1]*B.t_data[7] + A.t_data[3]*B.t_data[8];

  this->t_data[3] = A.t_data[1]*B.t_data[0] + A.t_data[2]*B.t_data[1] + A.t_data[4]*B.t_data[2];
  this->t_data[4] = A.t_data[1]*B.t_data[3] + A.t_data[2]*B.t_data[4] + A.t_data[4]*B.t_data[5];
  this->t_data[5] = A.t_data[1]*B.t_data[6] + A.t_data[2]*B.t_data[7] + A.t_data[4]*B.t_data[8];

  this->t_data[6] = A.t_data[3]*B.t_data[0] + A.t_data[4]*B.t_data[1] + A.t_data[5]*B.t_data[2];
  this->t_data[7] = A.t_data[3]*B.t_data[3] + A.t_data[4]*B.t_data[4] + A.t_data[5]*B.t_data[5];
  this->t_data[8] = A.t_data[3]*B.t_data[6] + A.t_data[4]*B.t_data[7] + A.t_data[5]*B.t_data[8];
//  }
//  else
//  {
//    GEOSX_WARNING("R2TensorT::AijBkj(R2SymTensorT,R2TensorT) not implemented for dimension > 3");
//  }
}

/**
 * @param[in] A rank-2 tensor
 * @param[in] B symmetric rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {AB}\f$ -or-
 *\f$A_{ij} B_{jk}\f$
 */
template< int T_dim >
inline void R2TensorT< T_dim >::AijBjk( const R2TensorT< T_dim > & A, const R2SymTensorT< T_dim > & B )
{

  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[1];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[1]*B.t_data[2];

    this->t_data[2] = A.t_data[2]*B.t_data[0] + A.t_data[3]*B.t_data[1];
    this->t_data[3] = A.t_data[2]*B.t_data[1] + A.t_data[3]*B.t_data[2];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[1]*B.t_data[1] + A.t_data[2]*B.t_data[3];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[1]*B.t_data[2] + A.t_data[2]*B.t_data[4];
    this->t_data[2] = A.t_data[0]*B.t_data[3] + A.t_data[1]*B.t_data[4] + A.t_data[2]*B.t_data[5];

    this->t_data[3] = A.t_data[3]*B.t_data[0] + A.t_data[4]*B.t_data[1] + A.t_data[5]*B.t_data[3];
    this->t_data[4] = A.t_data[3]*B.t_data[1] + A.t_data[4]*B.t_data[2] + A.t_data[5]*B.t_data[4];
    this->t_data[5] = A.t_data[3]*B.t_data[3] + A.t_data[4]*B.t_data[4] + A.t_data[5]*B.t_data[5];

    this->t_data[6] = A.t_data[6]*B.t_data[0] + A.t_data[7]*B.t_data[1] + A.t_data[8]*B.t_data[3];
    this->t_data[7] = A.t_data[6]*B.t_data[1] + A.t_data[7]*B.t_data[2] + A.t_data[8]*B.t_data[4];
    this->t_data[8] = A.t_data[6]*B.t_data[3] + A.t_data[7]*B.t_data[4] + A.t_data[8]*B.t_data[5];
  }
  else
  {
    GEOSX_WARNING( "R2TensorT::AijBjk(R2TensorT,R2SymTensorT) not implemented for dimension > 3" );
  }
}

/**
 * @param[in] A rank-2 tensor
 * @param[in] B symmetric rank-2 tensor
 * @return none
 *
 * This function performs matrix multiplication \f$\mathbf {A^TB}\f$ -or-
 *\f$A_{ji} B_{jk}\f$
 */
template< int T_dim >
inline void R2TensorT< T_dim >::AjiBjk( const R2TensorT< T_dim > & A, const R2SymTensorT< T_dim > & B )
{

  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[2]*B.t_data[1];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[2]*B.t_data[2];

    this->t_data[2] = A.t_data[1]*B.t_data[0] + A.t_data[3]*B.t_data[1];
    this->t_data[3] = A.t_data[1]*B.t_data[1] + A.t_data[3]*B.t_data[2];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0]*B.t_data[0] + A.t_data[3]*B.t_data[1] + A.t_data[6]*B.t_data[3];
    this->t_data[1] = A.t_data[0]*B.t_data[1] + A.t_data[3]*B.t_data[2] + A.t_data[6]*B.t_data[4];
    this->t_data[2] = A.t_data[0]*B.t_data[3] + A.t_data[3]*B.t_data[4] + A.t_data[6]*B.t_data[5];

    this->t_data[3] = A.t_data[1]*B.t_data[0] + A.t_data[4]*B.t_data[1] + A.t_data[7]*B.t_data[3];
    this->t_data[4] = A.t_data[1]*B.t_data[1] + A.t_data[4]*B.t_data[2] + A.t_data[7]*B.t_data[4];
    this->t_data[5] = A.t_data[1]*B.t_data[3] + A.t_data[4]*B.t_data[4] + A.t_data[7]*B.t_data[5];

    this->t_data[6] = A.t_data[2]*B.t_data[0] + A.t_data[5]*B.t_data[1] + A.t_data[8]*B.t_data[3];
    this->t_data[7] = A.t_data[2]*B.t_data[1] + A.t_data[5]*B.t_data[2] + A.t_data[8]*B.t_data[4];
    this->t_data[8] = A.t_data[2]*B.t_data[3] + A.t_data[5]*B.t_data[4] + A.t_data[8]*B.t_data[5];
  }
  else
  {
    GEOSX_WARNING( "R2TensorT::AjiBjk(R2TensorT,R2SymTensorT) not implemented for dimension > 3" );
  }
}

/**
 * @param[in] a rank-1 tensor
 * @param[in] b rank-1 tensor
 * @return none
 *
 * This function performs a dyadic product of two rank-1 tensors \f$\mathbf {a
 *\otimes b}\f$ -or- \f$a_i b_j\f$
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void R2TensorT< T_dim >::dyadic_ab( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b )
{
//  if( T_dim == 2 )
//  {
//    this->t_data[0] = a.t_data[0]*b.t_data[0];
//    this->t_data[1] = a.t_data[0]*b.t_data[1];
//
//    this->t_data[2] = a.t_data[1]*b.t_data[0];
//    this->t_data[3] = a.t_data[1]*b.t_data[1];
//  }
//  else if( T_dim == 3 )
//  {
  this->t_data[0] = a.t_data[0]*b.t_data[0];
  this->t_data[1] = a.t_data[0]*b.t_data[1];
  this->t_data[2] = a.t_data[0]*b.t_data[2];

  this->t_data[3] = a.t_data[1]*b.t_data[0];
  this->t_data[4] = a.t_data[1]*b.t_data[1];
  this->t_data[5] = a.t_data[1]*b.t_data[2];

  this->t_data[6] = a.t_data[2]*b.t_data[0];
  this->t_data[7] = a.t_data[2]*b.t_data[1];
  this->t_data[8] = a.t_data[2]*b.t_data[2];
//  }
//  else
//  {
//    GEOSX_WARNING("R2TensorT::dyadic_ab(R1TensorT,R1TensorT) not implemented for dimension > 3");
//  }
}


/**
 * @param[in] a rank-1 tensor
 * @return none
 *
 * This function performs a dyadic product of a rank-1 tensor with itself
 *  \f$\mathbf {a \otimes a}\f$ -or- \f$a_i a_j\f$
 */
template< int T_dim >
inline void R2TensorT< T_dim >::dyadic_aa( const R1TensorT< T_dim > & a )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = a.t_data[0]*a.t_data[0];
    this->t_data[1] = a.t_data[0]*a.t_data[1];

    this->t_data[2] = a.t_data[1]*a.t_data[0];
    this->t_data[3] = a.t_data[1]*a.t_data[1];
  }
  else if( T_dim == 3 )
  {

    this->t_data[0] = a.t_data[0]*a.t_data[0];
    this->t_data[1] = a.t_data[0]*a.t_data[1];
    this->t_data[2] = a.t_data[0]*a.t_data[2];

    this->t_data[3] = this->t_data[1];
    this->t_data[4] = a.t_data[1]*a.t_data[1];
    this->t_data[5] = a.t_data[1]*a.t_data[2];

    this->t_data[6] = this->t_data[2];
    this->t_data[7] = this->t_data[5];
    this->t_data[8] = a.t_data[2]*a.t_data[2];
  }
  else
  {
    GEOSX_WARNING( "R2TensorT::dyadic_ab(R1TensorT,R1TensorT) not implemented for dimension > 3" );
  }
}

/**
 * @param[in] a rank-1 tensor
 * @param[in] b rank-1 tensor
 * @return none
 *
 * This function adds the dyadic product of two rank-1 tensors to this
 * tensor.
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void R2TensorT< T_dim >::plus_dyadic_ab( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b )
{
//  if( T_dim == 2 )
//  {
//    this->t_data[0] += a.t_data[0]*b.t_data[0];
//    this->t_data[1] += a.t_data[0]*b.t_data[1];
//
//    this->t_data[2] += a.t_data[1]*b.t_data[0];
//    this->t_data[3] += a.t_data[1]*b.t_data[1];
//  }
//  else if( T_dim == 3 )
//  {
  this->t_data[0] += a.t_data[0]*b.t_data[0];
  this->t_data[1] += a.t_data[0]*b.t_data[1];
  this->t_data[2] += a.t_data[0]*b.t_data[2];

  this->t_data[3] += a.t_data[1]*b.t_data[0];
  this->t_data[4] += a.t_data[1]*b.t_data[1];
  this->t_data[5] += a.t_data[1]*b.t_data[2];

  this->t_data[6] += a.t_data[2]*b.t_data[0];
  this->t_data[7] += a.t_data[2]*b.t_data[1];
  this->t_data[8] += a.t_data[2]*b.t_data[2];
//  }
//  else
//  {
//    GEOSX_WARNING("R2TensorT::dyadic_ab(R1TensorT,R1TensorT) not implemented for dimension > 3");
//  }
}

/**
 * @return trace of (*this)
 *
 * This function returns the trace of the tensor that it is called from.
 */
template< int T_dim >
GEOSX_HOST_DEVICE
inline realT R2TensorT< T_dim >::Trace( void ) const
{
  realT trace=0;

  if( T_dim==2 )
  {
    trace = this->t_data[0] + this->t_data[3];
  }
  else if( T_dim==3 )
  {
    trace = this->t_data[0] + this->t_data[4] + this->t_data[8];
  }
  else
  {
    int c=0;

    for( int ii=1; ii<=T_dim; ++ii )
    {
      trace += this->t_data[c];
      c += T_dim;
      ++c;
    }
  }
  return trace;
}


/**
 * @return determinant of (*this)
 *
 * This function returns the determinate of the tensor that it is called from.
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
realT R2TensorT< T_dim >::Det( void ) const
{
  realT det=0;
  if( T_dim == 2 )
    det = this->t_data[0]*(this->t_data[3]) - this->t_data[1]*(this->t_data[2]);
  else if( T_dim == 3 )
    det = this->t_data[0]*( this->t_data[4]*(this->t_data[8]) - this->t_data[5]*(this->t_data[7]) )
          - this->t_data[3]*( this->t_data[1]*(this->t_data[8]) - this->t_data[2]*(this->t_data[7]) )
          + this->t_data[6]*( this->t_data[1]*(this->t_data[5]) - this->t_data[2]*(this->t_data[4]) );
  return det;
}

/**
 * @return inner product of (*this) with itself
 *
 * This function returns the inner product of the tensor that it is called from
 * with itself
 */
template< int T_dim >
inline realT R2TensorT< T_dim >::Inner( void ) const
{
  realT rval=0;
  if( T_dim == 2 )
    rval = this->t_data[0]*(this->t_data[0]) + this->t_data[1]*(this->t_data[1])
           + this->t_data[2]*(this->t_data[2]) + this->t_data[3]*(this->t_data[3]);
  else if( T_dim == 3 )
    rval = this->t_data[0]*(this->t_data[0]) + this->t_data[1]*(this->t_data[1])
           + this->t_data[2]*(this->t_data[2]) + this->t_data[3]*(this->t_data[3])
           + this->t_data[4]*(this->t_data[4]) + this->t_data[5]*(this->t_data[5])
           + this->t_data[6]*(this->t_data[6]) + this->t_data[7]*(this->t_data[7])
           + this->t_data[8]*(this->t_data[8]);
  else
  {
    GEOSX_WARNING( "R2TensorT::Inner() not implemented for dimension > 3" );
  }

  return rval;
}


/**
 * @return 1-Det(*this)
 *
 * This function returns 1 - det(F) assuming that (*this) = F-I where I is the
 * identity. This
 * is useful when you have a tensor \f$\mathbf{F}\f$ close to the identity that
 * you are storing as
 * \f$\mathbf{F-I}\f$.
 */
template< int T_dim >
inline realT R2TensorT< T_dim >::OneMinusDet( void ) const
{
  realT One_Det=0;
  if( T_dim == 2 )
    One_Det = -this->t_data[0] - this->t_data[3] - this->t_data[0]*(this->t_data[3]) + this->t_data[1]*(this->t_data[2]);
  else if( T_dim == 3 )
    One_Det = -( this->t_data[0] + this->t_data[4] + this->t_data[8])
              - ( this->t_data[0]*(this->t_data[4])*(this->t_data[8]) )
              - ( this->t_data[0]*(this->t_data[4]) + this->t_data[0]*(this->t_data[8]) + this->t_data[4]*(this->t_data[8]) )
              + ( this->t_data[0] + 1.0 )*(this->t_data[5])*(this->t_data[7])
              + ( this->t_data[4] + 1.0 )*(this->t_data[2])*(this->t_data[6])
              + ( this->t_data[8] + 1.0 )*(this->t_data[1])*(this->t_data[3])
              - ( this->t_data[2]*(this->t_data[3])*(this->t_data[7]) + this->t_data[1]*(this->t_data[5])*(this->t_data[6]) );
  else
  {
    GEOSX_WARNING( "R2TensorT::Det() not implemented for dimension > 3" );
  }

  return One_Det;
}


/**
 * @return Det(a)
 *
 * This function inverts a
 */
template< int T_dim >
//R2TensorT<T_dim>& R2TensorT<T_dim>::Inverse( R2TensorT<T_dim>& a )
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
realT R2TensorT< T_dim >::Inverse( R2TensorT< T_dim > & a )
{
  realT const o1 = a.t_data[4] * a.t_data[8] - a.t_data[5] * a.t_data[7];
  realT const o2 = a.t_data[2] * a.t_data[7] - a.t_data[1] * a.t_data[8];
  realT const o3 = a.t_data[1] * a.t_data[5] - a.t_data[2] * a.t_data[4];
  realT const o4 = a.t_data[5] * a.t_data[6] - a.t_data[3] * a.t_data[8];
  realT const o5 = a.t_data[0] * a.t_data[8] - a.t_data[2] * a.t_data[6];
  realT const o6 = a.t_data[2] * a.t_data[3] - a.t_data[0] * a.t_data[5];
  realT const o7 = a.t_data[3] * a.t_data[7] - a.t_data[4] * a.t_data[6];
  realT const o8 = a.t_data[1] * a.t_data[6] - a.t_data[0] * a.t_data[7];
  realT const o9 = a.t_data[0] * a.t_data[4] - a.t_data[1] * a.t_data[3];

  realT const det = a.t_data[0] * o1 + a.t_data[3] * o2 + a.t_data[6] * o3;

  realT const tol = 1.0e-14 * a.MaxVal();
  GEOSX_ERROR_IF( std::abs( det ) <= tol, "Tolerance error. |det| = " << std::abs( det ) << ", tol = " << tol );

  realT const o11 = 1.0 / det;

  this->t_data[0] = o1 * o11;
  this->t_data[1] = o2 * o11;
  this->t_data[2] = o3 * o11;
  this->t_data[3] = o4 * o11;
  this->t_data[4] = o5 * o11;
  this->t_data[5] = o6 * o11;
  this->t_data[6] = o7 * o11;
  this->t_data[7] = o8 * o11;
  this->t_data[8] = o9 * o11;

  return det;
}


/**
 * @return *this
 *
 * This function assumes that "a" refers to a Tensor that is close to
 * the Identity, thus is stored as F-I.
 *  a = F-I
 * Assigns Value of (*this) to Inverse(F) - I
 */
template< int T_dim >
inline R2TensorT< T_dim > & R2TensorT< T_dim >::Inverse_I( R2TensorT< T_dim > & a )
{

  if( T_dim == 2 )
  {
    /* temps - incase matrix is *this */
    realT A0 = a.t_data[0];
    realT A1 = a.t_data[1];
    realT A2 = a.t_data[2];
    realT A3 = a.t_data[3];

    realT idet = 1 / (A0*A3 - A1*A2);
    const realT one_det = a.OneMinusDet();
    this->t_data[0] =  ( A3 + one_det ) * idet;
    this->t_data[1] =-A1*idet;
    this->t_data[2] =-A2*idet;
    this->t_data[3] =  ( A0 + one_det ) * idet;
  }
  else if( T_dim == 3 )
  {
    realT A11, A12, A13,
          A21, A22, A23,
          A31, A32, A33;
    const realT one_det = a.OneMinusDet();
    const realT idet = 1.0/( 1.0 - one_det );

    A11 = a.t_data[4]*a.t_data[8] - a.t_data[5]*a.t_data[7] + a.t_data[4] + a.t_data[8];
    A12 = a.t_data[2]*a.t_data[7] - a.t_data[1]*a.t_data[8] - a.t_data[1];
    A13 = a.t_data[1]*a.t_data[5] - a.t_data[2]*a.t_data[4] - a.t_data[2];
    A21 = a.t_data[5]*a.t_data[6] - a.t_data[3]*a.t_data[8] - a.t_data[3];
    A22 = a.t_data[0]*a.t_data[8] - a.t_data[2]*a.t_data[6] + a.t_data[0] + a.t_data[8];
    A23 = a.t_data[2]*a.t_data[3] - a.t_data[0]*a.t_data[5] - a.t_data[5];
    A31 = a.t_data[3]*a.t_data[7] - a.t_data[4]*a.t_data[6] - a.t_data[6];
    A32 = a.t_data[1]*a.t_data[6] - a.t_data[0]*a.t_data[7] - a.t_data[7];
    A33 = a.t_data[0]*a.t_data[4] - a.t_data[1]*a.t_data[3] + a.t_data[0] + a.t_data[4];


    this->t_data[0] = ( A11 + one_det ) * idet;
    this->t_data[1] = A12*idet;
    this->t_data[2] = A13*idet;
    this->t_data[3] = A21*idet;
    this->t_data[4] = ( A22 + one_det ) * idet;
    this->t_data[5] = A23*idet;
    this->t_data[6] = A31*idet;
    this->t_data[7] = A32*idet;
    this->t_data[8] = ( A33 + one_det ) * idet;
  }
  else
  {
    GEOSX_WARNING( "R2TensorT::Inverse( R2TensorT ) not implemented for dimension > 3" );
  }
  return *this;
}


template< int T_dim >
inline void R2TensorT< T_dim >::FillColumn( const int col, const R1TensorT< T_dim > & a )
{
  this->t_data[0+col]       = a.t_data[0];
  this->t_data[T_dim+col]   = a.t_data[1];
  this->t_data[2*T_dim+col] = a.t_data[2];
}

template< int T_dim >
inline void R2TensorT< T_dim >::AddToColumn( const int col, const R1TensorT< T_dim > & a )
{
  this->t_data[0+col]       += a.t_data[0];
  this->t_data[T_dim+col]   += a.t_data[1];
  this->t_data[2*T_dim+col] += a.t_data[2];
}


template< int T_dim >
inline void R2TensorT< T_dim >::FillRow( const int row, const R1TensorT< T_dim > & a )
{
  this->t_data[row*T_dim+1] = a.t_data[0];
  this->t_data[row*T_dim+2] = a.t_data[1];
  this->t_data[row*T_dim+3] = a.t_data[2];
}

template< int T_dim >
inline void R2TensorT< T_dim >::AddToRow( const int row, const R1TensorT< T_dim > & a )
{
  this->t_data[row*T_dim+1] += a.t_data[0];
  this->t_data[row*T_dim+2] += a.t_data[1];
  this->t_data[row*T_dim+3] += a.t_data[2];
}


template< int T_dim >
inline void R2TensorT< T_dim >::RREF( void )
{
  int pivot_row;
  realT pivot_val;

  for( int i=1; i<=T_dim; ++i )
  {
    pivot_val = 0;
    for( int ii=1; ii<=T_dim; ++ii )
    {
      int index = (ii-1)*T_dim + (i-1);
      if( this->t_data[index] > pivot_val )
      {
        pivot_val = this->t_data[index];
        pivot_row = ii;
      }
    }
//    if( i != pivot_row )
  }

}

template< int T_dim >
inline void R2TensorT< T_dim >::Aijkl_to_Bmn( const R4minSymTensorT< 3 > & A )
{
  int n_dim = 6;

  if( T_dim!=6 )
    GEOSX_WARNING( "R2TensorT<T_dim>::Aijkl_to_Bmn not implemented for T_dim /= 6" );
  else
  {
    for( int ii=0, m=0, j=0; m<n_dim; m+=++ii )
    {
      for( int i=0, c=m; c<n_dim; c+=(++i)+1+ii, ++j )
      {
        this->t_data[j*n_dim ] = A.t_data[ 0+c];
        this->t_data[j*n_dim+1] = A.t_data[12+c];
        this->t_data[j*n_dim+2] = A.t_data[30+c];
        this->t_data[j*n_dim+3] = A.t_data[ 6+c];
        this->t_data[j*n_dim+4] = A.t_data[24+c];
        this->t_data[j*n_dim+5] = A.t_data[18+c];
      }
    }
  }
}



//**** META-PROGRAMS **********************************************************

template< int T_dim >
inline void R2TensorT< T_dim >::RowSwap( const int i1, const int i2 )
{
  realT temp;
  int index1 = (i1-1) * T_dim;
  int index2 = (i2-1) * T_dim;

  for( int i=0; i<T_dim; ++i )
  {
    temp = this->t_data[index1+i];
    this->t_data[index1+i] = this->t_data[index2+i];
    this->t_data[index2+i] = temp;
  }
}


#ifdef SIG_FIG_TRUNC

template< int T_dim >
void R2TensorT< T_dim >::SetMaxVal( const R2TensorT< T_dim > & A )
{
  if( T_dim == 3 )
  {
    if( fabs( A.t_data[0] ) > (this->t_data[0]) )
      this->t_data[0] = fabs( A.t_data[0] );
    if( fabs( A.t_data[1] ) > (this->t_data[1]) )
      this->t_data[1] = fabs( A.t_data[1] );
    if( fabs( A.t_data[2] ) > (this->t_data[2]) )
      this->t_data[2] = fabs( A.t_data[2] );
    if( fabs( A.t_data[3] ) > (this->t_data[3]) )
      this->t_data[3] = fabs( A.t_data[3] );
    if( fabs( A.t_data[4] ) > (this->t_data[4]) )
      this->t_data[4] = fabs( A.t_data[4] );
    if( fabs( A.t_data[5] ) > (this->t_data[5]) )
      this->t_data[5] = fabs( A.t_data[5] );
    if( fabs( A.t_data[6] ) > (this->t_data[6]) )
      this->t_data[6] = fabs( A.t_data[6] );
    if( fabs( A.t_data[7] ) > (this->t_data[7]) )
      this->t_data[7] = fabs( A.t_data[7] );
    if( fabs( A.t_data[8] ) > (this->t_data[8]) )
      this->t_data[8] = fabs( A.t_data[8] );

  }
  else
    GEOSX_WARNING( "R1TensorT not implemented for nsdof>3" );

}
#endif

#endif
