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
 * @brief This file contains the definition of the R2SymTensorT.h class
 * @file R2SymTensorT.h
 */

#ifndef R2_SYM_TENSOR_T_H_
#define R2_SYM_TENSOR_T_H_

#include "TensorBaseT.h"

// PRAGMAS to suppress somewhat incorrect out of bounds errors in xlc and icc
#ifdef __INTEL_COMPILER
#pragma warning push
#pragma warning disable 175
#endif

#ifdef __IBMC__
#pragma report(disable, "1540-2907")
#endif

/// Function to determine the size of the a symmetric rank-2 tensor
template< int N >
struct SymSize
{
  enum
  { value = N + SymSize< N - 1 >::value };
};

/// template specialization for SymSize
template<>
struct SymSize< 1 >
{
  enum { value = 1 };
};

template< int T_dim > class R2TensorT;
template< int T_dim > class R1TensorT;
template< int T_dim > class R4minSymTensorT;
template< int T_dim > class R6minSymTensorT;

/**
 * @brief R2SymTensorT is a symetic rank-2 tensor object type
 * @tparam T_dim length of tensor index
 *
 * R2symTensorT derives from TensorBaseT, and defines basic operations that can
 * be
 * done on a symmetric rank-2 tensor, as well as operations that result in a
 * symmetric rank-2 tensor.
 */
template< int T_dim >
class R2SymTensorT : public TensorBaseT< SymSize< T_dim >::value >
{
public:
  static constexpr int SIZE = SymSize< T_dim >::value;

  //**** CONSTRUCTORS AND DESTRUCTORS *****************************************
  /// default constructor
  GEOSX_HOST_DEVICE
  R2SymTensorT( void );

  /**

   * @param[in] data use for initialization of t_data
   */
  explicit R2SymTensorT( const realT data ): TensorBaseT< SIZE >( data ) {}

  /// default destructor
  ~R2SymTensorT( void ) = default;

  /// copy constructor
  R2SymTensorT( const R2SymTensorT & rhs ) = default;

  explicit R2SymTensorT( const TensorBaseT< SIZE > & rhs ): TensorBaseT< SIZE >()
  { TensorBaseT< SIZE >::operator=( rhs ); }

  template< int USD >
  GEOSX_HOST_DEVICE
  R2SymTensorT( LvArray::ArraySlice< realT const, 1, USD > const & src ):
    TensorBaseT< SIZE >()
  { *this = src; }

  template< int USD >
  GEOSX_HOST_DEVICE
  R2SymTensorT( LvArray::ArraySlice< realT, 1, USD > const & src ):
    TensorBaseT< SIZE >()
  { *this = src; }

  //***** ASSIGNMENT OPERATORS
  // **************************************************

  /// assignment of all data to an integer
  GEOSX_HOST_DEVICE
  R2SymTensorT & operator=( const int & rhs );

  /// assignment of all data to a realT
  GEOSX_HOST_DEVICE
  R2SymTensorT & operator=( const realT & rhs );

  /// assignment to another R2SymTensorT
  R2SymTensorT & operator=( const R2SymTensorT & rhs ) = default;

  template< int USD >
  GEOSX_HOST_DEVICE inline
  R2SymTensorT & operator=( LvArray::ArraySlice< realT const, 1, USD > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), SIZE );

    this->t_data[ 0 ] = src[ 0 ];
    this->t_data[ 2 ] = src[ 1 ];
    this->t_data[ 5 ] = src[ 2 ];
    this->t_data[ 4 ] = src[ 3 ];
    this->t_data[ 3 ] = src[ 4 ];
    this->t_data[ 1 ] = src[ 5 ];

    return *this;
  }

  template< int USD >
  GEOSX_HOST_DEVICE inline
  R2SymTensorT & operator=( LvArray::ArraySlice< realT, 1, USD > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), SIZE );

    this->t_data[ 0 ] = src[ 0 ];
    this->t_data[ 2 ] = src[ 1 ];
    this->t_data[ 5 ] = src[ 2 ];
    this->t_data[ 4 ] = src[ 3 ];
    this->t_data[ 3 ] = src[ 4 ];
    this->t_data[ 1 ] = src[ 5 ];

    return *this;
  }

  GEOSX_HOST_DEVICE
  R2SymTensorT & operator+=( const R2SymTensorT & rhs );

  template< int USD >
  GEOSX_HOST_DEVICE inline
  R2SymTensorT & operator+=( LvArray::ArraySlice< realT const, 1, USD > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), SIZE );

    this->t_data[ 0 ] += src[ 0 ];
    this->t_data[ 2 ] += src[ 1 ];
    this->t_data[ 5 ] += src[ 2 ];
    this->t_data[ 4 ] += src[ 3 ];
    this->t_data[ 3 ] += src[ 4 ];
    this->t_data[ 1 ] += src[ 5 ];

    return *this;
  }

  template< int USD >
  GEOSX_HOST_DEVICE inline
  R2SymTensorT & operator+=( LvArray::ArraySlice< realT, 1, USD > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), SIZE );

    this->t_data[ 0 ] += src[ 0 ];
    this->t_data[ 2 ] += src[ 1 ];
    this->t_data[ 5 ] += src[ 2 ];
    this->t_data[ 4 ] += src[ 3 ];
    this->t_data[ 3 ] += src[ 4 ];
    this->t_data[ 1 ] += src[ 5 ];

    return *this;
  }

  //***** ACCESS OPERATORS ****************************************************
  /// const access to data
  realT operator()( const int i, const int j ) const;

  /// non-const access to data
  realT & operator()( const int i, const int j );

  //***** MULTIPLICATION OPERATIONS *******************************************
  realT AijBij( const R2SymTensorT & A, const R2SymTensorT & B );
  void AijBjk( const R2SymTensorT & A, const R2SymTensorT & B );
  void AijAkj( const R2TensorT< T_dim > & A );
  void AjiAjk( const R2TensorT< T_dim > & A );
  void AijAjk( const R2SymTensorT & A );

  void AijAkj_plus_Aik_plus_Aki( const R2TensorT< T_dim > & A );
  void AjiAjk_plus_Aik_plus_Aki( const R2TensorT< T_dim > & A );
  void AijAkj_m_Aik_m_Aki( const R2TensorT< T_dim > & A );

  GEOSX_HOST_DEVICE
  void QijAjkQlk( const R2SymTensorT & A, const R2TensorT< T_dim > & Q );

  void dyadic_aa( const R1TensorT< T_dim > & a );
  void dyadic_ab_plus_ba( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b );

  void AijklBkl( const R4minSymTensorT< T_dim > & A, const R2SymTensorT & B );
  void AijklBij( const R4minSymTensorT< T_dim > & A, const R2SymTensorT & B );

  //**** Overloaded arithmetic operators
  //  ******************************************

  // Scalar product
  friend R2SymTensorT operator*( realT k, const R2SymTensorT & V ){ return R2SymTensorT( V*k ); }
  friend R2SymTensorT operator*( R2SymTensorT V, realT k ){return static_cast< R2SymTensorT >(V*=k); }

  //****** TENSOR OPERATIONS **************************************************

  realT Inner( void ) const;
  GEOSX_HOST_DEVICE
  realT Trace( void ) const;
  realT Det( void ) const;
  realT AijAij( void ) const;
  void EigenVals( realT eigenvals[T_dim], const realT limArg ) const;
  void EigenVals( realT eigenvals[T_dim] ) const { EigenVals( eigenvals, 0.0 ); } // limArg
                                                                                  // is
                                                                                  // not
                                                                                  // used
                                                                                  // -
                                                                                  // causes
                                                                                  // issues
                                                                                  // with
                                                                                  // cray
                                                                                  // compiler
                                                                                  // -
                                                                                  // SW

  //  void EigenSystem( realT eigenvals[T_dim] , R1TensorT<T_dim> v[T_dim] )
  // const;
  void EigenVecs( const realT eigenvals[T_dim], R1TensorT< T_dim > v[T_dim] ) const;
  void EigenVector( const realT eigenval, R1TensorT< T_dim > & v ) const;

  void Sqrt();
  void Pow( const realT & r );

  /// inverse of *this
  inline realT Inverse( void )
  {
    return Inverse( *this );
  }


  inline realT
  Inverse( R2SymTensorT & tensor );

  /// add identity
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void PlusIdentity( const realT rhs )
  {
//    int c = 0;
//    for (int ii = 0 ; ii < T_dim ; c+=(++ii)+1)
//      this->t_data[c] += rhs;
//
    this->t_data[0] += rhs;
    this->t_data[2] += rhs;
    this->t_data[5] += rhs;
  }

  void
  print( std::ostream & os ) const;

  friend class R2TensorT< T_dim >;
  friend class R1TensorT< T_dim >;
  friend class R4minSymTensorT< T_dim >;
  friend class R6minSymTensorT< T_dim >;


private:
//  R2SymTensorT(R2SymTensorT< T_dim >&);

//  static constexpr int map[SIZE] = {{0,1,2,0,1,2,0,1,2}};


};

/**
 * @param os output stream
 *
 * This function prints the contents of the tensor
 */
template< int T_dim >
void R2SymTensorT< T_dim >::print( std::ostream & os ) const
{
  //  if( ouput_format_flag == 0 )
  //    R2TensorBaseT< T_dim ,
  //                   SymSize<T_dim>::value ,
  //                   R2SymTensorT<T_dim> >::print(os);
  //  else if( ouput_format_flag == 1 )
  //  {
  for( int i = 1; i <= T_dim; ++i )
    os << (*this)( i, i ) << '\t';

  for( int i = 1; i <= T_dim; ++i )
  {
    for( int j = i; j <= T_dim; ++j )
      if( !d_ij( i, j ))
        os << (*this)( i, j ) << '\t';
  }
  //  }
}

template< int T_dim >
GEOSX_HOST_DEVICE
R2SymTensorT< T_dim >::R2SymTensorT( void ):
  TensorBaseT< SIZE >()
{}

//***** ASSIGNMENT OPERATORS **************************************************

// Assigns all components to an integer
template< int T_dim >
GEOSX_HOST_DEVICE
R2SymTensorT< T_dim > &
R2SymTensorT< T_dim >::operator=( const int & rhs )
{
  TensorBaseT< SIZE >::operator=( rhs );
  return *this;
}

// Assigns all components to a realT
template< int T_dim >
GEOSX_HOST_DEVICE
R2SymTensorT< T_dim > &
R2SymTensorT< T_dim >::operator=( const realT & rhs )
{
  TensorBaseT< SIZE >::operator=( rhs );
  return *this;
}


template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
R2SymTensorT< T_dim > & R2SymTensorT< T_dim >::operator+=( const R2SymTensorT< T_dim > & rhs )
{
  TensorBaseT< SIZE >::operator+=( rhs );
  return *this;
}

template< int T_dim >
inline realT R2SymTensorT< T_dim >::operator()( const int i, const int j ) const
{
  int i_sym = i;
  int j_sym = j;

  if( j > i )
  {
    i_sym = j;
    j_sym = i;
  }

  int index = 0;
  for( int ii = 1; ii <= i_sym; ++ii )
    index += ii;

  index += j_sym;

  return this->t_data[index];
}

template<>
inline realT R2SymTensorT< 2 >::operator()( const int i, const int j ) const
{
  static constexpr int map[2][2] = {
    {0, 1},
    {1, 2}
  };

  return this->t_data[map[i][j]];
}

template<>
inline realT R2SymTensorT< 3 >::operator()( const int i, const int j ) const
{
  static constexpr int map[3][3] = {
    {0, 1, 3},
    {1, 2, 4},
    {3, 4, 5}
  };

  return this->t_data[map[i][j]];
}

template< int T_dim >
inline realT & R2SymTensorT< T_dim >::operator()( const int i, const int j )
{
  int i_sym = i;
  int j_sym = j;

  if( j > i )
  {
    i_sym = j;
    j_sym = i;
  }

  int index = 0;
  if( T_dim==3 )
    index = i_sym == 2 ? 3 : i_sym;
  else
    for( int ii = 1; ii <= i_sym; ++ii )
      index += ii;
  index += j_sym;

  return this->t_data[index];
}

/**
 * @return trace of (*this)
 *
 * This function returns the trace of the tensor that it is called from.
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
realT R2SymTensorT< T_dim >::Trace( void ) const
{
  realT trace = 0;
  int c = 0;

  for( int ii = 0; ii < T_dim; c+=(++ii)+1 )
  {
    trace += this->t_data[c];
  }

  return trace;
}

template<>
inline realT R2SymTensorT< 3 >::Trace( void ) const
{ return this->t_data[0] + this->t_data[2] + this->t_data[5]; }


/**
 * @return determinant of (*this)
 *
 * This function returns the determinate of the tensor that it is called from.
 */
template< int T_dim >
inline realT R2SymTensorT< T_dim >::Det( void ) const
{
  realT det = 0;
  if( T_dim == 2 )
    det = this->t_data[0] * (this->t_data[2]) - this->t_data[1] * (this->t_data[1]);
  else if( T_dim == 3 )
    det = this->t_data[0] * (this->t_data[2] * (this->t_data[5]) - this->t_data[4] * (this->t_data[4])) - this->t_data[1] * (this->t_data[1]
                                                                                                                             * (this->t_data[5]) -
                                                                                                                             this->t_data[3] *
                                                                                                                             (this->t_data[4])) +
          this->t_data[3] * (this->t_data[1] * (this->t_data[4]) - this->t_data[2]
                             * (
                               this
                                 ->
                                 t_data[3]));
  else
  {
    GEOSX_WARNING( "R2TensorT::Det() not implemented for dimension > 3" );
  }

  return det;
}

template< int T_dim >
inline realT R2SymTensorT< T_dim >::AijAij( void ) const
{
  int n_dim =  SymSize< T_dim >::value;
  realT results = 0;
  for( int i=0; i<n_dim; i++ )
  {
    results += this->t_data[i]*this->t_data[i];
  }
  return results;
}

/**
 * @return inner product of (*this) with itself
 *
 * This function returns the inner product of the tensor that it is called from
 * with itself
 */
template< int T_dim >
inline realT R2SymTensorT< T_dim >::Inner( void ) const
{
  realT rval = 0;
  if( T_dim == 2 )
    rval = this->t_data[0] * (this->t_data[0]) + 2 * (this->t_data[1]) * (this->t_data[1]) + this->t_data[2] * (this->t_data[2]);
  else if( T_dim == 3 )
    rval = this->t_data[0] * (this->t_data[0]) + 2 * (this->t_data[1]) * (this->t_data[1]) + this->t_data[2] * (this->t_data[2]) + 2 * (this->t_data[3])
           * (this->t_data[3]) + 2 * (this->t_data[4]) * (this->t_data[4]) + this->t_data[5] * (this->t_data[5]);
  else
  {
    GEOSX_WARNING( "R2TensorT::Inner() not implemented for dimension > 3" );
  }

  return rval;
}

/**
 * @param tensor tensor to invert
 * @return determinant of tensor
 */
template< int T_dim >
realT R2SymTensorT< T_dim >::Inverse( R2SymTensorT< T_dim > & tensor )
{
  realT det; // o10;

  if( T_dim == 2 )
  {
    // temps - incase matrix is *this
    const realT A0 = tensor.t_data[0];
    const realT & A1 = tensor.t_data[1];
    const realT & A2 = tensor.t_data[2];

    det = A0*A2-A1*A1;

    realT idet = 1/det;   // 1 / (A0*A2 - A1*A1);

    this->t_data[0] =  A2*idet;
    this->t_data[1] = -A1*idet;
    this->t_data[2] =  A0*idet;
  }
  else if( T_dim == 3 )
  {
#ifdef __INTEL_COMPILER
#pragma warning push
#pragma warning disable 175
#endif
#ifdef  __IBMC__
#pragma report(disable, "1540-2907")
#endif


    const realT o1 = tensor.t_data[2]*tensor.t_data[5] - tensor.t_data[4]*tensor.t_data[4];
    const realT o2 = tensor.t_data[3]*tensor.t_data[4] - tensor.t_data[1]*tensor.t_data[5];
    const realT o3 = tensor.t_data[0]*tensor.t_data[5] - tensor.t_data[3]*tensor.t_data[3];
    const realT o4 = tensor.t_data[1]*tensor.t_data[4] - tensor.t_data[2]*tensor.t_data[3];
    const realT o5 = tensor.t_data[1]*tensor.t_data[3] - tensor.t_data[0]*tensor.t_data[4];
    const realT o6 = tensor.t_data[0]*tensor.t_data[2] - tensor.t_data[1]*tensor.t_data[1];

    det = tensor.Det();
    const realT o11 = 1.0/det;

    this->t_data[0] = o1*o11;
    this->t_data[1] = o2*o11;
    this->t_data[2] = o3*o11;
    this->t_data[3] = o4*o11;
    this->t_data[4] = o5*o11;
    this->t_data[5] = o6*o11;
#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
#ifdef  __IBMC__
#pragma report(enable, "1540-2907")
#endif

  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::Inverse( R2TensorT ) not implemented for dimension > 3" );
  }
  return det;
}

/**
 * @param[out] eigenvals naked array that holds the Eigenvalues
 *
 * This function calculates the eigenvalues of *this.
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::EigenVals( realT eigenvals[T_dim],
                                              const realT /*limArg*/ ) const
{
  if( T_dim != 3 )
  {
    GEOSX_WARNING( "R2SymTensorT::EigenVals not implemented for dimension != 3" );
    return;
  }

#ifdef __INTEL_COMPILER
#pragma warning push
#pragma warning disable 175
#endif
#ifdef  __IBMC__
#pragma report(disable, "1540-2907")
#endif

  const realT Inv3 = 1 / 3.0;
  const realT sigma = Inv3 * Trace();
  const realT scale = this->MaxVal();
  if( !(scale > 0.0 ) )
  {
    for( int i =0; i< T_dim; i++ )
      eigenvals[i]=0;
    return;
  }
  const realT iscale = 1 / scale;
  const realT A11 = (this->t_data[0] - sigma) * iscale;
  const realT A22 = (this->t_data[2] - sigma) * iscale;
  const realT A33 = (this->t_data[5] - sigma) * iscale;
  const realT A23 = (this->t_data[4]) * iscale;
  const realT A13 = (this->t_data[3]) * iscale;
  const realT A12 = (this->t_data[1]) * iscale;

  const realT A12A12 = A12 * A12;
  const realT A13A13 = A13 * A13;
  const realT A23A23 = A23 * A23;
  const realT A12A13A23 = A12 * A13 * A23;

  const realT II = 0.5 * (A11 * A11 + A22 * A22 + A33 * A33) + (A12A12 + A13A13 + A23A23);

  const realT III = A11 * (A22 * A33 - A23A23) - (A12A12 * A33 - A12A13A23) + (A12A13A23 - A13A13 * A22);

  //  realT max_eigen=0;
  realT theta[T_dim];

  if( II > 0.0 )
  {
    realT temp = 3.0 / II;
    const realT sqrt_3divII = sqrt( temp );
    realT arg = 0.5 * III * temp * sqrt_3divII;

    if( arg > 1.0 )
      arg = 1.0;
    if( arg < -1.0 )
      arg = -1.0;

    temp = 2 * M_PI * Inv3;
    theta[0] = acos( arg ) * Inv3;
    theta[1] = theta[0] + temp;
    theta[2] = theta[0] - temp;

    temp = 2.0 / sqrt_3divII * scale;
    for( int i = 0; i < T_dim; ++i )
      eigenvals[i] = temp * cos( theta[i] ) + sigma;
  }
  else
    for( int i = 0; i < T_dim; ++i )
      eigenvals[i] = sigma;

  for( int i = 0; i < T_dim - 1; ++i )
    for( int j = 0; j < T_dim - 1; ++j )
      if( eigenvals[j] < eigenvals[j + 1] )
      {
        realT temp = eigenvals[j];
        eigenvals[j] = eigenvals[j + 1];
        eigenvals[j + 1] = temp;
      }
#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
#ifdef  __IBMC__
#pragma report(enable, "1540-2907")
#endif
}


/**
 * @param[in] lambda naked array that holds the Eigenvalues
 * @param[out] v a naked array of R1TensorT's that hold the eigenvectors
 *
 * This function calculates the eigenvector of *this for each eigenvalue
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::EigenVecs( const realT lambda[T_dim], R1TensorT< T_dim > v[T_dim] ) const
{
  if( T_dim != 3 )
  {
    GEOSX_WARNING( "R2SymTensorT::EigenVecs not implemented for dimension != 3" );
    return;
  }

  v[0] = 0;
  v[1] = 0;
  v[2] = 0;

  const int num_digits_4equal = 14;
  // compare eigenvalues for multiplicity
  // due to the fact that eigenvalues must be ordered
  // either e1=e2, e2=e3 or e1=e2=e3
  int multiplicity_flag = 0;
  // if e1=e2 then multiplicity flag = 1;
  if( is_equal( lambda[0], lambda[1], num_digits_4equal ))
    multiplicity_flag = 1;

  // if e2=e3 then multiplicity flag += 2;
  if( is_equal( lambda[1], lambda[2], num_digits_4equal ))
    multiplicity_flag += 2;

  // if all eigenvalues are the same, then eigenvectors are any 3 orthoganal
  // vectors.
  if( multiplicity_flag == 3 )
  {
    v[0]( 0 ) = 1;
    v[1]( 1 ) = 1;
    v[2]( 2 ) = 1;
  }
  else
  {
    R2SymTensorT< T_dim > M;
    R1TensorT< T_dim > vtemp[T_dim];

    realT v_mag[T_dim];
    realT maxmag;

    //    int i_start = 0;
    //    int i_end = T_dim;

    if( multiplicity_flag > 0 )
    {
      M = (*this);
      M.PlusIdentity( -lambda[1] );

      maxmag = 0.0;
      int imax = -1;
      for( int i = 0; i < SIZE; ++i )
        if( maxmag < fabs( M.t_data[i] ))
        {
          maxmag = fabs( M.t_data[i] );
          imax = i;
        }

      if( imax == 0 || imax == 1 )
      {
        vtemp[0].t_data[0] = -M.t_data[1];
        vtemp[0].t_data[1] = M.t_data[0];
        vtemp[0].t_data[2] = 0;

        vtemp[1].t_data[0] = -M.t_data[3] * M.t_data[0];
        vtemp[1].t_data[1] = -M.t_data[3] * M.t_data[1];
        vtemp[1].t_data[2] = M.t_data[0] * M.t_data[0] + M.t_data[1] * M.t_data[1];

        vtemp[0] /= vtemp[0].L2_Norm();
        vtemp[1] /= vtemp[1].L2_Norm();
      }
      else if( imax == 3 )
      {
        vtemp[0].t_data[0] = M.t_data[3];
        vtemp[0].t_data[1] = 0;
        vtemp[0].t_data[2] = -M.t_data[0];

        vtemp[1].t_data[0] = -M.t_data[1] * M.t_data[0];
        vtemp[1].t_data[1] = M.t_data[0] * M.t_data[0] + M.t_data[3] * M.t_data[3];
        vtemp[1].t_data[2] = -M.t_data[1] * M.t_data[3];

        vtemp[0] /= vtemp[0].L2_Norm();
        vtemp[1] /= vtemp[1].L2_Norm();
      }
      else if( imax == 2 || imax == 4 )
      {
        vtemp[0].t_data[0] = 0;
        vtemp[0].t_data[1] = -M.t_data[4];
        vtemp[0].t_data[2] = M.t_data[2];

        vtemp[1].t_data[0] = M.t_data[2] * M.t_data[2] + M.t_data[4] * M.t_data[4];
        vtemp[1].t_data[1] = -M.t_data[1] * M.t_data[2];
        vtemp[1].t_data[2] = -M.t_data[1] * M.t_data[4];

        vtemp[0] /= vtemp[0].L2_Norm();
        vtemp[1] /= vtemp[1].L2_Norm();
      }
      else if( imax == 5 )
      {
        vtemp[0].t_data[0] = 0;
        vtemp[0].t_data[1] = -M.t_data[5];
        vtemp[0].t_data[2] = M.t_data[4];

        vtemp[1].t_data[0] = M.t_data[4] * M.t_data[4] + M.t_data[5] * M.t_data[5];
        vtemp[1].t_data[1] = -M.t_data[3] * M.t_data[4];
        vtemp[1].t_data[2] = -M.t_data[3] * M.t_data[5];

        vtemp[0] /= vtemp[0].L2_Norm();
        vtemp[1] /= vtemp[1].L2_Norm();
      }

      if( multiplicity_flag == 1 )
      {
        v[0] = vtemp[0];
        v[1] = vtemp[1];
        //        i_start = 2;
        //        i_end = 2;
        v[2].Cross( v[0], v[1] );
      }
      else if( multiplicity_flag == 2 )
      {
        v[1] = vtemp[0];
        v[2] = vtemp[1];
        //        i_start = 0;
        //        i_end = 0;
        v[0].Cross( v[1], v[2] );
      }

    }
    else
    {
      for( int i = 0; i < T_dim; ++i )
      {
        maxmag = 0.0;
        M = (*this);
        M.PlusIdentity( -lambda[i] );

        vtemp[0].t_data[0] = M.t_data[2] * M.t_data[5] - M.t_data[4] * M.t_data[4];
        vtemp[0].t_data[1] = M.t_data[3] * M.t_data[4] - M.t_data[1] * M.t_data[5];
        vtemp[0].t_data[2] = M.t_data[1] * M.t_data[4] - M.t_data[3] * M.t_data[2];

        vtemp[1].t_data[0] = vtemp[0].t_data[1];
        vtemp[1].t_data[1] = M.t_data[0] * M.t_data[5] - M.t_data[3] * M.t_data[3];
        vtemp[1].t_data[2] = M.t_data[3] * M.t_data[1] - M.t_data[4] * M.t_data[0];

        vtemp[2].t_data[0] = vtemp[0].t_data[2];
        vtemp[2].t_data[1] = vtemp[1].t_data[2];
        vtemp[2].t_data[2] = M.t_data[0] * M.t_data[2] - M.t_data[1] * M.t_data[1];

        for( int j = 0; j < T_dim; ++j )
        {
          v_mag[j] = vtemp[j].L2_Norm();
          if( maxmag < v_mag[j] )
          {
            maxmag = v_mag[j];
            v[i] = vtemp[j];
            v[i] /= v_mag[j];
          }
        }
      }

      const realT TOL = pow( 10.0, -14 );
      if( fabs( Dot( v[0], v[1] ) ) > TOL )
        v[1].Cross( v[0], v[2] );

      else if( fabs( Dot( v[0], v[2] ) ) > TOL )
        v[0].Cross( v[1], v[2] );

      else if( fabs( Dot( v[1], v[2] ) ) > TOL )
        v[2].Cross( v[0], v[1] );

    }
  }
}


/**
 * @param[in] lambda the eigenvalues
 * @param[out] v R1TensorT that holds the eigenvector
 *
 * This function calculates the eigenvector of *this for an eigenvalue
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::EigenVector( const realT lambda, R1TensorT< T_dim > & v ) const
{
  if( T_dim != 3 )
  {
    GEOSX_WARNING( "R2SymTensorT::EigenVector not implemented for dimension != 3" );
    return;
  }

  R2SymTensorT< T_dim > M;
  v = 0;

  M = (*this);
  M.PlusIdentity( -lambda );

  const realT det11 = M.t_data[2] * M.t_data[5] - M.t_data[4] * M.t_data[4];
  const realT det22 = M.t_data[0] * M.t_data[5] - M.t_data[3] * M.t_data[3];
  const realT det33 = M.t_data[0] * M.t_data[2] - M.t_data[1] * M.t_data[1];

  const realT det12 = M.t_data[1] * M.t_data[5] - M.t_data[3] * M.t_data[4];
  const realT det13 = M.t_data[1] * M.t_data[4] - M.t_data[3] * M.t_data[2];
  const realT det23 = M.t_data[4] * M.t_data[0] - M.t_data[3] * M.t_data[1];

  if( fabs( det11 ) >= fabs( det22 ) && fabs( det11 ) >= fabs( det33 ))
  {
    v.t_data[0] = det11;
    v.t_data[1] = det23;
    v.t_data[2] = -det13;
  }
  else if( fabs( det22 ) >= fabs( det11 ) && fabs( det22 ) >= fabs( det33 ))
  {
    v.t_data[0] = det12;
    v.t_data[1] = det22;
    v.t_data[2] = det23;
  }
  else
  {
    v.t_data[0] = det12;
    v.t_data[1] = det22;
    v.t_data[2] = det23;
  }
}

/**
 *
 * This function takes the square root of *this
 */
template< int T_dim >
void R2SymTensorT< T_dim >::Sqrt()
{
  if( T_dim != 3 )
  {
    GEOSX_WARNING( "R2SymTensorT::Sqrt not implemented for dimension != 3" );
    return;
  }

  realT lambda[3];
  R1TensorT< T_dim > v[3];
  R2SymTensorT< T_dim > dyad;


  (*this).EigenVals( lambda );
  (*this).EigenVecs( lambda, v );

  (*this) = 0.0;
  for( int j = 0; j < 3; ++j )
  {
    lambda[j] = sqrt( lambda[j] );

    dyad.dyadic_aa( v[j] );
    dyad *= lambda[j];

    (*this) += dyad;
  }
}

/**
 * @param[in] r exponent
 *
 * This function takes power of *this
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::Pow( const realT & r )
{
  realT lambda[3];
  R1TensorT< T_dim > v[3];
  R2SymTensorT< T_dim > dyad;

  (*this).EigenVals( lambda );
  (*this).EigenVecs( lambda, v );

  (*this) = 0.0;
  for( int j = 0; j < 3; ++j )
  {
    lambda[j] = pow( lambda[j], r );

    dyad.dyadic_aa( v[j] );
    dyad *= lambda[j];

    (*this) += dyad;
  }
}


template< int T_dim >
inline realT R2SymTensorT< T_dim >::AijBij( const R2SymTensorT< T_dim > & A, const R2SymTensorT< T_dim > & B )
{
  int n_dim =  SymSize< T_dim >::value;
  realT results = 0;
  for( int i=0; i<n_dim; i++ )
  {
    results += A.t_data[i]*B.t_data[i];
  }
  return results;

}

/**
 * @param[in] A symmetric rank-2 tensor
 * @param[in] B symmetric rank-2 tensor
 *
 * This function performs matrix multiplication \f$\mathbf {AB}\f$ -or-
 *\f$A_{ij} B_{jk}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AijBjk( const R2SymTensorT< T_dim > & A, const R2SymTensorT< T_dim > & B )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1];
    this->t_data[1] = A.t_data[1] * B.t_data[0] + A.t_data[2] * B.t_data[1];
    this->t_data[2] = A.t_data[1] * B.t_data[1] + A.t_data[2] * B.t_data[2];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1] + A.t_data[3] * B.t_data[3];

    this->t_data[1] = A.t_data[1] * B.t_data[0] + A.t_data[2] * B.t_data[1] + A.t_data[4] * B.t_data[3];
    this->t_data[2] = A.t_data[1] * B.t_data[1] + A.t_data[2] * B.t_data[2] + A.t_data[4] * B.t_data[4];

    this->t_data[3] = A.t_data[3] * B.t_data[0] + A.t_data[4] * B.t_data[1] + A.t_data[5] * B.t_data[3];
    this->t_data[4] = A.t_data[3] * B.t_data[1] + A.t_data[4] * B.t_data[2] + A.t_data[5] * B.t_data[4];
    this->t_data[5] = A.t_data[3] * B.t_data[3] + A.t_data[4] * B.t_data[4] + A.t_data[5] * B.t_data[5];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AijBjk(R2TensorT) not implemented for dimension > 3 " );
  }

}

#include "R2TensorT.h"
#include "R4minSymTensorT.h"

/**
 * @param[in] A rank-2 tensor
 *
 * This function performs matrix multiplication \f$\mathbf {AA^T}\f$ -or-
 *\f$A_{ij} A_{kj}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AijAkj( const R2TensorT< T_dim > & A )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0] * A.t_data[0] + A.t_data[1] * A.t_data[1];
    this->t_data[1] = A.t_data[0] * A.t_data[2] + A.t_data[1] * A.t_data[3];
    this->t_data[2] = A.t_data[2] * A.t_data[2] + A.t_data[3] * A.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0] * A.t_data[0] + A.t_data[1] * A.t_data[1] + A.t_data[2] * A.t_data[2];

    this->t_data[1] = A.t_data[3] * A.t_data[0] + A.t_data[4] * A.t_data[1] + A.t_data[5] * A.t_data[2];
    this->t_data[2] = A.t_data[3] * A.t_data[3] + A.t_data[4] * A.t_data[4] + A.t_data[5] * A.t_data[5];

    this->t_data[3] = A.t_data[6] * A.t_data[0] + A.t_data[7] * A.t_data[1] + A.t_data[8] * A.t_data[2];
    this->t_data[4] = A.t_data[6] * A.t_data[3] + A.t_data[7] * A.t_data[4] + A.t_data[8] * A.t_data[5];
    this->t_data[5] = A.t_data[6] * A.t_data[6] + A.t_data[7] * A.t_data[7] + A.t_data[8] * A.t_data[8];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AijAkj(R2TensorT) not implemented for dimension > 3 " );
  }

}

/**
 * @param[in] A rank-2 tensor
 *
 * This function performs matrix multiplication \f$\mathbf {A^TA}\f$ -or-
 *\f$A_{ji} A_{jk}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AjiAjk( const R2TensorT< T_dim > & A )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0] * A.t_data[0] + A.t_data[2] * A.t_data[2];
    this->t_data[1] = A.t_data[0] * A.t_data[1] + A.t_data[2] * A.t_data[3];
    this->t_data[2] = A.t_data[1] * A.t_data[1] + A.t_data[3] * A.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = A.t_data[0] * A.t_data[0] + A.t_data[3] * A.t_data[3] + A.t_data[6] * A.t_data[6];

    this->t_data[1] = A.t_data[0] * A.t_data[1] + A.t_data[3] * A.t_data[4] + A.t_data[6] * A.t_data[7];
    this->t_data[2] = A.t_data[1] * A.t_data[1] + A.t_data[4] * A.t_data[4] + A.t_data[7] * A.t_data[7];

    this->t_data[3] = A.t_data[0] * A.t_data[2] + A.t_data[3] * A.t_data[5] + A.t_data[6] * A.t_data[8];
    this->t_data[4] = A.t_data[1] * A.t_data[2] + A.t_data[4] * A.t_data[5] + A.t_data[7] * A.t_data[8];
    this->t_data[5] = A.t_data[2] * A.t_data[2] + A.t_data[5] * A.t_data[5] + A.t_data[8] * A.t_data[8];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AijAkj(R2TensorT) not implemented for dimension > 3 " );
  }

}

/**
 * @param[in] A symmetric rank-2 tensor
 *
 * This function performs matrix multiplication \f$\mathbf {AA}\f$ -or-
 *\f$A_{ij} A_{jk}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AijAjk( const R2SymTensorT< T_dim > & A )
{

  if( T_dim == 2 )
  {
    this->t_data[0] = A.t_data[0] * A.t_data[0] + A.t_data[2] * A.t_data[2];
    this->t_data[1] = A.t_data[0] * A.t_data[1] + A.t_data[2] * A.t_data[3];
    this->t_data[2] = A.t_data[1] * A.t_data[1] + A.t_data[3] * A.t_data[3];
  }
  else if( T_dim == 3 )
  {
    /*    this->t_data[0] = A.t_data[0]*A.t_data[0] + A.t_data[3]*A.t_data[3] +
       A.t_data[6]*A.t_data[6];

       this->t_data[1] = A.t_data[0]*A.t_data[1] + A.t_data[3]*A.t_data[4] +
          A.t_data[6]*A.t_data[7];
       this->t_data[2] = A.t_data[1]*A.t_data[1] + A.t_data[4]*A.t_data[4] +
          A.t_data[7]*A.t_data[7];

       this->t_data[3] = A.t_data[0]*A.t_data[2] + A.t_data[3]*A.t_data[5] +
          A.t_data[6]*A.t_data[8];
       this->t_data[4] = A.t_data[1]*A.t_data[2] + A.t_data[4]*A.t_data[5] +
          A.t_data[7]*A.t_data[8];
       this->t_data[5] = A.t_data[2]*A.t_data[2] + A.t_data[5]*A.t_data[5] +
          A.t_data[8]*A.t_data[8];*/

    const realT o1 = A.t_data[1] * A.t_data[1];
    const realT o2 = A.t_data[3] * A.t_data[3];
    const realT o3 = A.t_data[4] * A.t_data[4];

    this->t_data[0] = o1 + o2 + A.t_data[0] * A.t_data[0];
    this->t_data[1] = A.t_data[0] * A.t_data[1] + A.t_data[1] * A.t_data[2] + A.t_data[3] * A.t_data[4];
    this->t_data[2] = o1 + o3 + A.t_data[2] * A.t_data[2];
    this->t_data[3] = A.t_data[0] * A.t_data[3] + A.t_data[1] * A.t_data[4] + A.t_data[3] * A.t_data[5];
    this->t_data[4] = A.t_data[1] * A.t_data[3] + A.t_data[2] * A.t_data[4] + A.t_data[4] * A.t_data[5];
    this->t_data[5] = o2 + o3 + A.t_data[5] * A.t_data[5];

  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AijAkj(R2TensorT) not implemented for dimension > 3" );
  }

}

/**
 * @param[in] A rank-2 tensor
 *
 * This function performs compound matrix operation \f$\mathbf {AA^T+A+A^T}\f$
 *-or- \f$A_{ij} A_{kj} +  A_{ik} +  A_{ki}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AijAkj_plus_Aik_plus_Aki( const R2TensorT< T_dim > & A )
{
  AijAkj( A );
  if( T_dim == 2 )
  {
    this->t_data[0] += 2 * A.t_data[0];
    this->t_data[1] += A.t_data[1] + A.t_data[2];
    this->t_data[2] += 2 * A.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] += 2 * A.t_data[0];

    this->t_data[1] += A.t_data[1] + A.t_data[3];
    this->t_data[2] += 2 * A.t_data[4];

    this->t_data[3] += A.t_data[2] + A.t_data[6];
    this->t_data[4] += A.t_data[5] + A.t_data[7];
    this->t_data[5] += 2 * A.t_data[8];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AijAkj_plus_Aik_plus_Aki not implemented for dimension > 3" );
    return;
  }
}

/**
 * @param[in] A rank-2 tensor
 *
 * This function performs compound matrix operation \f$\mathbf {A^TA+A+A^T}\f$
 *-or- \f$A_{ji} A_{jk} +  A_{ik} +  A_{ki}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AjiAjk_plus_Aik_plus_Aki( const R2TensorT< T_dim > & A )
{
  AjiAjk( A );
  if( T_dim == 2 )
  {
    this->t_data[0] += 2 * A.t_data[0];
    this->t_data[1] += A.t_data[1] + A.t_data[2];
    this->t_data[2] += 2 * A.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] += 2 * A.t_data[0];

    this->t_data[1] += A.t_data[1] + A.t_data[3];
    this->t_data[2] += 2 * A.t_data[4];

    this->t_data[3] += A.t_data[2] + A.t_data[6];
    this->t_data[4] += A.t_data[5] + A.t_data[7];
    this->t_data[5] += 2 * A.t_data[8];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AjiAjk_plus_Aik_plus_Aki not implemented for dimension > 3" );
    return;
  }
}


/**
 * @param[in] A rank-2 tensor
 *
 * This function performs compound matrix operation \f$\mathbf {AA^T-A-A^T}\f$
 *-or- \f$A_{ij} A_{kj} -  A_{ik} -  A_{ki}\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::AijAkj_m_Aik_m_Aki( const R2TensorT< T_dim > & A )
{
  this->AijAkj( A );
  if( T_dim == 2 )
  {
    this->t_data[0] -= 2 * A.t_data[0];
    this->t_data[1] -= A.t_data[1] + A.t_data[2];
    this->t_data[2] -= 2 * A.t_data[3];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] -= 2 * A.t_data[0];

    this->t_data[1] -= A.t_data[1] + A.t_data[3];
    this->t_data[2] -= 2 * A.t_data[4];

    this->t_data[3] -= A.t_data[2] + A.t_data[6];
    this->t_data[4] -= A.t_data[5] + A.t_data[7];
    this->t_data[5] -= 2 * A.t_data[8];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::AijAkj_m_Aik_m_Aki not implemented for dimension > 3" );
    return;
  }
}


/**
 * @param[in] A symmetric rank-2 tensor
 * @param[in] Q rank-2 tensor
 *
 * This function performs compound matrix operation \f$\mathbf {QAQ^T}\f$ -or-
 *\f$Q_{ij} A_{jk} Q_{lk}\f$. This is
 * inteded to be a rotationAxis of a R2SymTensor with Q being an orthonormal
 * matrix.
 */
template< int T_dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void R2SymTensorT< T_dim >::QijAjkQlk( const R2SymTensorT< T_dim > & A, const R2TensorT< T_dim > & Q )
{
//  if (T_dim == 2)
//  {
//    this->t_data[0] = Q.t_data[0] * Q.t_data[0] * A.t_data[0] + 2 * Q.t_data[0] * Q.t_data[1] * A.t_data[1] +
// Q.t_data[1] * Q.t_data[1] * A.t_data[2];
//    this->t_data[1] = Q.t_data[0] * Q.t_data[3] * A.t_data[0] + Q.t_data[0] * Q.t_data[4] * A.t_data[1] + Q.t_data[1]
// * Q.t_data[3] * A.t_data[1] + Q.t_data[1]
//                      * Q.t_data[4] * A.t_data[2];
//    this->t_data[2] = Q.t_data[3] * Q.t_data[3] * A.t_data[0] + 2 * Q.t_data[3] * Q.t_data[4] * A.t_data[1] +
// Q.t_data[4] * Q.t_data[4] * A.t_data[2];
//  }
//  else if (T_dim == 3)
//  {
  realT o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13, o14, o15, o16, o17, o18, o19, o20, o21, o22, o23, o24;
  o1 = A.t_data[0] * Q.t_data[3];
  o2 = A.t_data[1] * Q.t_data[4];
  o3 = A.t_data[3] * Q.t_data[5];
  o4 = o1 + o2 + o3;
  o5 = A.t_data[1] * Q.t_data[3];
  o6 = A.t_data[2] * Q.t_data[4];
  o7 = A.t_data[4] * Q.t_data[5];
  o8 = o5 + o6 + o7;
  o9 = A.t_data[3] * Q.t_data[3];
  o10 = A.t_data[4] * Q.t_data[4];
  o11 = A.t_data[5] * Q.t_data[5];
  o12 = o10 + o11 + o9;
  o13 = A.t_data[0] * Q.t_data[6];
  o14 = A.t_data[1] * Q.t_data[7];
  o15 = A.t_data[3] * Q.t_data[8];
  o16 = o13 + o14 + o15;
  o17 = A.t_data[1] * Q.t_data[6];
  o18 = A.t_data[2] * Q.t_data[7];
  o19 = A.t_data[4] * Q.t_data[8];
  o20 = o17 + o18 + o19;
  o21 = A.t_data[3] * Q.t_data[6];
  o22 = A.t_data[4] * Q.t_data[7];
  o23 = A.t_data[5] * Q.t_data[8];
  o24 = o21 + o22 + o23;

  this->t_data[0] = Q.t_data[0] * (A.t_data[0] * Q.t_data[0] + A.t_data[1] * Q.t_data[1] + A.t_data[3] * Q.t_data[2])
                    + Q.t_data[1] * (A.t_data[1] * Q.t_data[0] + A.t_data[2] * Q.t_data[1] + A.t_data[4] * Q.t_data[2])
                    + Q.t_data[2] * (A.t_data[3] * Q.t_data[0] + A.t_data[4] * Q.t_data[1] + A.t_data[5] * Q.t_data[2]);

  this->t_data[1] = o4 * Q.t_data[0] + o8 * Q.t_data[1] + o12 * Q.t_data[2];
  this->t_data[2] = o4 * Q.t_data[3] + o8 * Q.t_data[4] + o12 * Q.t_data[5];
  this->t_data[3] = o16 * Q.t_data[0] + o20 * Q.t_data[1] + o24 * Q.t_data[2];
  this->t_data[4] = o16 * Q.t_data[3] + o20 * Q.t_data[4] + o24 * Q.t_data[5];
  this->t_data[5] = o16 * Q.t_data[6] + o20 * Q.t_data[7] + o24 * Q.t_data[8];
//  }
//  else
//  {
//    GEOSX_WARNING("R2SymTensorT::QijAjkQlk(R2TensorT) not implemented for dimension > 3");
//  }

}


/**
 * @param[in] a rank-1 tensor
 *
 * This function performs a dyadic product of a rank-1 tensor with itself
 *  \f$\mathbf {a \otimes a}\f$ -or- \f$a_i a_j\f$
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::dyadic_aa( const R1TensorT< T_dim > & a )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = a.t_data[0] * a.t_data[0];
    this->t_data[1] = a.t_data[0] * a.t_data[1];
    this->t_data[2] = a.t_data[1] * a.t_data[1];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = a.t_data[0] * a.t_data[0];
    this->t_data[1] = a.t_data[1] * a.t_data[0];
    this->t_data[2] = a.t_data[1] * a.t_data[1];
    this->t_data[3] = a.t_data[2] * a.t_data[0];
    this->t_data[4] = a.t_data[2] * a.t_data[1];
    this->t_data[5] = a.t_data[2] * a.t_data[2];
  }
  else
  {
    for( int i = 0; i < T_dim; ++i )
    {
      int index = 0;
      for( int ii = 1; ii <= i; ++ii )
        index += ii;
      for( int j = 0; j <= i; ++j, ++index )
        this->t_data[index] = a.t_data[i] * a.t_data[j];
    }
  }
}

/**
 * @param[in] a rank-1 tensor
 * @return none
 *
 * This function performs a dyadic product of a rank-1 tensor with itself
 *  \f$\mathbf {a \otimes a}\f$ -or- \f$a_i a_j\f$
 */
template< int T_dim >
R2SymTensorT< T_dim > DyadicProduct( const R1TensorT< T_dim > & a )
{
  R2SymTensorT< T_dim > c;
  c.dyadic_aa( a );
  return c;
}

/**
 * @param[in] a rank-1 tensor
 * @param[in] b rank-1 tensor
 *
 * This function adds the dyadic product of two rank-1 tensors to this
 * tensor.
 */
template< int T_dim >
inline void R2SymTensorT< T_dim >::dyadic_ab_plus_ba( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b )
{
  if( T_dim == 2 )
  {
    this->t_data[0] = a.t_data[0] * b.t_data[0] + b.t_data[0] * a.t_data[0];
    this->t_data[1] = a.t_data[0] * b.t_data[1] + b.t_data[0] * a.t_data[1];
    this->t_data[2] = a.t_data[1] * b.t_data[1] + b.t_data[1] * a.t_data[1];
  }
  else if( T_dim == 3 )
  {
    this->t_data[0] = a.t_data[0] * b.t_data[0] + b.t_data[0] * a.t_data[0];
    this->t_data[1] = a.t_data[1] * b.t_data[0] + b.t_data[1] * a.t_data[0];
    this->t_data[2] = a.t_data[1] * b.t_data[1] + b.t_data[1] * a.t_data[1];
    this->t_data[3] = a.t_data[2] * b.t_data[0] + b.t_data[2] * a.t_data[0];
    this->t_data[4] = a.t_data[2] * b.t_data[1] + b.t_data[2] * a.t_data[1];
    this->t_data[5] = a.t_data[2] * b.t_data[2] + b.t_data[2] * a.t_data[2];
  }
  else
  {
    GEOSX_WARNING( "R2SymTensorT::dyadic_ab(R1TensorT,R1TensorT) not implemented for dimension )> 3" );
  }
}

template< int T_dim >
inline void R2SymTensorT< T_dim >::AijklBkl( const R4minSymTensorT< T_dim > & A, const R2SymTensorT< T_dim > & B )
{
  int n_dim =  SymSize< T_dim >::value;

  for( int i=0; i<n_dim; ++i )
  {
    this->t_data[i] = 0.;
    for( int j=0; j<n_dim; ++j )
    {
      this->t_data[i] += A.t_data[i+j*n_dim]*B.t_data[j];
    }
  }
}

template< int T_dim >
inline void R2SymTensorT< T_dim >::AijklBij( const R4minSymTensorT< T_dim > & A, const R2SymTensorT< T_dim > & B )
{
  int n_dim =  SymSize< T_dim >::value;

  for( int i=0; i<n_dim; ++i )
  {
    this->t_data[i] = 0.;
    for( int j=0; j<n_dim; ++j )
    {
      this->t_data[i] += A.t_data[i*n_dim+j]*B.t_data[j];
    }
  }
}


#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
#ifdef  __IBMC__
#pragma report(enable, "1540-2907")
#endif

#endif
