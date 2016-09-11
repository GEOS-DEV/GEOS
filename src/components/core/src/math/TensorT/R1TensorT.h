//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security,
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief This file contains the definition of the R1TensorT class
 * @file R1TensorT.h
 * @author Randolph Settgast
 */

#ifndef R1_TENSOR_T_H_
#define R1_TENSOR_T_H_

#include "TensorBaseT.h"

#include <cstdlib>

template<int T_dim> class R2SymTensorT;
template<int T_dim> class R2TensorT;


/**
 * @brief R1TensorT is a rank-1 tensor object type
 * @author Randolph Settgast
 * @tparam T_dim length of tensor index
 *
 * R1TensorT derives from TensorBaseT, and defines basic operations that can be
 * done on a rank-1 tensor, as well as operations that result in a rank-1 tensor.
 */
template<int T_dim>
class R1TensorT : public TensorBaseT< T_dim >
{

//**** Overloaded arithmetic operators  ******************************************

// Scalar product
  friend R1TensorT<T_dim> operator*(realT k, const R1TensorT<T_dim> &V){ return V*k; }
  friend R1TensorT<T_dim> operator*(R1TensorT<T_dim> V, realT k){return V*=k; }
// Dot product *
  friend realT operator*(const R1TensorT<T_dim> &Va, const R1TensorT<T_dim> &Vb){return Dot(Va,Vb); }

// Division by scalar
  friend R1TensorT<T_dim> operator/(R1TensorT<T_dim> V, realT k){return V/=k; }


public:
  //**** CONSTRUCTORS AND DESTRUCTORS ******************************************

  /**
   * @author Randolph Settgast
   */
  R1TensorT( void ) : TensorBaseT< T_dim > () {}

  /**
   * @author Randolph Settgast
   * @param[in] data use for initialization of t_data
   */
  explicit R1TensorT( const realT data ) : TensorBaseT< T_dim >(data) {}

  /**
   * @author Randolph Settgast
   * @param[in] data naked array used for initialization of t_data
   */
  explicit R1TensorT( realT* const data ) : TensorBaseT< T_dim >(data) {}

  /**
   * @author walsh24
   * @param[in] data use for initialization of t_data
   */
  explicit R1TensorT( const int data ) : TensorBaseT< T_dim >( realT(data) ) {}

  //**** CONSTRUCTORS AND DESTRUCTORS *******************************************

  /**
   * @author Randolph Settgast
   * @param[in] rhs reference to R1TensorT object to use in initialization
   */
  R1TensorT( const R1TensorT< T_dim >& rhs ) : TensorBaseT< T_dim > ()
  { TensorBaseT< T_dim >::operator=( rhs ); }

  R1TensorT( const TensorBaseT< T_dim >& rhs ) : TensorBaseT< T_dim > ()
  { TensorBaseT< T_dim >::operator=( rhs ); }

  /**
   * @author Stuart Walsh
   *
   * Explicit constructors - will throw compile-time errors if not called with the correct dimension
   */
  R1TensorT(realT x,realT y);  //2D only
  R1TensorT(realT x,realT y, realT z); //3D only

  /// non-virtual destructor
  ~R1TensorT( void ) {}

  //***** ASSIGNMENT OPERATORS *************************************************
  /// assignment of all data to an integer
  R1TensorT< T_dim >& operator=( const int& rhs );

  /// assignment to all data to a realT
  R1TensorT< T_dim >& operator=( const realT& rhs );

  /// assignment to another R1TensorT
  R1TensorT< T_dim >& operator=( const R1TensorT< T_dim >& rhs );

  //***** ACCESS OPERATORS ****************************************************
  /// const access to data
  inline const realT& operator()( const int i ) const { return this->t_data[i];  }

  /// non-const access to data
  inline realT& operator()( const int i )             { return this->t_data[i];  }

  /// const access to data
  inline const realT& operator[]( const int i ) const { return this->t_data[i];  }

  /// non-const access to data
  inline realT& operator[]( const int i )             { return this->t_data[i];  }

  //***** MULTIPLICATION OPERATIONS *******************************************
  /// multiply (inner product) Rank2 tensor with Rank 1 tensor
  void AijBj( const R2TensorT< T_dim >& A, const R1TensorT< T_dim >& B );

  realT ProductOfSquares() const;

  /// subtract inner product of Rank2 tensor with Rank 1 tensor
  void minusAijBj( const R2TensorT< T_dim >& A, const R1TensorT< T_dim >& B );

  /// subtract inner product of Rank2 tensor with Rank 1 tensor
  void minusAijBj( const R2SymTensorT< T_dim >& A, const R1TensorT< T_dim >& B );

  /// multiply (inner product) transpose Rank2 tensor with Rank 1 tensor

  void AijBi( const R2TensorT< T_dim >& A, const R1TensorT< T_dim >& B );

  /// multiply (inner product) Symmetric Rank2 tensor with Rank 1 tensor
  void AijBj( const R2SymTensorT< T_dim >& A, const R1TensorT< T_dim >& B );

  /// Hadamard product between two Rank1 tensors
  void AiBi( const R1TensorT< T_dim >& A, const R1TensorT< T_dim >& B );

  /// permutation operator contracted on a Rank2 tensor
  void eijkAjk( const R2TensorT< T_dim >& A );

  /// cross product of 2 rank1 tensors
  void Cross( const R1TensorT< T_dim >& a, const R1TensorT< T_dim >& b );

  /// get a row from a symmetric rank2 tensor
  void GetRow( const R2SymTensorT< T_dim >& A, const int row );

  /// get a column from a symmetric rank2 tensor (same as GetRow())
  void GetCol( const R2SymTensorT< T_dim >& A, const int col );


  //****** TENSOR OPERATIONS **************************************************
  /// take the L2 norm of the tensor
  realT L2_Norm( void ) const;

  /// get the unit vector
  R1TensorT< T_dim > UnitVector( void ) const
  { realT n = this->L2_Norm(); return (n>0.0) ? (*this/n) : *this; }

  /// Normalize the vector
  realT Normalize( void )
  { realT n = this->L2_Norm(); if(n>0.0) *this /= n; return n; }

  /// sum the components of the tensor
  inline realT Sum( void ) const;

  //***** OUTPUT **************************************************************
  /// output
  void print( std::ostream& os ) const;

  //***** FRIEND DECLARATIONS *************************************************
  /// declare R2SymTensorT a friend so that it can access t_data directly
  friend class R2SymTensorT< T_dim >;

  /// declare R2TensorT a friend so that it can access t_data directly
  friend class R2TensorT< T_dim >;

  // define cross product
  friend inline
  R1TensorT< T_dim > Cross( const R1TensorT< T_dim >& a, const R1TensorT< T_dim >& b )
  {
    R1TensorT< T_dim > c;
    c.Cross(a,b);
    return c;
  }

private:
};


//*****************************************************************************
//***** END DECLARATION *******************************************************
//*****************************************************************************

/// Explicit 2D constructor
///
/// Template specialisation - if templated on another dimension constructor will throw a compile time error.
template<>
inline R1TensorT<2>::R1TensorT(realT x,realT y) :
  TensorBaseT< 2 >()
{
  this->t_data[0] = x;
  this->t_data[1] = y;
}


/// Explicit 3D constructor
///
/// Template specialisation - if templated on another dimension constructor will throw a compile time error.
template<>
inline R1TensorT<3>::R1TensorT(realT x,realT y,realT z) :
  TensorBaseT< 3 >()
{
  this->t_data[0] = x;
  this->t_data[1] = y;
  this->t_data[2] = z;
}



template<int T_dim>
void R1TensorT< T_dim >::print( std::ostream& os ) const
{
  for (int i = 0 ; i < T_dim ; ++i)
    os << (*this)( i ) << '\t';
}

#include "R2SymTensorT.h"
#include "R2TensorT.h"

//*****************************************************************************
//***** R1TensorT Member Function Definition **********************************
//*****************************************************************************



//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @author Randolph Settgast
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template<int T_dim>
inline R1TensorT< T_dim >& R1TensorT< T_dim >::operator=( const int& rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template<int T_dim>
inline R1TensorT< T_dim >& R1TensorT< T_dim >::operator=( const realT& rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

/**
 * @author Randolph Settgast
 * @param[in] rhs tensor to copy
 * @return reference to *this
 */
template<int T_dim>
inline R1TensorT< T_dim >& R1TensorT< T_dim >::operator=( const R1TensorT< T_dim >& rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

//***** TENSOR OPERATORS ******************************************************

/**
 * @author Randolph Settgast
 * @return L2 norm of tensor
 */
template<int T_dim>
inline realT R1TensorT< T_dim >::L2_Norm( void ) const
{
  realT norm = 0.0;

  for (int i = 1 ; i <= T_dim ; ++i)
    norm += this->t_data[i - 1] * this->t_data[i - 1];
  norm = sqrt( norm );

  return norm;
}


/**
 * @author Randolph Settgast
 * @return sum of the tensor components
 */
template<int T_dim>
inline realT R1TensorT< T_dim >::Sum( void ) const
{
  realT sum = 0.0;

  for (int i = 1 ; i <= T_dim ; ++i)
    sum += this->t_data[i - 1];

  return sum;
}

//***** MULTIPLICATION OPERATORS **********************************************
/**
 * @author Randolph Settgast
 * @param[in] A rank2 tensor
 * @param[in] B rank1 tensor
 * @return none
 *
 * this function contracts the input tensors and places the result into
 * this->t_data.
 */
template<int T_dim>
inline void R1TensorT< T_dim >::AijBj( const R2TensorT< T_dim >& A, const R1TensorT< T_dim >& B )
{
  if (T_dim == 1)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0];
  }
  if (T_dim == 2)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1];
    this->t_data[1] = A.t_data[2] * B.t_data[0] + A.t_data[3] * B.t_data[1];
  }
  else if (T_dim == 3)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1] + A.t_data[2] * B.t_data[2];
    this->t_data[1] = A.t_data[3] * B.t_data[0] + A.t_data[4] * B.t_data[1] + A.t_data[5] * B.t_data[2];
    this->t_data[2] = A.t_data[6] * B.t_data[0] + A.t_data[7] * B.t_data[1] + A.t_data[8] * B.t_data[2];
  }
  else
    std::cout << "R1TensorT::ProductOfSquares not implemented for nsdof>3";
}

/**
 * @author Scott Johnson
 * @param[in] A rank2 tensor
 * @param[in] B rank1 tensor
 * @return none
 *
 * Hadamard product
 */
template<int T_dim>
inline void R1TensorT< T_dim >::AiBi( const R1TensorT< T_dim >& A, const R1TensorT< T_dim >& B )
{
  for(int i = 0 ; i < T_dim ; i++)
    this->t_data[i] = A.t_data[i] * B.t_data[i];
}

/**
 * @author Scott Johnson
 * @return Product of squares
 */
template<int T_dim>
inline realT R1TensorT< T_dim >::ProductOfSquares() const
{
  if (T_dim == 1)
  {
    return this->t_data[0]*this->t_data[0];
  }
  if (T_dim == 2)
  {
    return (this->t_data[0]*this->t_data[0]) *
           (this->t_data[1]*this->t_data[1]);
  }
  else if (T_dim == 3)
  {
    return (this->t_data[0]*this->t_data[0]) *
           (this->t_data[1]*this->t_data[1]) *
           (this->t_data[2]*this->t_data[2]);
  }
  else
  {
    std::cout << "R1TensorT::ProductOfSquares not implemented for nsdof>3";
    exit(1);
  }

  return 0.0; // never get here
}

/**
 * @author Randolph Settgast
 * @param[in] A rank2 tensor
 * @param[in] B rank1 tensor
 * @return none
 *
 * this function contracts the input tensors and subtracts the result into
 * this->t_data.
 */
template<int T_dim>
inline void R1TensorT< T_dim >::minusAijBj( const R2TensorT< T_dim >& A, const R1TensorT< T_dim >& B )
{
  if (T_dim == 1)
  {
    this->t_data[0] -= A.t_data[0] * B.t_data[0];
  }
  if (T_dim == 2)
  {
    this->t_data[0] -= A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1];
    this->t_data[1] -= A.t_data[2] * B.t_data[0] + A.t_data[3] * B.t_data[1];
  }
  else if (T_dim == 3)
  {
    this->t_data[0] -= A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1] + A.t_data[2] * B.t_data[2];
    this->t_data[1] -= A.t_data[3] * B.t_data[0] + A.t_data[4] * B.t_data[1] + A.t_data[5] * B.t_data[2];
    this->t_data[2] -= A.t_data[6] * B.t_data[0] + A.t_data[7] * B.t_data[1] + A.t_data[8] * B.t_data[2];
  }
  else
    std::cout << "R1TensorT::ProductOfSquares not implemented for nsdof>3";

}


template<int T_dim>
inline void R1TensorT< T_dim >::minusAijBj( const R2SymTensorT< T_dim >& A, const R1TensorT< T_dim >& B )
{
  if (T_dim == 1)
  {
    this->t_data[0] -= A.t_data[0] * B.t_data[0];
  }
  if (T_dim == 2)
  {
    this->t_data[0] -= A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1];
    this->t_data[1] -= A.t_data[1] * B.t_data[0] + A.t_data[2] * B.t_data[1];
  }
  else if (T_dim == 3)
  {
    this->t_data[0] -= A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1] + A.t_data[3] * B.t_data[2];
    this->t_data[1] -= A.t_data[1] * B.t_data[0] + A.t_data[2] * B.t_data[1] + A.t_data[4] * B.t_data[2];
    this->t_data[2] -= A.t_data[3] * B.t_data[0] + A.t_data[4] * B.t_data[1] + A.t_data[5] * B.t_data[2];
  }
  else
    std::cout << "R1TensorT::ProductOfSquares not implemented for nsdof>3";

}
/**
 * @author Randolph Settgast
 * @param[in] A rank2 tensor
 * @param[in] B rank1 tensor
 * @return none
 *
 * this function contracts the input tensors and places the result into
 * this->t_data.
 */
template<int T_dim>
inline void R1TensorT< T_dim >::AijBi( const R2TensorT< T_dim >& A, const R1TensorT< T_dim >& B )
{
  if (T_dim == 2)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[2] * B.t_data[1];
    this->t_data[1] = A.t_data[1] * B.t_data[0] + A.t_data[3] * B.t_data[1];
  }
  else if (T_dim == 3)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[3] * B.t_data[1] + A.t_data[6] * B.t_data[2];
    this->t_data[1] = A.t_data[1] * B.t_data[0] + A.t_data[4] * B.t_data[1] + A.t_data[7] * B.t_data[2];
    this->t_data[2] = A.t_data[2] * B.t_data[0] + A.t_data[5] * B.t_data[1] + A.t_data[8] * B.t_data[2];
  }
  else
    std::cout << "R1TensorT not implemented for nsdof>3";
}

/**
 * @author Randolph Settgast
 * @param[in] A symmetric rank2 tensor
 * @param[in] B rank1 tensor
 * @return none
 *
 * this function contracts the input tensors and places the result into
 * this->t_data.
 */
template<int T_dim>
inline void R1TensorT< T_dim >::AijBj( const R2SymTensorT< T_dim >& A, const R1TensorT< T_dim >& B )
{
#ifdef __INTEL_COMPILER
#pragma warning push
#pragma warning disable 175
#endif
#ifdef  __IBMC__
#pragma report(disable, "1540-2907")
#endif
  if (T_dim == 2)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1];
    this->t_data[1] = A.t_data[1] * B.t_data[0] + A.t_data[2] * B.t_data[1];
  }
  else if (T_dim == 3)
  {
    this->t_data[0] = A.t_data[0] * B.t_data[0] + A.t_data[1] * B.t_data[1] + A.t_data[3] * B.t_data[2];
    this->t_data[1] = A.t_data[1] * B.t_data[0] + A.t_data[2] * B.t_data[1] + A.t_data[4] * B.t_data[2];
    this->t_data[2] = A.t_data[3] * B.t_data[0] + A.t_data[4] * B.t_data[1] + A.t_data[5] * B.t_data[2];
  }
  else
    std::cout << "R1TensorT not implemented for nsdof>3";
}

/**
 * @author Randolph Settgast
 * @param[in] A rank2 tensor
 * @return none
 *
 * this function contracts the permutation operator on two indicies of a
 * rank2 tensor and places the result in this->tdata
 */
template<int T_dim>
inline void R1TensorT< T_dim >::eijkAjk( const R2TensorT< T_dim >& A )
{
  if (T_dim == 3)
  {
    this->t_data[0] = A.t_data[5] - A.t_data[7];
    this->t_data[1] = A.t_data[6] - A.t_data[2];
    this->t_data[2] = A.t_data[1] - A.t_data[3];
  }
  else
    std::cout << "R1TensorT not implemented for nsdof>3";

}

/**
 * @author Randolph Settgast
 * @param[in] a rank1 tensor
 * @param[in] b rank1 tensor
 * @return none
 *
 * this function takes the cross product of two rank1 tensors and places the
 * result into this->tdata
 */
template<int T_dim>
inline void R1TensorT< T_dim >::Cross( const R1TensorT< T_dim >& a, const R1TensorT< T_dim >& b )
{
  if (T_dim == 3)
  {
    this->t_data[0] = a.t_data[1] * b.t_data[2] - a.t_data[2] * b.t_data[1];
    this->t_data[1] = -(a.t_data[0] * b.t_data[2] - a.t_data[2] * b.t_data[0]);
    this->t_data[2] = a.t_data[0] * b.t_data[1] - a.t_data[1] * b.t_data[0];
  }
  else
    std::cout << "R1TensorT not implemented for nsdof>3";

}

/**
 * @author Randolph Settgast
 * @param[in] A symmetric rank2 tensor
 * @param[in] row row number to extract
 * @return none
 *
 * this function extracts a row from A and places it in this->tdata
 */
template<int T_dim>
inline void R1TensorT< T_dim >::GetRow( const R2SymTensorT< T_dim >& A, const int row )
{

  if (T_dim == 3)
  {
    if (row == 1)
    {
      this->t_data[0] = A.t_data[0];
      this->t_data[1] = A.t_data[1];
      this->t_data[2] = A.t_data[3];
    }
    else if (row == 2)
    {
      this->t_data[0] = A.t_data[1];
      this->t_data[1] = A.t_data[2];
      this->t_data[2] = A.t_data[4];
    }
    else if (row == 3)
    {
      this->t_data[0] = A.t_data[3];
      this->t_data[1] = A.t_data[4];
      this->t_data[2] = A.t_data[5];
    }
  }
  else
    std::cout << "R1TensorT not implemented for nsdof>3";
}

/**
 * @author Randolph Settgast
 * @param[in] A symmetric rank2 tensor
 * @param[in] col col number to extract
 * @return none
 *
 * this function extracts a column from A and places it in this->tdata
 */
template<int T_dim>
inline void R1TensorT< T_dim >::GetCol( const R2SymTensorT< T_dim >& A, const int col )
{

  if (T_dim == 3)
  {
    if (col == 1)
    {
      this->t_data[0] = A.t_data[0];
      this->t_data[1] = A.t_data[1];
      this->t_data[2] = A.t_data[3];
    }
    else if (col == 2)
    {
      this->t_data[0] = A.t_data[1];
      this->t_data[1] = A.t_data[2];
      this->t_data[2] = A.t_data[4];
    }
    else if (col == 3)
    {
      this->t_data[0] = A.t_data[3];
      this->t_data[1] = A.t_data[4];
      this->t_data[2] = A.t_data[5];
    }
  }
  else
    std::cout << "R1TensorT not implemented for nsdof>3";
}

#endif
