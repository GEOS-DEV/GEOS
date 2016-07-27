/*
 * R6minSymTensorT.h
 *
 *  Created on: Aug 26, 2013
 *
 */

#ifndef _R6MINTENSORT_H_
#define _R6MINTENSORT_H_

#include "TensorBaseT.h"
//#include "R4minSymTensorT.h"
//#include "R2SymTensorT.h"
//#include "R2TensorT.h"


template<int T_dim> class R2SymTensorT;
template<int T_dim> class R2TensorT;
template<int T_dim> class R4minSymTensorT;
template<int T_dim> struct SymSize;

//*****************************************************************************
//*****  R6minTensorT Declaration  **********************************************
//*****************************************************************************
template<int T_dim>
class R6minSymTensorT :
public TensorBaseT<SymSize< T_dim >::value*SymSize< T_dim >::value*SymSize< T_dim >::value>
{
  public:
    //**** CONSTRUCTORS AND DESTRUCTORS *****************************************
    /// default constructor
    R6minSymTensorT(void);

    // destructor
    ~R6minSymTensorT(void);


    /// copy constructor
    //R6minSymTensorT(const R6minSymTensorT<T_dim>& rhs);

    //***** ACCESS OPERATORS ****************************************************
    /// const access to data
    const realT& operator()( const int i , const int j , const int k , const int l ,
                             const int m , const int n ) const;

    /// non-const access to data
    realT& operator()( const int i , const int j , const int k , const int l ,
                       const int m , const int n ) ;


    friend class R2TensorT<T_dim> ;
    friend class R2SymTensorT<T_dim> ;
    friend class R4minSymTensorT<T_dim> ;

  private:

};
//*****************************************************************************
//***** END DECLARATIONS ******************************************************
//*****************************************************************************

#include "R2SymTensorT.h"
#include "R2TensorT.h"
//#include "R4minSymTensorT.h"

template< int T_dim >
R6minSymTensorT<T_dim>::R6minSymTensorT(void):
TensorBaseT<SymSize< T_dim >::value*SymSize< T_dim >::value*SymSize< T_dim >::value>()
{
}

template< int T_dim >
R6minSymTensorT<T_dim>::~R6minSymTensorT(void)
{
}

    //***** ACCESS OPERATORS ******************************************************


template< int T_dim >
inline const realT& R6minSymTensorT<T_dim>::operator()( const int i , const int j ,
                                                        const int k , const int l ,
                                                        const int m , const int n ) const
{
  int n_dim = (SymSize< T_dim >::value);
  int i_sym = i, k_sym = k, m_sym = m;
  int j_sym = j, l_sym = l, n_sym = n;
  int index_0 = 0, index_1 = 0, index_2 = 0;

  if (j > i)
  {
    i_sym = j;
    j_sym = i;
  }

  if (l > k)
  {
    k_sym = l;
    l_sym = k;
  }

  if (n > m)
  {
    m_sym = n;
    n_sym = m;
  }

  for (int ii = 1; ii <= i_sym; ++ii)
    {index_0 += ii;}
  index_0 += j_sym;

  for (int jj = 1; jj <= k_sym; ++jj)
    {index_1 += jj;}
  index_1 += l_sym;

  for (int kk = 1; kk <= m_sym; ++kk)
    {index_2 += kk;}
  index_2 += n_sym;

  return this->t_data[index_0+index_1*n_dim+index_2*n_dim*n_dim];
}
//
//
template< int T_dim >
inline realT& R6minSymTensorT<T_dim>::operator()( const int i , const int j ,
                                                  const int k , const int l ,
                                                  const int m , const int n )
{
  int n_dim = (SymSize< T_dim >::value);
  int i_sym = i, k_sym = k, m_sym = m;
  int j_sym = j, l_sym = l, n_sym = n;

  if (j > i)
  {
    i_sym = j;
    j_sym = i;
  }

  if (l > k)
  {
    k_sym = l;
    l_sym = k;
  }

  if (n > m)
  {
    m_sym = n;
    n_sym = m;
  }

  int index_0 = 0, index_1 = 0, index_2 = 0;

  for(int ii = 1; ii <= i_sym; ++ii)
    {index_0 += ii;}
  index_0 += j_sym;

  for(int jj = 1; jj <= k_sym; ++jj)
    {index_1 += jj;}
  index_1 += l_sym;

  for(int kk = 1; kk <= m_sym; ++kk)
    {index_2 += kk;}
  index_2 += n_sym;

  return this->t_data[index_0+index_1*n_dim+index_2*n_dim*n_dim];
}


#endif
