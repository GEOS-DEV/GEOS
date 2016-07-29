/*
 * R4minSymTensorT.h
 *
 *  Created on: Aug 26, 2013
 *      Author: xu11
 */

#ifndef _R4MINSYMTENSORT_H_
#define _R4MINSYMTENSORT_H_

#include "TensorBaseT.h"
//#include "R6minSymTensorT.h"
//#include "R2SymTensorT.h"
#include "R2TensorT.h"



template<int T_dim> class R2SymTensorT;
template<int T_dim> class R2TensorT;
template<int T_dim> class R6minSymTensorT;
template<int T_dim> struct SymSize;

//*****************************************************************************
//*****  R4minSymTensorT Declaration  **********************************************
//*****************************************************************************
template<int T_dim>
class R4minSymTensorT :
public TensorBaseT <SymSize< T_dim >::value*SymSize< T_dim >::value>
{
  public:
    //**** CONSTRUCTORS AND DESTRUCTORS *****************************************
    /// default constructor
    R4minSymTensorT(void);

    // destructor
    ~R4minSymTensorT(void);

    /// copy constructor
    R4minSymTensorT(const R4minSymTensorT<T_dim>& rhs);

    //***** ASSIGNMENT OPERATORS **************************************************
    //R4minSymTensorT<T_dim>& operator=(const int& rhs);
    //R4minSymTensorT<T_dim>& operator=(const realT& rhs);


    R4minSymTensorT<T_dim>& operator=(const R4minSymTensorT<T_dim>& rhs);

    //***** ACCESS OPERATORS ****************************************************
    /// const access to data
    const realT&
    operator()(const int i, const int j, const int k, const int l) const;

    /// non-const access to data
    realT&
    operator()(const int i, const int j, const int k, const int l);

    //***** MULTIPLICATION OPERATIONS *******************************************
    void AijmnBmnkl( const R4minSymTensorT<T_dim>& A , const R4minSymTensorT<T_dim>& B );
    void AijklmnBmn( const R6minSymTensorT<T_dim>& A , const R2SymTensorT<T_dim>& B );
    void AijBkl( const R2SymTensorT<T_dim>& A , const R2SymTensorT<T_dim>& B );

    //****** TENSOR TRSNSFORMATION **********************************************

    void Aij_to_Bklmn( const R2TensorT<6>& A);

    //****** INVERSE OF TENSOR **********************************************

    void Inverse_4(void) { return Inverse_4(*this); }

    void Inverse_4( R4minSymTensorT<T_dim>& A );

    friend class R2TensorT<T_dim> ;
    friend class R2TensorT<6>;
    friend class R2SymTensorT<T_dim> ;
    friend class R6minSymTensorT<T_dim> ;

  private:
  //  R4minSymTensorT(R4minSymTensorT<T_dim>&);

};
//*****************************************************************************
//***** END DECLARATIONS ******************************************************
//*****************************************************************************

#include "R2SymTensorT.h"
#include "R2TensorT.h"
#include "R6minSymTensorT.h"


template<int T_dim>
R4minSymTensorT<T_dim>::R4minSymTensorT(void):
TensorBaseT<SymSize< T_dim >::value*SymSize< T_dim >::value> ()
{
}

template<int T_dim>
R4minSymTensorT<T_dim>::~R4minSymTensorT(void)
{
}

template<int T_dim>
R4minSymTensorT<T_dim>::R4minSymTensorT(const R4minSymTensorT<T_dim>& rhs) :
TensorBaseT<SymSize< T_dim >::value*SymSize< T_dim >::value> ()
{
  TensorBaseT<SymSize< T_dim >::value*SymSize< T_dim >::value>::operator=(rhs);
}

//***** ASSIGNMENT OPERATORS **************************************************

// Assigns all components to an integer
/*
template<int T_dim>
  R4minSymTensorT<T_dim>&
  R4minSymTensorT<T_dim>::operator=(const int& rhs)
  {
    TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
    return *this;
  }

// Assigns all components to a realT
template<int T_dim>
  R4minSymTensorT<T_dim>&
  R4minSymTensorT<T_dim>::operator=(const realT& rhs)
  {
    TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
    return *this;
  }
*/

// Assigns all components to another TensorBaseT's (Copy Constructor)
template<int T_dim>
R4minSymTensorT<T_dim>&
R4minSymTensorT<T_dim>::operator=(const R4minSymTensorT<T_dim>& rhs)
{
  TensorBaseT<SymSize< T_dim >::value*SymSize< T_dim >::value>::operator=(rhs);
  return *this;
}

template<int T_dim>
inline const realT&
R4minSymTensorT<T_dim>::operator()(const int i, const int j, const int k, const int l) const
{
  int n_dim = (SymSize< T_dim >::value);
  int i_sym = i, k_sym = k;
  int j_sym = j, l_sym = l;

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

  int index0 = 0, index1 = 0;

  for(int ii = 1; ii <= i_sym; ++ii)
    index0 += ii;
  index0 += j_sym;

  for(int ii = 1; ii <= k_sym; ++ii)
    index1 += ii;
  index1 += l_sym;

  return this->t_data[index0+index1*n_dim];
}

template<int T_dim>
inline realT&
R4minSymTensorT<T_dim>::operator()(const int i, const int j, const int k, const int l)
{
  int n_dim = (SymSize< T_dim >::value);
  int i_sym = i, k_sym = k;
  int j_sym = j, l_sym = l;

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

  int index0 = 0, index1 = 0;


  for(int ii = 1; ii <= i_sym; ++ii)
    index0 += ii;
  index0 += j_sym;

  for(int ii = 1; ii <= k_sym; ++ii)
    index1 += ii;
  index1 += l_sym;

  return this->t_data[index0+index1*n_dim];
}

template<int T_dim>
inline void R4minSymTensorT< T_dim >::AijmnBmnkl(const R4minSymTensorT< T_dim >& A,
                                                 const R4minSymTensorT< T_dim >& B)
{
  int n_dim = (SymSize< T_dim >::value);

  for(int i=0; i<n_dim*n_dim;++i)
  {
    int x=(i % n_dim);
    int y=i-x;

    for(int j=0; j<n_dim; ++j)
    {
      this->t_data[i]+=A.t_data[j*n_dim+x]*B.t_data[j+y]*2;
    }

    for(int k=0,ii=0;k<n_dim;k+=(++ii)+1 )
    {
      this->t_data[i]-=A.t_data[k*n_dim+x]*B.t_data[k+y];
    }
  }

}

template<int T_dim>
inline void R4minSymTensorT< T_dim >::AijklmnBmn( const R6minSymTensorT<T_dim>& A ,
                                                  const R2SymTensorT<T_dim>& B )
{
  int n_dim = (SymSize< T_dim >::value);

  for(int i=0; i<n_dim*n_dim;++i)
  {



    for(int j=0; j<n_dim; ++j)
    {
      this->t_data[i]+=A.t_data[j*n_dim*n_dim+i]*B.t_data[j]*2;
    }

    for(int k=0,ii=0;k<n_dim;k+=(++ii)+1 )
    {
      this->t_data[i]-=A.t_data[k*n_dim*n_dim+i]*B.t_data[k];
    }
  }
}

template<int T_dim>
inline void R4minSymTensorT< T_dim >::AijBkl( const R2SymTensorT<T_dim>& A ,
                                              const R2SymTensorT<T_dim>& B )
{
  int n_dim = (SymSize< T_dim >::value);

  for(int i=0; i<n_dim;++i)
  {

    for(int j=0; j<n_dim; ++j)
    {
      this->t_data[i+j*n_dim]=A.t_data[i]*B.t_data[j];
    }
  }
}

template<int T_dim>
void R4minSymTensorT<T_dim>::Aij_to_Bklmn( const R2TensorT<6>& A)
{
  int n_dim = 6;
  if (T_dim!=3)
    std::cout<< "R4minSymTensorT<T_dim>::Aij_to_Bklmn not implemented for T_dim /= 3";
  else
  {

    for(int ii=0,m=0,j=0; m<n_dim; m+=++ii)
    {
      for(int i=0, c=m; c<n_dim;c+=(++i)+1+ii, ++j)
      {
        this->t_data[0+c] = A.t_data[j*n_dim];
        this->t_data[12+c] = A.t_data[j*n_dim+1];
        this->t_data[30+c] = A.t_data[j*n_dim+2];
        this->t_data[6+c] = A.t_data[j*n_dim+3];
        this->t_data[24+c] = A.t_data[j*n_dim+4];
        this->t_data[18+c] = A.t_data[j*n_dim+5];
      }
    }
  }
}

template<int T_dim>
void R4minSymTensorT<T_dim>::Inverse_4( R4minSymTensorT<T_dim>& A )
{
  if (T_dim!=3)
    std::cout<< "R4minSymTensorT<T_dim>::Inverse_4 not implemented for T_dim /= 3";
  else
  {
  R2TensorT<6> B;
  B.Aijkl_to_Bmn(A);
  B.Inverse();
  this->Aij_to_Bklmn(B);
  }
}

#endif /* R4MINSYMTENSORT_H_ */

