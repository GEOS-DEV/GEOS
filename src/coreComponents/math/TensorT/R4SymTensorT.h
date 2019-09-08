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

/******************************************************************************
 * R2SymTensorT.h - Rank 2 Symmetric Tensor Class
 *
 * created: RRS (01/18/2002)
 *****************************************************************************/

#ifndef _R4_SYM_TENSOR_T_H_
#define _R4_SYM_TENSOR_T_H_

#include "TensorBaseT.h"

template<int N>
struct SymSize
{
  enum
  {
    value = N + SymSize<N - 1>::value
  };
};

template<>
struct SymSize<1>
{
  enum
  {
    value = 1
  };
};

template<int T_dim>
class R2SymTensorT;

//*****************************************************************************
//***** R4SymTensorT Declaration **********************************************
//*****************************************************************************
template<int T_dim>
class R4SymTensorT : public TensorBaseT<SymSize<T_dim>::value>
{
public:
  //**** CONSTRUCTORS AND DESTRUCTORS *****************************************
  R4SymTensorT(void);
  ~R4SymTensorT(void);
  R4SymTensorT(const R4SymTensorT<T_dim>& rhs);

  //***** ASSIGNMENT OPERATORS
  // **************************************************
  R4SymTensorT<T_dim>&
  operator=(const int& rhs);
  R4SymTensorT<T_dim>&
  operator=(const realT& rhs);
#if __LONG_real
  R4SymTensorT<T_dim>& operator=( const realT& rhs )
  {
    operator=(static_cast<realT>(rhs));
    return *this;
  }
#endif

  R4SymTensorT<T_dim>&
  operator=(const R4SymTensorT<T_dim>& rhs);

  //***** ACCESS OPERATORS ****************************************************
  inline realT
  operator()(const int i, const int j) const;
  inline realT&
  operator()(const int i, const int j);

  //***** MULTIPLICATION OPERATIONS *******************************************


  //****** TENSOR OPERATIONS **************************************************


  friend class R2TensorT<T_dim>;

private:
  R4SymTensorT(R4SymTensorT<T_dim>&);

};


template<int T_dim>
R4SymTensorT<T_dim>::R4SymTensorT(void):
  TensorBaseT<SymSize<T_dim>::value> ()
{}

template<int T_dim>
R4SymTensorT<T_dim>::~R4SymTensorT(void)
{}

template<int T_dim>
R4SymTensorT<T_dim>::R4SymTensorT(const R4SymTensorT<T_dim>& rhs):
  TensorBaseT<SymSize<T_dim>::value> ()
{
  TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
}

//***** ASSIGNMENT OPERATORS **************************************************

// Assigns all components to an integer
template<int T_dim>
R4SymTensorT<T_dim>&
R4SymTensorT<T_dim>::operator=(const int& rhs)
{
  TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
  return *this;
}

// Assigns all components to a realT
template<int T_dim>
R4SymTensorT<T_dim>&
R4SymTensorT<T_dim>::operator=(const realT& rhs)
{
  TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
  return *this;
}

// Assigns all components to another TensorBaseT's (Copy Constructor)
template<int T_dim>
R4SymTensorT<T_dim>&
R4SymTensorT<T_dim>::operator=(const R4SymTensorT<T_dim>& rhs)
{
  TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
  return *this;
}

template<int T_dim>
realT
R4SymTensorT<T_dim>::operator()(const int i, const int j) const
{
  int index = 0;
  int i_sym = i;
  int j_sym = j;

  if (j > i)
  {
    i_sym = j;
    j_sym = i;
  }

  for (int k = 1 ; k < i_sym ; ++k)
    index += k;
  index += (j_sym - 1);

  return this->t_data[index];
}

template<int T_dim>
inline realT&
R4SymTensorT<T_dim>::operator()(const int i, const int j)
{
  int index = 0;
  int i_sym = i;
  int j_sym = j;

  if (j > i)
  {
    i_sym = j;
    j_sym = i;
  }

  for (int k = 1 ; k < i_sym ; ++k)
    index += k;
  index += (j_sym - 1);

  return this->t_data[index];
}

#endif
