/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PythonFunction.cpp
 */

#include "PythonFunction.hpp"

namespace geos
{

template< typename INDEX_T >
PythonFunction<INDEX_T>::PythonFunction(const string & name,
                                        Group * const parent ): 
                                  dataRepository::Group( name, parent ),
                                  py_evaluate_func(nullptr), 
                                  py_derivative_func(nullptr)
{
  // Initialize Python if needed (ensure Python interpreter is initialized)
  if (!Py_IsInitialized()) 
    Py_Initialize();
}

template< typename INDEX_T >
PythonFunction<INDEX_T>::~PythonFunction()
{
  // Decrement reference count for the evaluation function if set
  if (py_evaluate_func)
    Py_XDECREF(py_evaluate_func);

  // Decrement reference count for the derivative function if set (optional)
  if (py_derivative_func)
    Py_XDECREF(py_derivative_func);
}

template< typename INDEX_T >
void
PythonFunction<INDEX_T>::setEvaluateFunction(PyObject* pyFunc, PyObject* pyDerivativeFunc)
{
  // Set the evaluation function (required)
  if (PyCallable_Check(pyFunc)) 
  {
    Py_XINCREF(pyFunc);  // Increment reference count for the new function
    Py_XDECREF(py_evaluate_func);  // Decrement reference count of the old function
    py_evaluate_func = pyFunc;  // Assign the new function
  } 
  else 
  {
    GEOS_THROW( GEOS_FMT( "{}: Provided Python function is not callable.", 
                        getName()), InputError );
  }

  // Set the derivative function (optional)
  if (pyDerivativeFunc) 
  {
    if (PyCallable_Check(pyDerivativeFunc)) 
    {
      Py_XINCREF(pyDerivativeFunc);  // Increment reference count for the new derivative function
      Py_XDECREF(py_derivative_func); // Decrement reference count of the old derivative function
      py_derivative_func = pyDerivativeFunc;  // Assign the new derivative function
    } 
    else 
    {
      GEOS_THROW( GEOS_FMT( "{}: Provided Python derivative function is not callable.", 
                          getName()), InputError );
    }
  }
}

template< typename INDEX_T >
void
PythonFunction<INDEX_T>::setDimensions( integer numberOfDimesions,
                                        integer numberOfOperators )
{
  numDims = numberOfDimesions;
  numOps = numberOfOperators;
  numVerts = 1 << numDims;
}

template< typename INDEX_T >
void
PythonFunction<INDEX_T>::setAxes( integer numberOfDimesions,
                                  integer numberOfOperators,
                                  real64_array axisMinimums,
                                  real64_array axisMaximums,
                                  integer_array axisPoints )
{
  numDims = numberOfDimesions;
  numOps = numberOfOperators;
  numVerts = 1 << numDims;

  m_axisMinimums = axisMinimums;
  m_axisMaximums = axisMaximums;
  m_axisPoints = axisPoints;

  m_axisSteps.resize(numDims);
  m_axisStepInvs.resize(numDims);
  m_axisHypercubeMults.resize(numDims);

  for( integer dim = 0; dim < numDims; dim++ )
  {
    m_axisSteps[dim] = (m_axisMaximums[dim] - m_axisMinimums[dim]) / (m_axisPoints[dim] - 1);
    m_axisStepInvs[dim] = 1 / m_axisSteps[dim];
  }

  m_axisHypercubeMults[numDims - 1] = 1;
  for( integer dim = numDims - 2; dim >= 0; --dim )
  {
    m_axisHypercubeMults[dim] = m_axisHypercubeMults[dim + 1] * (m_axisPoints[dim + 1] - 1);
  }
}


template class PythonFunction<__uint128_t>;
// REGISTER_CATALOG_ENTRY( dataRepository::Group, PythonFunction<__uint128_t>, string const &, dataRepository::Group * const )

} // end of namespace geos

