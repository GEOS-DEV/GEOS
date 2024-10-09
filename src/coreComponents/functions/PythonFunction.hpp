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
 * @file PythonFunction.hpp
 */

#ifndef GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_
#define GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_

#include <Python.h>
#include "common/DataTypes.hpp"

namespace geos
{

template< integer NUM_DIMS, integer NUM_OPS >
class PythonFunction
{
public:
  /// Compile time value for the number of table dimensions (inputs)
  static constexpr integer numDims = NUM_DIMS;

  /// Compile time value for the number of operators (interpolated functions, outputs)
  static constexpr integer numOps = NUM_OPS;

  /// Python function reference to be called in evaluate
  mutable PyObject* py_evaluate_func = nullptr;

  /// Python function reference to be called for evaluation of derivatives
  mutable PyObject* py_derivative_func = nullptr;

  /**
   * @brief Construct a new wrapper for python function
   */
  PythonFunction(): py_evaluate_func(nullptr), 
                    py_derivative_func(nullptr)
  {
    // Initialize Python if needed (ensure Python interpreter is initialized)
    if (!Py_IsInitialized()) 
      Py_Initialize();
  }

  /**
   * @brief Desctructor of a wrapper of python function
   */
  ~PythonFunction() 
  {
    // Decrement reference count for the evaluation function if set
    if (py_evaluate_func)
      Py_XDECREF(py_evaluate_func);

    // Decrement reference count for the derivative function if set (optional)
    if (py_derivative_func)
      Py_XDECREF(py_derivative_func);
  }

  /**
   * @brief Set Python functions for values and (optionally) for their derivatives
   * 
   * @param[in] pyFunc pointer to the function
   * @param[in] pyDerivativeFunc pointer to the function evaluating derivatives
   */
  void setEvaluateFunction(PyObject* pyFunc, PyObject* pyDerivativeFunc = nullptr) 
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
      throw std::runtime_error("Provided Python evaluation function is not callable.");
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
        throw std::runtime_error("Provided Python derivative function is not callable.");
      }
    }
  }

  /**
   * @brief interpolate all operators values at a given point
   * 
   * @tparam IN_ARRAY type of input array of coordinates
   * @tparam OUT_ARRAY type of output array of values
   * @param[in] state function arguments
   * @param[out] values function values
   */
  template< typename IN_ARRAY, typename OUT_ARRAY >
  GEOS_HOST_DEVICE   
  inline  
  void 
  evaluate( IN_ARRAY const & state, 
            OUT_ARRAY && values ) const
  {
    if (py_evaluate_func && PyCallable_Check(py_evaluate_func))
    {
      // Prepare Python arguments (state and values as lists)
      PyObject* py_state = PyList_New(numDims);
      PyObject* py_values = PyList_New(numOps);

      for (integer i = 0; i < numDims; ++i)
        PyList_SetItem(py_state, i, PyFloat_FromDouble(state[i]));

      for (integer i = 0; i < numOps; ++i)
        PyList_SetItem(py_values, i, PyFloat_FromDouble(values[i]));

      // Pack arguments into a tuple
      PyObject* args = PyTuple_New(2);
      PyTuple_SetItem(args, 0, py_state);
      PyTuple_SetItem(args, 1, py_values);

      // Call the Python function
      PyObject* result = PyObject_CallObject(py_evaluate_func, args);

      if (result == nullptr) 
      {
        PyErr_Print();  // Print Python error if any
        throw std::runtime_error("Python function call failed.");
      } 
      else 
      {
        // Update values with results from Python function
        for (integer i = 0; i < numOps; ++i)
          values[i] = PyFloat_AsDouble(PyList_GetItem(py_values, i));

        Py_DECREF(result);  // Decrement reference count of result
      }

      // Clean up references
      Py_DECREF(args);
    }
    else
    {
      throw std::runtime_error("No valid Python function set or the function is not callable.");
    }
  }

  /**
   * @brief interpolate all operators values at a given point
   * 
   * @tparam IN_ARRAY type of input array of coordinates
   * @tparam OUT_ARRAY type of output array of values
   * @tparam OUT_2D_ARRAY type of output array of derivatives
   * @param[in] state function arguments
   * @param[out] values function values
   * @param[out] derivatives derivatives of function values w.r.t. arguments
   */
  template< typename IN_ARRAY, typename OUT_ARRAY, typename OUT_2D_ARRAY >
  GEOS_HOST_DEVICE   
  inline  
  void 
  evaluate(IN_ARRAY const & state,
           OUT_ARRAY && values, 
           OUT_2D_ARRAY && derivatives) const
  {
    // Values
    if (py_evaluate_func && PyCallable_Check(py_evaluate_func))
    {
      // Prepare Python arguments (state, values, and derivatives as lists)
      PyObject* py_state = PyList_New(numDims);
      PyObject* py_values = PyList_New(numOps);

      // Populate Python lists with initial values for state 
      for (integer i = 0; i < numDims; ++i)
        PyList_SetItem(py_state, i, PyFloat_FromDouble(state[i]));

      for (integer i = 0; i < numOps; ++i)
        PyList_SetItem(py_values, i, PyFloat_FromDouble(values[i]));

      // Pack arguments into a tuple
      PyObject* args = PyTuple_New(2);
      PyTuple_SetItem(args, 0, py_state);
      PyTuple_SetItem(args, 1, py_values);

      // Call the Python function
      PyObject* result = PyObject_CallObject(py_evaluate_func, args);

      if (result == nullptr) 
      {
        PyErr_Print();  // Print Python error if any
        throw std::runtime_error("Python function call failed.");
      } 
      else 
      {
        // Update values with results from Python function
        for (integer i = 0; i < numOps; ++i)
          values[i] = PyFloat_AsDouble(PyList_GetItem(py_values, i));

        Py_DECREF(result);  // Decrement reference count of result
      }

      // Clean up references
      Py_DECREF(args);
    }
    else
    {
      throw std::runtime_error("No valid Python function set or the function is not callable.");
    }

    // Derivatives
    if (py_derivative_func && PyCallable_Check(py_derivative_func))
    {
      // Prepare Python arguments (state and derivatives as lists)
      PyObject* py_state = PyList_New(numDims);
      PyObject* py_derivatives = PyList_New(numOps);

      for (integer i = 0; i < numDims; ++i)
        PyList_SetItem(py_state, i, PyFloat_FromDouble(state[i]));

      for (integer i = 0; i < numOps; ++i)
      {
        PyObject* py_derivative = PyList_New(numDims);
        for (integer j = 0; j < numDims; ++j)
        {
          PyList_SetItem(py_derivative, j, PyFloat_FromDouble(derivatives(i, j)));
        }
        PyList_SetItem(py_derivatives, i, py_derivative);
      }

      // Pack arguments into a tuple
      PyObject* args = PyTuple_New(2);  // Only two arguments for derivatives
      PyTuple_SetItem(args, 0, py_state);
      PyTuple_SetItem(args, 1, py_derivatives);

      // Call the Python derivative function
      PyObject* result = PyObject_CallObject(py_derivative_func, args);

      if (result == nullptr)
      {
        PyErr_Print();
        throw std::runtime_error("Python derivative function call failed.");
      }
      else
      {
        // Update derivatives with results from Python function
        for (integer i = 0; i < numOps; ++i)
        {
          PyObject* py_derivative = PyList_GetItem(py_derivatives, i);
          for (integer j = 0; j < numDims; ++j)
          {
            derivatives(i, j) = PyFloat_AsDouble(PyList_GetItem(py_derivative, j));
          }
        }
        Py_DECREF(result);
      }

      // Clean up references
      Py_DECREF(args);
    }
  }
};

} /* namespace geos */

#endif /* GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_ */