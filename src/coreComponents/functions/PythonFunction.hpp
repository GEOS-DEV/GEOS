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

  PythonFunction(): py_evaluate_func(nullptr)
  {
    // Initialize Python if needed (ensure Python interpreter is initialized)
    if (!Py_IsInitialized()) 
      Py_Initialize();
  }

  ~PythonFunction() 
  {
    Py_XDECREF(py_evaluate_func);  // Decrement reference count of the function
    if (Py_IsInitialized()) 
    {
      Py_Finalize();
    }
  }

  void setEvaluateFunction(PyObject* pyFunc) 
  {
    if (PyCallable_Check(pyFunc)) 
    {
      Py_XINCREF(pyFunc);  // Increment reference count for the new function
      Py_XDECREF(py_evaluate_func);  // Decrement reference count of the old function
      py_evaluate_func = pyFunc;  // Assign the new function
    } 
    else 
    {
      throw std::runtime_error("Provided Python object is not callable.");
    }
  }

  GEOS_HOST_DEVICE   
  inline  
  void 
  evaluate(const stackArray1d<real64, numDims>& state, stackArray1d<real64, numOps>& values) const
  {
    if (py_evaluate_func && PyCallable_Check(py_evaluate_func))
    {
      // Prepare Python arguments (state and values as lists)
      PyObject* py_state = PyList_New(numDims);
      PyObject* py_values = PyList_New(numOps);

      for (integer i = 0; i < numDims; ++i) {
        PyList_SetItem(py_state, i, PyFloat_FromDouble(state[i]));
      }

      for (integer i = 0; i < numOps; ++i) {
        PyList_SetItem(py_values, i, PyFloat_FromDouble(values[i]));
      }

      // Pack arguments into a tuple
      PyObject* args = PyTuple_New(2);
      PyTuple_SetItem(args, 0, py_state);
      PyTuple_SetItem(args, 1, py_values);

      // Call the Python function
      PyObject* result = PyObject_CallObject(py_evaluate_func, args);

      if (result == nullptr) {
        PyErr_Print();  // Print Python error if any
        throw std::runtime_error("Python function call failed.");
      } else {
        // Update values with results from Python function
        for (integer i = 0; i < numOps; ++i) {
          values[i] = PyFloat_AsDouble(PyList_GetItem(py_values, i));
        }
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
};

} /* namespace geos */

#endif /* GEOS_FUNCTIONS_PYTHONFUNCTION_HPP_ */