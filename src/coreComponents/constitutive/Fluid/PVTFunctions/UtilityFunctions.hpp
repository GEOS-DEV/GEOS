
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

/**
  * @file TableFunction.hpp
  */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_UTILITYFUNCTION_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_UTILITYFUNCTION_HPP

#include "common/DataTypes.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"

namespace geosx
{

namespace PVTProps
{

  template< typename T >
  using array1dT = std::vector<T>;
  
  using real64_vector = std::vector<geosx::real64>;  
  using string_vector = std::vector<std::string>;  

  static constexpr geosx::localIndex MAX_VAR_DIM = 10;  

  template<typename T, int Dim>
  class EvalArgs 
  {
  public:

    EvalArgs() 
      {
        m_var = 0.0;
        for (int i = 0; i < Dim; ++i)
        {
          m_der[i] = 0;
        }
      }

    ~EvalArgs() = default;

    EvalArgs(const EvalArgs& arg)
      {
        m_var = arg.m_var;
	      for (int i = 0; i < Dim; ++i)
        {
          m_der[i] = arg.m_der[i];
        }
      }

    EvalArgs(const T& var)
      {
        m_var = var;
        for (int i = 0; i < Dim; ++i)
        {
          m_der[i] = 0;
        }
      }

    EvalArgs& operator+=(const EvalArgs& arg)
      {
        this->m_var += arg.m_var;
        for (localIndex i = 0; i < Dim; ++i)
            this->m_der[i] += arg.m_der[i];

        return *this;
      }


    EvalArgs& operator-=(const EvalArgs& arg)
      {
        this->m_var -= arg.m_var;
        for (localIndex i = 0; i < Dim; ++i)
            this->m_der[i] -= arg.m_der[i];

        return *this;
      }

    EvalArgs& operator*=(const EvalArgs& arg)
      {
        const T& u = this->m_var;
        const T& v = arg.m_var;
        for(localIndex i = 0; i < Dim; ++i) 
	  {
            const T& uDer = this->m_der[i];
            const T& vDer = arg.m_der[i];

            this->m_der[i] = (v * uDer + u * vDer);
	  }

        this->m_var *= v;

        return *this;
    }


    EvalArgs& operator/=(const EvalArgs& arg)
      {
        const T& u = this->m_var;
        const T& v = arg.m_var;
        for(localIndex i = 0; i < Dim; ++i) {
	  const T& uDer = this->m_der[i];
	  const T& vDer = arg.m_der[i];

	  this->m_der[i] = (v * uDer - u * vDer) /(v * v);
        }
      
        this->m_var /= v;

        return *this;
      }

    EvalArgs& Exponent(const EvalArgs& arg)
      {
        const T& v = arg.m_var;
        this->m_var = exp(v);
        for(localIndex i = 0; i < Dim; ++i) {
	  const T& vDer = arg.m_der[i];
	  this->m_der[i] = this->m_var * vDer;
        }

        return *this;
      }

    EvalArgs operator+(const EvalArgs& arg) const
      {
        EvalArgs result(*this);
        result += arg;
        return result;
      }

    EvalArgs operator-(const EvalArgs& arg) const
      {
        EvalArgs result(*this);
        result -= arg;
        return result;
      }

    EvalArgs operator-() const
      {
	EvalArgs result;
        result.m_var = -this->m_var;
        for(localIndex i = 0; i < Dim; ++i)
	  result.m_der[i] = - this->m_der[i];

        return result;
      }

    EvalArgs operator*(const EvalArgs& arg) const
      {
        EvalArgs result(*this);
        result *= arg;
        return result;
      }

    EvalArgs operator/(const EvalArgs& arg) const
      {
        EvalArgs result(*this);
        result /= arg;
        return result;
      }

    EvalArgs& operator=(const EvalArgs& arg)
      {
        this->m_var = arg.m_var;
	      for (int i = 0; i < Dim; ++i)
        {
          m_der[i] = arg.m_der[i];
        }

        return *this;
      }

    bool operator==(const EvalArgs& arg) const
      {
	if(this->m_var != arg._m_var)
	  return 0;

        for(localIndex i = 0; i < Dim; ++i)
	  if(this->m_der[i] != arg.m_der[i])
	    return 0;

        return 1;
      }

    bool operator!=(const EvalArgs& arg) const
      {
	return !operator==(arg);

      }

    bool operator>(const EvalArgs& arg) const
      {
	return this->m_var > arg.m_var;
      }

    bool operator<(const EvalArgs& arg) const
      {
	return this->m_var < arg.m_var;
      }

    bool operator>=(const EvalArgs& arg) const
      {
	return this->m_var >= arg.m_var;
      }
    bool operator<=(const EvalArgs& arg) const
      {
	return this->m_var <= arg.m_var;
      }

    EvalArgs& operator+=(const T& arg)
      {
        this->m_var += arg;
        return *this;
      }


    EvalArgs& operator-=(const T& arg)
      {
        this->m_var -= arg;

        return *this;
      }

    EvalArgs& operator*=(const T& arg)
      {
        for(localIndex i = 0; i < Dim; ++i) 
	  {
            this->m_der[i] *= arg;
	  }

        this->m_var *= arg;

        return *this;
    }

    EvalArgs& operator/=(const T& arg)
      {
        for(localIndex i = 0; i < Dim; ++i) 
	  {
            this->m_der[i] /= arg;
	  }

        this->m_var /= arg;

        return *this;
    }

    EvalArgs operator+(const T& arg) const
      {
        EvalArgs result(*this);
        result += arg;
        return result;
      }

    EvalArgs operator-(const T& arg) const
      {
        EvalArgs result(*this);
        result -= arg;
        return result;
      }

    EvalArgs operator*(const T& arg) const
      {
        EvalArgs result(*this);
        result *= arg;
        return result;
      }

    EvalArgs operator/(const T& arg) const
      {
        EvalArgs result(*this);
        result /= arg;
        return result;
      }

    EvalArgs& operator=(const T& arg)
      {
        m_var = arg;
        for (int i = 0; i < Dim; ++i)
        {
          m_der[i] = 0;
        }
        return *this;
      }

    bool operator==(const T& arg) const
      {
	if(this->m_var != arg)
	  return 0;

        return 1;
      }

    bool operator!=(const T& arg) const
      {
	return !operator==(arg);

      }

    bool operator>(const T& arg) const
      {
	return this->m_var > arg;
      }

    bool operator<(const T& arg) const
      {
	return this->m_var < arg;
      }

    bool operator>=(const T& arg) const
      {
	return this->m_var >= arg;
      }
    bool operator<=(const T& arg) const
      {
	return this->m_var <= arg;
      }

    T m_var;
    T m_der[Dim];
    
  };

template<class T, int Dim>
inline EvalArgs<T, Dim> operator+(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    EvalArgs<T, Dim> result(arg2);

    result += arg1;

    return result;
  }

template<class T, int Dim>
inline EvalArgs<T, Dim> operator-(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    EvalArgs<T, Dim> result(arg2);

    result.m_var = arg1 - result.m_var;

    for(localIndex i = 0; i < Dim; ++i)
      {
	result.m_der[i] = -result.m_der[i];
      }

    return result;

  }

template<class T, int Dim>
inline EvalArgs<T, Dim> operator*(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    EvalArgs<T, Dim> result(arg2);

    result *= arg1;

    return result;
  }

template<class T, int Dim>
inline EvalArgs<T, Dim> operator/(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    EvalArgs<T, Dim> result(arg2);

    T coef = - (arg1 / result.m_var / result.m_var);
    result.m_var = arg1 / result.m_var;

    for(localIndex i = 0; i < Dim; ++i)
      {
	result.m_der[i] = coef * result.m_der[i];
      }

    return result;
  }


template<class T, int Dim>
inline bool operator==(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    return (arg2 == arg1);

  }

template<class T, int Dim>
inline bool operator!=(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    return (arg2 != arg1);

  }

template<class T, int Dim>
inline bool operator>(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    return arg1 > arg2.m_var;
  }

template<class T, int Dim>
inline bool operator<(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    return arg1 < arg2.m_var;
  }

template<class T, int Dim>
inline bool operator>=(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    return arg1 >= arg2.m_var;
  }


template<class T, int Dim>
inline bool operator<=(const T& arg1, const EvalArgs<T, Dim> &arg2)
  {
    return arg1 <= arg2.m_var;
  }

  typedef EvalArgs<geosx::real64, 1> EvalArgs1D;
  typedef EvalArgs<geosx::real64, 2> EvalArgs2D;
  typedef EvalArgs<geosx::real64, 3> EvalArgs3D;
  typedef EvalArgs<geosx::real64, MAX_VAR_DIM> EvalVarArgs;

  class TableFunctionBase 
  {
  public:
    virtual const string& TableName() const = 0;
    virtual ~TableFunctionBase(){}
    
    virtual EvalArgs1D Value(const EvalArgs1D& x) const = 0;

    virtual EvalArgs2D Value(const EvalArgs2D& x) const = 0;

    virtual EvalArgs2D Value(const EvalArgs2D& x, const EvalArgs2D& y) const = 0;

    virtual void Print() const = 0;

  };  


  typedef std::shared_ptr<TableFunctionBase> TableFunctionPtr;

  class XYTable : public TableFunctionBase
  {
  public:

    XYTable(const std::string& tableName, const real64_vector& x, const real64_vector& y, const array1dT<real64_vector> &value) : m_tableName(tableName), m_x(x), m_y(y), m_value(value) {}

    ~XYTable(){}

    real64_vector &XArray() {
      return m_x;
    }

    real64_vector &YArray() 
    {
      return m_y;
    }

    array1dT<real64_vector> &ValueArray() 
    {
      return m_value;
    }


    virtual const string& TableName() const 
    {
      return m_tableName;
    }


    virtual EvalArgs1D Value(const EvalArgs1D& x) const
    {
      return 0;
    }

    virtual EvalArgs2D Value(const EvalArgs2D& x) const 
    {
      return 0;
    }

    virtual EvalArgs2D Value(const EvalArgs2D& x, const EvalArgs2D& y) const;

    virtual void Print() const
    {
    }

  private:

    std::string m_tableName;
    real64_vector m_x;
    real64_vector m_y;
    array1dT<real64_vector> m_value;
  
  };

  class XTable : public TableFunctionBase
  {
  public:

    XTable(const string& tableName, const real64_vector& x, const real64_vector &value) : m_tableName(tableName), m_x(x), m_value(value) {}
    ~XTable(){}

    real64_vector &XArray()
    {
      return m_x;
    }

    real64_vector &ValueArray()
    {
      return m_value;
    }

    virtual const string& TableName() const 
    {
      return m_tableName;
    }

    virtual EvalArgs1D Value(const EvalArgs1D& x) const 
    {

      return GetValue<EvalArgs1D>(x);

    }

    virtual EvalArgs2D Value(const EvalArgs2D& x) const 
    {

      return GetValue<EvalArgs2D>(x);

    }

    virtual EvalArgs2D Value(const EvalArgs2D& x, const EvalArgs2D& y) const  
    {
      return 0;
    }

  private:

    template<class T> 
    T GetValue(const T& x) const;

    virtual void Print() const
    {
    }

    string m_tableName;
    real64_vector m_x;
    real64_vector m_value;

  };

}

}  
#endif
