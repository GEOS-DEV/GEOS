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
 * @file sfinae.hpp
 */

#ifndef SRC_CODINGUTILITIES_SFINAE_HPP_
#define SRC_CODINGUTILITIES_SFINAE_HPP_

#include "cxx-utilities/src/Macros.hpp"


/** This macro creates a struct that has a static member "value" which returns
   true if the typename "TT" has a
 *  datamember "NAME".
 */
#define HAS_MEMBER_DATA( NAME ) \
  template<typename TT> \
  struct has_datamember_ ## NAME \
  { \
private: \
    template<typename U> static constexpr auto test( int )->decltype( std::declval<U>().NAME, bool() ) \
    { \
      return std::is_member_object_pointer<decltype(&U::NAME)>::value; \
    } \
    template<typename U> static constexpr auto test( ... )->bool \
    { \
      return false; \
    } \
public: \
    static constexpr bool value = test<TT>( 0 ); \
  };

#define HAS_STATIC_MEMBER_DATA( NAME ) \
  template<typename TT> \
  struct has_staticdatamember_ ## NAME \
  { \
private: \
    template<typename U> static constexpr auto test( int )->decltype( U::NAME, bool() ) \
    { \
      return !std::is_function< decltype(U::NAME) >::value && !std::is_member_object_pointer<decltype(&U::NAME)>::value; \
    } \
    template<typename U> static constexpr auto test( ... )->bool \
    { \
      return false; \
    } \
public: \
    static constexpr bool value = test<TT>( 0 ); \
  };



#define HAS_MEMBER_FUNCTION_NAME( NAME ) \
  template<typename TT > \
  struct has_memberfunction_name_ ## NAME \
  { \
private: \
    template<typename U> static constexpr auto test( int )->decltype( std::is_member_function_pointer<decltype(&U::NAME)>::value, bool() ) \
    { \
      return std::is_member_function_pointer<decltype(&U::NAME)>::value; \
    } \
    template<typename U> static constexpr auto test( ... )->bool \
    { \
      return false; \
    } \
public: \
    static constexpr bool value = test<TT>( 0 ); \
  };


#define HAS_MEMBER_FUNCTION_VARIANT( NAME, VARIANT, RTYPE, CONST, PARAMS, ARGS ) \
  template<typename TT > \
  struct has_memberfunction_v ## VARIANT ## _ ## NAME \
  { \
private: \
    template<typename U> static constexpr auto test( int )->decltype( static_cast<RTYPE (U::*)( PARAMS ) CONST>(&U::NAME), bool() ) \
    { \
      return std::is_same< decltype( std::declval<U>().NAME( ARGS ) ), RTYPE>::value; \
    } \
    template<typename U> static constexpr auto test( ... )->bool \
    { \
      return false; \
    } \
public: \
    static constexpr bool value = test<TT>( 0 ); \
  };

#define HAS_MEMBER_FUNCTION( NAME, RTYPE, CONST, PARAMS, ARGS ) \
  template<typename TT > \
  struct has_memberfunction_ ## NAME \
  { \
private: \
    template<typename U> static constexpr auto test( int )->decltype( static_cast<RTYPE (U::*)( PARAMS ) CONST>(&U::NAME), bool() ) \
    { \
      return std::is_same< decltype( std::declval<U>().NAME( ARGS ) ), RTYPE>::value; \
    } \
    template<typename U> static constexpr auto test( ... )->bool \
    { \
      return false; \
    } \
public: \
    static constexpr bool value = test<TT>( 0 ); \
  };


#define HAS_STATIC_MEMBER_FUNCTION( NAME, RTYPE, ... ) \
  template<typename TT> \
  struct has_staticmemberfunction_ ## NAME \
  { \
private: \
    template<typename U> static constexpr auto test( int )->decltype( U::NAME( __VA_ARGS__ ), bool() ) \
    { \
      return std::is_same< decltype( U::NAME( __VA_ARGS__ ) ), RTYPE>::value; \
    } \
    template<typename U> static constexpr auto test( ... )->bool \
    { \
      return false; \
    } \
public: \
    static constexpr bool value = test<TT>( 0 ); \
  };


#define HAS_ENUM( NAME ) \
  template<typename TT > \
  struct has_enum_ ## NAME \
  { \
private: \
    template<typename U> static auto test( int )->typename std::enable_if<std::is_enum< typename U::NAME >::value, std::true_type>::type; \
    template<typename U> static auto test( ... )->std::false_type; \
public: \
    static constexpr bool value = decltype(test<TT>( 0 ))::value; \
  };

#define HAS_ALIAS( NAME ) \
  template<typename TT> \
  struct has_alias_ ## NAME \
  { \
private: \
    template<typename U> static auto test( int )->typename std::enable_if<!std::is_enum< typename U::NAME >::value, std::true_type>::type; \
    template<typename U> static auto test( ... )->std::false_type; \
public: \
    static constexpr bool value = !(std::is_same<decltype(test<TT>( 0 )), std::false_type>::value); \
  };



#define CONDITIONAL_VIRTUAL_FUNCTION( CLASSNAME, FUNCNAME, RTYPE, CONST, PARAMS, ARGS ) \
  template<typename U=T, bool has = has_memberfunction_ ## FUNCNAME<U>::value > \
  struct wrapper ## FUNCNAME \
  { \
    RTYPE f( ... ) {return RTYPE(); } \
  }; \
  template<typename U> \
  struct wrapper ## FUNCNAME<U, true> \
  { \
    RTYPE f( CLASSNAME * const obj, PARAMS ) \
    { \
      return (*obj).m_data->FUNCNAME( ARGS ); \
    } \
  }; \
  virtual RTYPE FUNCNAME( PARAMS ) CONST override final \
  { \
    wrapper ## FUNCNAME<T> temp; \
    return temp.f( this, ARGS ); \
  }

#define CONDITIONAL_VIRTUAL_FUNCTION0( CLASSNAME, FUNCNAME, RTYPE, CONST ) \
  template<typename U=T, bool has = has_memberfunction_ ## FUNCNAME<U>::value > \
  struct wrapper ## FUNCNAME \
  { \
    RTYPE f( ... ) {return RTYPE(); } \
  }; \
  template<typename U> \
  struct wrapper ## FUNCNAME<U, true> \
  { \
    RTYPE f( CLASSNAME CONST * const obj ) \
    { \
      return (*obj).m_data->FUNCNAME(); \
    } \
  }; \
  virtual RTYPE FUNCNAME() CONST override final \
  { \
    wrapper ## FUNCNAME<T> temp; \
    return temp.f( this ); \
  }

#endif /* SRC_CODINGUTILITIES_SFINAE_HPP_ */
