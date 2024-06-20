
/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ConstexprConstructs.hpp
 */

#ifndef GEOS_CODINGUTILITIES_CONSTEXPRCONSTRUCTS_HPP_
#define GEOS_CODINGUTILITIES_CONSTEXPRCONSTRUCTS_HPP_

#include <type_traits>
#include <utility>

namespace geos
{

namespace compileTime
{

namespace internal
{
/**
 * @struct to_string_t
 * @brief Provides the ability to convert any integral to a string at compile-time.
 * @tparam N Number to convert
 * @tparam base Desired base, can be from 2 to 36
 */
template< auto N, unsigned int base >
struct to_string_t
{
    static_assert(base > 1 && base < 37, "Base must be between 2 and 36");
    static_assert(std::is_integral_v<decltype(N)>, "N must be integral");
    // The lambda calculates what the string length of N will be, so that `buf`
    // fits to the number perfectly including the null terminator.
    char buf[( []
    {
      unsigned int len = N > 0 ? 1 : 2;
      for (auto n = N < 0 ? -N : N; n; len++, n /= base);
      return len;
    } () ) ] = {};

    /**
     * Constructs the object, filling `buf` with the string representation of N.
     */
    constexpr to_string_t() noexcept
    {
        char * ptr = buf + sizeof(buf) - 1;
        *ptr = '\0'; // Set the last character as null terminator
        auto n = N < 0 ? -N : N;
        do {
            *--ptr = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[n % base];
            n /= base;
        } while (n > 0);

        if (N < 0) {
            *--ptr = '-';
        }
    }

    /**
     * Allows implicit conversion of this object to a `char *`.
     */
    constexpr operator char *()
    {
      return buf;
    }

    /**
     * Allows implicit conversion of this object to a `const char *`.
     */
    constexpr operator const char *() const
    {
      return buf;
    }
};

} // namespace internal

/**
 * Simplifies use of `to_string_t` from `to_string_t<N>()` to `to_string<N>`.
 */
template<auto N, unsigned int base = 10>
constexpr internal::to_string_t<N, base> to_string;

namespace internal
{
template<std::size_t N, std::size_t M>
struct strcat_result
{
  char buf[N + M - 1]; // -1 to not double count null terminator
  /**
   * Allows implicit conversion of this object to a `char *`.
   */
  constexpr operator char *()
  {
    return buf;
  }
  /**
   * Allows implicit conversion of this object to a `const char *`.
  */
  constexpr operator const char *() const
  {
    return buf;
  }
};
}

template<std::size_t N, std::size_t M>
constexpr internal::strcat_result<N, M> strcat(const char (&a)[N], const char (&b)[M])
{
  internal::strcat_result<N, M> result{};
  for (std::size_t i = 0; i < N - 1; ++i)
  {
    result.buf[i] = a[i];
  }
  for (std::size_t j = 0; j < M; ++j)
  {
    result.buf[N - 1 + j] = b[j];
  }
  return result;
}

template<int Start, int End, int Inc = 1, typename Func>
constexpr void static_for(Func && func)
{
  if constexpr (Start < End)
  {
    constexpr std::integral_constant<int, Start> idx{};
    func( idx );
    static_for<Start + Inc, End, Inc>(std::forward<Func>(func));
  }
}

} // namespace compileTime

} // namespace geos

#endif