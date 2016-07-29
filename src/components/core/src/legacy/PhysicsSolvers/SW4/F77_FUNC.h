/*
 * F77_FUNC.h
 *
 *  Created on: Apr 28, 2014
 *      Author: walsh24
 */

#ifndef F77_FUNC_H_
#define F77_FUNC_H_

/*--------------------------------------------------------------------------*/

/* F77_FUNC - attempt a uniform naming of FORTRAN functions which
 *          - gets around loader naming conventions
 *          - F77_FUNC(foo, FOO)(x, y, z)
 */

#ifndef F77_FUNC

/* MACOSX predefines __APPLE__ */
#ifdef __APPLE__
#  define F77_FUNC(x, X) x
#endif

# ifdef ANSI_F77
#  define F77_FUNC(x, X)  X
# endif /* ANSI_F77 */

# ifndef __GNUC__

#  ifdef __xlC__
#   define F77_FUNC(x, X)  x
#  endif /* IBM XL compiler */

#  ifdef HPUX
#   define F77_FUNC(x, X)  x
#  endif /* HPUX */

# endif /* __GNUC__ */

# ifndef F77_FUNC
#  define F77_FUNC(x, X)  x ## _
# endif /* F77_FUNC */

#endif /* F77_FUNC */


#endif /* F77_FUNC_H_ */
