/*
 * stackTrace.hpp
 *
 *  Created on: May 31, 2016
 *      Author: settgast
 */

#ifndef SRC_CODINGUTILITIES_STACKTRACE_HPP_
#define SRC_CODINGUTILITIES_STACKTRACE_HPP_

namespace geosx
{
namespace stacktrace
{

void handler(int sig, int exitFlag=1, int exitCode=1 );

inline void handler0(int sig)
{
  handler(sig,0,1);
}

inline void handler1(int sig)
{
  handler(sig,1,1);
}


}
}




#endif /* SRC_CODINGUTILITIES_STACKTRACE_HPP_ */
