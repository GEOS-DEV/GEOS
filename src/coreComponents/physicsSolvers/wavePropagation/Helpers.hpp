#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_HELPERS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_HELPERS_HPP_

#include <string>

namespace helpers
{
using namespace std;

string evar( string const & key );
int ienv( string const & key );
bool benv( string const & key );
void print_stacktrace( void );
}

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_HELPERS_ */
