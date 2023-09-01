#include "Helpers.hpp"

#include <execinfo.h> /* backtrace, backtrace_symbols_fd */
#include <unistd.h> /* STDOUT_FILENO */
#include <cstring>

namespace helpers
{
using namespace std;

string evar( string const & key )
{
  char const * val = getenv( key.c_str());
  return val == NULL ? string() : string( val );
}
int ienv( string const & key )
{
  string val = evar( key );
  return atoi( val.c_str());
}
bool benv( string const & key )
{
  string val = evar( key );
  char const *ptr = val.c_str();
  return strcmp( ptr, "true" ) == 0 || strcmp( ptr, "1" ) == 0 || strcmp( ptr, "T" ) == 0;
}

void print_stacktrace( void )
{
  size_t size;
  enum Constexpr { MAX_SIZE = 1024 };
  void *array[MAX_SIZE];
  size = backtrace( array, MAX_SIZE );
  backtrace_symbols_fd( array, size, STDOUT_FILENO );
}

}
