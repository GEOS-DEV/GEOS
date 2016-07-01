#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <cxxabi.h>
#include <sys/ucontext.h>

#include "stackTrace.hpp"

namespace geosx
{
namespace stacktrace
{

void handler(int sig, int exitFlag, int exitCode )
{
  void *array[20];
  int size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 20);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);

  fprintf(stderr,"attempting unmangled trace: \n");

//  char** mangled_name = backtrace_symbols(array, size);
  char ** messages    = backtrace_symbols(array, size);

  // skip first stack frame (points here)
  for ( int i = 1; i < size && messages != NULL; ++i)
  {
      char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;

      // find parantheses and +address offset surrounding mangled name
      for (char *p = messages[i]; *p; ++p)
      {
          if (*p == '(')
          {
              mangled_name = p;
          }
          else if (*p == '+')
          {
              offset_begin = p;
          }
          else if (*p == ')')
          {
              offset_end = p;
              break;
          }
      }

      // if the line could be processed, attempt to demangle the symbol
      if (mangled_name && offset_begin && offset_end &&
          mangled_name < offset_begin)
      {
          *mangled_name++ = '\0';
          *offset_begin++ = '\0';
          *offset_end++ = '\0';

          int status;
          char * real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);

          // if demangling is successful, output the demangled function name
          if (status == 0)
          {
              std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
                        << real_name << "+" << offset_begin << offset_end
                        << std::endl;

          }
          // otherwise, output the mangled function name
          else
          {
              std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
                        << mangled_name << "+" << offset_begin << offset_end
                        << std::endl;
          }
          free(real_name);
      }
      // otherwise, print the whole line
      else
      {
          std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
      }
  }
  std::cerr << std::endl;

  free(messages);
  if( exitFlag == 1)
    exit(exitCode);

}



}


}

