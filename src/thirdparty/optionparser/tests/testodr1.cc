/* Written 2012 by Matthias S. Benkmann
 *
 * The author hereby waives all copyright and related rights to the contents
 * of this example file (testodr1.cc) to the extent possible under the law.
 */

/**
 * @file
 * @brief Test for multiple definition errors.
 *
 * @note
 * This program is for developing TLMC++OP. It is neither an example nor a functionality test.
 * Do not worry if it doesn't compile or run on your platform.
 *
 * @ref testodr1.cc and @ref testodr2.cc test optionparser.h for
 * violations of the one definition rule, both at compile-time and at
 * link-time. IOW, they test if optionparser.h can be included
 * multiple times as well as that multiple object files that include
 * it can be linked together.
 *
 */

#include <cstdio>
#include "../../optionparser/src/optionparser.h" //intentionally included twice
#include "../../optionparser/src/optionparser.h"

using option::Option;
using option::Descriptor;

extern const Descriptor usage[];

extern bool bar(int argc, const char* argv[])
{
  printUsage(std::fwrite, stdout, usage);
  option::Stats stats(usage, argc, argv);
  option::Option buffer [stats.buffer_max];
  option::Option options[stats.options_max];
  option::Parser parse(usage, argc, argv, options, buffer);
  return parse.error();
}

int main()
{
  Descriptor d = usage[0];
  std::printf("%s",d.shortopt);
}



