/* Written 2012 by Matthias S. Benkmann
 *
 * The author hereby waives all copyright and related rights to the contents
 * of this example file (testodr2.cc) to the extent possible under the law.
 */

/**
 * @file
 * @brief Test for multiple definition errors.
 * @copydoc testodr1.cc
 */

#include <cstdio>
#include "../../optionparser/src/optionparser.h"

using option::Descriptor;
using option::Arg;
enum OptionIndex {CREATE};
enum OptionType {DISABLE, ENABLE, OTHER};

extern const Descriptor usage[] = {
   { CREATE, OTHER,
     "c", "create",
     Arg::None,
     "--create\t\t\tTells the program to create something."
   }
 };

extern bool foo(int argc, const char* argv[])
{
  printUsage(std::fwrite, stdout, usage);
  option::Stats stats(usage, argc, argv);
  option::Option buffer [stats.buffer_max];
  option::Option options[stats.options_max];
  option::Parser parse(usage, argc, argv, options, buffer);
  return parse.error();
}

