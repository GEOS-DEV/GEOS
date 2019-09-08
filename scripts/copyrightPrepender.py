#!/usr/bin/env python

# ----------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2018-2019 Total, S.A
# Copyright (c) 2019-     GEOSX Contributors
# All right reserved
#
# See toplevel LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
# ----------------------------------------------------------------------------


# Python script to add GEOSX copyright notice at top of files in a directory
#
# The script takes a directory and checks all files in the directory.
# If the second line in the file does not match the copyright string, we add it.
#
# Modified from an initial script by P. Sinha
# Modified on 2019-SEP-06 by Herve Gross

import os
import sys
import argparse
import numpy as np


new_copyright_str = \
    """/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

"""


old_copyright_str = \
    """/*
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

"""


old_copyright_str = \
    """/*
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

"""

new_copyright_str_arr = new_copyright_str.split("\n")[:-1]

old_copyright_str_arr = old_copyright_str.split("\n")[:-1]

max_copyright_lines = max(len(new_copyright_str_arr), len(old_copyright_str_arr))


def striplist(my_list):
    return([x.strip() for x in my_list])


def is_pattern_in_list(pattern, my_list):
    return True
    return any(my_list[i:i + len(pattern)] == pattern for i in xrange(len(my_list) - len(pattern) + 1))


def find_pattern(pattern, my_list):
    matches = []
    # Left here for debugging purposes
    # for i in range(0,1):
    # for j,p in enumerate(pattern):
    # res = my_list[j] == p
    # print('{} \t {} | {}'.format(res, my_list[j], p))

    for i in range(len(my_list)):
        if my_list[i] == pattern[0] and my_list[i:i + len(pattern)] == pattern:
            matches.append(i)
    return matches


def locate_block_of_text(pattern, file_content):
    line_numbers = []
    if is_pattern_in_list(pattern, file_content):
        line_number_where_pattern_starts = find_pattern(pattern, file_content)
        for starting_line in line_number_where_pattern_starts:
            line_numbers.extend(range(starting_line, starting_line + len(pattern)))
    return line_numbers


def delete_indices_from_list(my_list_of_indices, indices_to_delete):
    arr = np.array(my_list_of_indices, dtype='int32')
    return list(np.delete(arr, indices_to_delete))


def locate_string(pattern, file_content):
    ids = []
    for id, line in enumerate(file_content):
        if pattern in line:
            ids.append(id)
    return ids


def process_file(file_name, test_only):
    """
    Return process code:
    0: old copyright found and removed, new copyright not found, new copyright prepended [intended purpose]
    1: old copyright found and removed, new copyright not found, new copyright not prepended
    2: old copyright found and kept, new copyright not found, new copyright prepended
    3: old copyright found and kept, new copyright not found, new copyright not prepended [FILE UNMODIFIED]
    4: old copyright not found, new copyright not found, new copyright prepended
    5: old copyright not found, new copyright not found, new copyright not prepended [FILE UNMODIFIED]
    6: old copyright not found, new copyright found [FILE UNMODIFIED] - file is good already
    -1: error somewhere, file returned unmodified
    """

    process_code = -1

    with open(file_name, 'r') as f:
        indented_content = f.read().splitlines()

    file_content = striplist(indented_content)
    old_copyright = striplist(old_copyright_str_arr)
    new_copyright = striplist(new_copyright_str_arr)

    old_copyright_lines = locate_block_of_text(old_copyright, file_content)
    new_copyright_lines = locate_block_of_text(new_copyright, file_content)

    want_to_remove_old_copyright = True
    want_to_prepend_new_copyright = True
    found_old_copyright = bool(old_copyright_lines)
    found_new_copyright = bool(new_copyright_lines)

    if want_to_remove_old_copyright:
        lines_to_keep = delete_indices_from_list(range(len(file_content)), old_copyright_lines)
    else:
        lines_to_keep = range(len(file_content))  # Keep all lines from the original file

    # Find abnormal patterns (like a Copyright mention outside of the original header)
    abnormal_word = 'Copyright'
    found_abnormal_word = False
    lines_containing_abnormal_word = locate_string(abnormal_word, file_content)
    for line in lines_containing_abnormal_word:
        if line in lines_to_keep:
            found_abnormal_word = True

    # Get process code:
    if found_old_copyright and want_to_remove_old_copyright and not found_new_copyright and want_to_prepend_new_copyright:
        process_code = 0
    elif found_old_copyright and want_to_remove_old_copyright and not found_new_copyright and not want_to_prepend_new_copyright:
        process_code = 1
    elif found_old_copyright and not want_to_remove_old_copyright and not found_new_copyright and want_to_prepend_new_copyright:
        process_code = 2
    elif found_old_copyright and not want_to_remove_old_copyright and not found_new_copyright and not want_to_prepend_new_copyright:
        process_code = 3  # file unmodified
    elif not found_old_copyright and not found_new_copyright and want_to_prepend_new_copyright:
        process_code = 4
    elif not found_old_copyright and not found_new_copyright and not want_to_prepend_new_copyright:
        process_code = 5  # file unmodified
    elif not found_old_copyright and found_new_copyright:
        process_code = 6  # file unmodified
    elif found_abnormal_word:
        process_code = 7

    # Now do the file rewrite, if desired
    if not test_only:
        if process_code in [3, 5, 6]:
            pass
        else:
            with open(file_name, 'w') as f_out:
                f_out.seek(0)
                if process_code in [0, 2, 4]:  # prepending copyrights
                    f_out.write(new_copyright_str)
                for i in lines_to_keep:
                    f_out.write(indented_content[i] + '\n')

    return process_code


def fileNameGenerator(rootDir, validExtensions, isRecursive=False):
    """
    Generator function for file names whose extensions are in the validExtensions tuple.
    Files are rooted in rootDir, and process is recursive if isRecursive==True.
    """
    if isRecursive:
        for path, dirlist, filelist in os.walk(rootDir):
            for f in (f for f in filelist if f.lower().endswith(validExtensions)):
                yield os.path.join(path, f)
    else:
        for f in os.listdir(rootDir):
            if f.lower().endswith(validExtensions):
                yield os.path.join(rootDir, f)


if __name__ == "__main__":

    # -----------------
    # Argument parsing
    # -----------------

    # Setup the argument parser, dir is required first argument
    parser = argparse.ArgumentParser(description="Append LLNL copyright message to files.")
    # TODO -- should we accept multiple directories?
    parser.add_argument("dir", type=str, help="specify directory containing files on which we want to operate.")

    # Option to recursively search for files
    parser.add_argument("-r", "--recursive", dest='isRecursive', action='store_true',
                        help="add flag to recursively descend to subdirectories.")
    parser.add_argument("--no-recursive", dest='isRecursive', action='store_false')
    parser.set_defaults(isRecursive=False)

    # Test run to see which files might require a copyright notice.
    parser.add_argument("-t", "--test", action='store_true',
                        help="add flag if we only want to see which files are missing copyright statements, but not to modify any files.")

    # Additional possible featues to be implemented
    #
    # Specify a collection of file types.
    # parser.add_argument("-f", "--filetypes", type=str, nargs="*", help="add file types on which to operate ")    # default should be *.hpp and *.cpp
    # For now, we hardcode to c/c++ header and source files
    valid_extensions = (".hpp", ".cpp", ".h", ".c", ".cxx", ".cc")
    args = parser.parse_args()

    # -----------------
    # Processing
    # -----------------

    # Search for files to inspect in the base directory (recursively if using the -r command line token)
    print("Work directory {}".format(args.dir))
    file_name_generator = fileNameGenerator(args.dir, valid_extensions, args.isRecursive)

    # Make a list of files from the generator (helpful to restrict the number of files while debugging)
    file_list = list(file_name_generator)
    print('Inspecting {} files'.format(len(file_list)))

    # Store information for analysis purposes
    process_codes = []  # Integer value that reflects the changes proposed or done to the file
    file_names = []  # File name

    for file_name in file_list:
        process_code = process_file(file_name, args.test)
        process_codes.append(process_code)
        file_names.append(file_name)

    for list_this_code in [0, 1, 2, 3, 4, 5, 6, -1]:
        for p, f in zip(process_codes, file_names):
            if p == list_this_code:
                print('{} | {}'.format(p, f))

    print('Done')
