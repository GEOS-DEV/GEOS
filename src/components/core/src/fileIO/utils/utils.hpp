/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifndef SRC_COMPONENTS_CORE_SRC_FILEIO_UTILS_UTILS_HPP
#define SRC_COMPONENTS_CORE_SRC_FILEIO_UTILS_UTILS_HPP

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <vector>
#include <regex>

namespace geosx
{

/* Taken from http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html */
inline void readDirectory(const std::string& name, std::vector<std::string>& v)
{
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != nullptr)
  {
    v.push_back(dp->d_name);
  }
  closedir(dirp);
}


inline void getAbsolutePath(const std::string & path, std::string & absolute_path )
{
  char abs_file_path[PATH_MAX + 1];
  if (realpath(path.data(), abs_file_path))
  {
    absolute_path = abs_file_path;
  }
  else
  {
    getcwd(abs_file_path, PATH_MAX + 1);
    GEOS_ERROR("Could not get the absolute path for " << path << " from "  << abs_file_path);
  }
}


inline void splitPath(const std::string& path, std::string& dirname, std::string& basename)
{
  size_t pos = path.find_last_of('/');
  if (pos == string::npos)
  {
    dirname = std::string(".");
    basename = path;
  }

  if (pos == 0)
  {
    dirname = std::string("/");
    basename = path.substr(1);
  }

  dirname = path.substr(0, pos);
  basename = path.substr(pos + 1);
}


template <typename REGEX>
inline bool regexMatch(const std::string & str, REGEX regex)
{
  return std::regex_match(str, regex);
}

} /* end namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FILEIO_UTILS_UTILS_HPP */
