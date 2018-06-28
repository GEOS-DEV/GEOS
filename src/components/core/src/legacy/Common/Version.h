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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file Version.h
 * @author Stuart Walsh
 * @date created on Apr 25, 2014
 *
 */

#ifndef _VERSION_H_
#define _VERSION_H_

#include <cstdlib>
#include <iostream>
#include <string>

// This allows the makefile to set the version during compilation
#ifndef REPO_VERSION
#define REPO_VERSION "Unknown Repository Version"
#endif

#ifndef EXTERNAL_REPO_VERSION
#define EXTERNAL_REPO_VERSION ""
#endif

#ifndef INTERNAL_REPO_VERSION
#define INTERNAL_REPO_VERSION ""
#endif


void DisplayVersion(int rank);
void DisplaySplash(int rank);

void DisplayVersionHistory(int rank);

std::string GetRepoVersionString();


#endif /* _VERSION_H_ */
