// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#ifndef RESTARTFLAGS_H_
#define RESTARTFLAGS_H_

namespace geosx
{
namespace dataRepository
{

enum class RestartFlags : unsigned char
{
  NO_WRITE,
  WRITE,
  WRITE_AND_READ
};

}   /* namespace dataRepository */
}   /* namespace geosx */

#endif  /* RESTARTFLAGS_H_ */
