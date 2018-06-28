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

/*
 * RestartFile.h
 *
 *  Created on: Jul 18, 2012
 *      Author: settgast1
 */

#ifndef RESTARTFILE_H_
#define RESTARTFILE_H_

#include "PMPIOBase.h"
#include "BinStream.h"

class RestartFile : public PMPIOBase<BinStream>
{
public:
  RestartFile();
  ~RestartFile();

  void Initialize(const PMPIO_iomode_t readwrite);

};

#endif /* RESTARTFILE_H_ */
