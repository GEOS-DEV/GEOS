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
 * StableTimeStep.h
 *
 *  Created on: Jan 24, 2012
 *      Author: scottjohnson
 */

#ifndef STABLETIMESTEP_H_
#define STABLETIMESTEP_H_

struct StableTimeStep
{
//  std::string m_region;
//  localIndex m_index;
  realT m_maxdt;

  inline bool SetIfSmaller(const realT dt)
  {
    if(dt < m_maxdt)
    {
      this->m_maxdt = dt;
      return true;
    }
    else
    {
      return false;
    }
  }
};

#endif /* STABLETIMESTEP_H_ */
