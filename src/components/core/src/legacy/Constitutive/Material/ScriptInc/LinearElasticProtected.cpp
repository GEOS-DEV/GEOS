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

//FUNCTION_BEGIN_PARSE
virtual_void
LinearElastic::PreSetValues(const array<string>& names)
{
  //reset mechanical constants
  for (localIndex a = 0 ; a < m_parameterData.size() ; ++a)
  {
    m_parameterData[a].E = 0;
    m_parameterData[a].Lame = 0;
    m_parameterData[a].Nu = 0;
    m_parameterData[a].init_bulkModulus = 0;
    m_parameterData[a].init_shearModulus = 0;
  }
}

//FUNCTION_BEGIN_PARSE
virtual_void
LinearElastic::PostSetValues(const array<string>& names)
{
  //recalculate mechanical parameters that were not already explicitly set
  for(localIndex a = 0 ; a < m_parameterData.size() ; a++)
  {
    realT M = 0;
    FillLinearElasticModuli(m_parameterData[a].init_bulkModulus,
                            m_parameterData[a].init_shearModulus,
                            m_parameterData[a].E,
                            m_parameterData[a].Nu,
                            m_parameterData[a].Lame,
                            M);
  }
}
