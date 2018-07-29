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
 * SiloOutput.hpp
 *
 *  Created on: Jan 26, 2018
 *      Author: sherman
 */

#ifndef SRC_COMPONENTS_CORE_SRC_SILOOUTPUT_HPP_
#define SRC_COMPONENTS_CORE_SRC_SILOOUTPUT_HPP_

#include "OutputBase.hpp"


namespace geosx
{

class SiloOutput : public OutputBase
{
public:
  SiloOutput( std::string const & name,
              ManagedGroup * const parent );

  virtual ~SiloOutput() override;

  static string CatalogName() { return "Silo"; }

  virtual void FillDocumentationNode() override;

  virtual void Execute( real64 const & time_n,
                        real64 const & dt,
                        int const cycleNumber,
                        dataRepository::ManagedGroup * domain ) override;

  struct viewKeyStruct
  {
    dataRepository::ViewKey plotFileRoot = { "plotFileRoot" };
    dataRepository::ViewKey writeFEMFaces = { "writeFEMFaces" };
  } viewKeys;

};


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_SILOOUTPUT_HPP_ */
