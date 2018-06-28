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
 * @file ParallelPlateFluidModelBase.h
 * @author walsh24
 * @date March 10, 2014
 */

#ifndef PARALLELPLATEFLUIDMODELBASE_H_
#define PARALLELPLATEFLUIDMODELBASE_H_

#include "Utilities/StringUtilities.h"
#include "Common/typedefs.h"
#include "IO/ticpp/TinyXMLParser.h"

#include <map>
#include <string>
#include <vector>

//class ProblemManagerT;

//////////////////////////

// Fluid Model Base class


class ParallelPlateFluidModelBase
{

public:
  ParallelPlateFluidModelBase(){ /** empty **/};
  virtual ~ParallelPlateFluidModelBase(){ /** empty **/};

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn){};
  // two sided permeability
  virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                                      const realT apb,const realT w, const realT qMag, const realT SHP_FCT) = 0;
  // one sided permeability
  virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                      const realT qMag, const realT SHP_FCT) = 0; // one
                                                                                  // sided

};


//////////////////////////

// Fluid Model Factory
//
// Consists of the following parts:
//   * The function to generate new parallelPlateFluidModel pointers:
// "newParallelPlateFluidModel"
//   * A base class to derive the functions to generate parallelPlateFluidModel
// pointers: "ParallelPlateFluidModelInitializer"
//   * A String-to-ParallelPlateFluidModel-Intializer map hidden behind the
// getParallelPlateFluidModelCatalogue function
//   * A template to create parallelPlateFluidModel initializers:
// "ParallelPlateFluidModelRegistrator"
//   * A compiler directive to simplify autoregistration:
// "REGISTER_PARALLEL_PLATE_FLUID_MODEL"
//
// Most ParallelPlateFluidModels will only need to use one or two of the parts:
//   * To register a new ParallelPlateFluidModel in the factory:
// REGISTER_PARALLEL_PLATE_FLUID_MODEL( ModelClassName )
//   * To load a ParallelPlateFluidModel pointer from the factory:
//       ParallelPlateFluidModelBase* aParallelPlateFluidModelPtr =
// ParallelPlateFluidModelFactory::NewParallelPlateFluidModel(ParallelPlateFluidModelString,
// args );

/// Base class to generate new ParallelPlateFluidModel pointers
class ParallelPlateFluidModelInitializer
{
public:
  virtual ParallelPlateFluidModelBase* InitializeParallelPlateFluidModel(TICPP::HierarchicalDataNode* const hdn) = 0;

  virtual ~ParallelPlateFluidModelInitializer()
  {}
};

typedef std::map<std::string, ParallelPlateFluidModelInitializer*> ParallelPlateFluidModelCatalogueType;

class ParallelPlateFluidModel
{
public:
  /// The ParallelPlateFluidModel Factory.
  static ParallelPlateFluidModelBase* NewParallelPlateFluidModel(const std::string& ParallelPlateFluidModelName,
                                                                 TICPP::HierarchicalDataNode* const hdn);

  /// Interface to the ParallelPlateFluidModel name -> ParallelPlateFluidModel
  // initializer map
  static ParallelPlateFluidModelCatalogueType& GetParallelPlateFluidModelCatalogue();

  /// Return a list of supported ParallelPlateFluidModel names
  static void GetParallelPlateFluidModelNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from
// ParallelPlateFluidModelInitializer
template<class ParallelPlateFluidModelType>
class ParallelPlateFluidModelRegistrator : public ParallelPlateFluidModelInitializer
{
public:
  ParallelPlateFluidModelRegistrator(void)
  {
    const std::string parallelPlateFluidModelName(ParallelPlateFluidModelType::FluidModelName());
    ParallelPlateFluidModelFactory::GetParallelPlateFluidModelCatalogue()[parallelPlateFluidModelName] = this;
  }

  ParallelPlateFluidModelBase* InitializeParallelPlateFluidModel(TICPP::HierarchicalDataNode* const hdn)
  {
    if(!hdn)
      throw GPException("Need to specify a valid HierarchicalDataNode to InitializeParallelPlateFluidModel");

    ParallelPlateFluidModelBase* ret = new ParallelPlateFluidModelType(pm);
    ret->ReadXML(hdn);
    return ret;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_PARALLEL_PLATE_FLUID_MODEL( ClassName ) namespace { ParallelPlateFluidModelRegistrator<ClassName> reg_ppfm_ ## ClassName; }

#endif
