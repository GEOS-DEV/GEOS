//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
