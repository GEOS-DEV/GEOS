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
 * @file CohesiveZoneFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef COHESIVEZONEFACTORY_H_
#define COHESIVEZONEFACTORY_H_

#include "CohesiveZoneBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// CohesiveZone Factory
//
// Consists of the following parts:
//   * The function to generate new cohesive zone pointers: "NewCohesiveZone"
//   * A base class to derive the functions to generate cohesive zone pointers:
// "CohesiveZoneInitializer"
//   * A String-to-CohesiveZone-Intializer map hidden behind the
// GetCohesiveZoneCatalogue function
//   * A template to create cohesive zone initializers:
// "CohesiveZoneRegistrator"
//   * A compiler directive to simplify autoregistration:
// "REGISTER_COHESIVEZONE"
//
// Most cohesive zones will only need to use one or two of the parts:
//   * To register a new cohesive zone in the factory: REGISTER_COHESIVEZONE(
// CohesiveZoneClassName )
//   * To load a cohesive zone pointer from the factory:       CohesiveZoneBase*
// aCohesiveZonePtr = CohesiveZoneFactory::NewCohesiveZone(cohesiveZoneString,
// args );

/// Base class to generate new CohesiveZone pointers
class CohesiveZoneInitializer
{
public:
  virtual
#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase>
#else
  CohesiveZoneBase*
#endif
  InitializeCohesiveZone(TICPP::HierarchicalDataNode* hdn) = 0;

  virtual ~CohesiveZoneInitializer()
  {}
};

typedef std::map<std::string, CohesiveZoneInitializer*> CohesiveZoneCatalogueType;

class CohesiveZoneFactory
{
public:
  /// The CohesiveZone Factory.
  static
#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase>
#else
  CohesiveZoneBase*
#endif
  NewCohesiveZone(const std::string& cohesiveZoneName,
                  TICPP::HierarchicalDataNode* hdn = 0);

  /// CohesiveZone to the CohesiveZone name -> CohesiveZone initializer map
  static CohesiveZoneCatalogueType& GetCohesiveZoneCatalogue();

  /// Return a list of supported cohesive zone names
  static void GetCohesiveZoneNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from CohesiveZoneInitializer
template<class CohesiveZoneType>
class CohesiveZoneRegistrator : public CohesiveZoneInitializer
{
public:
  CohesiveZoneRegistrator(void)
  {
    CohesiveZoneFactory::GetCohesiveZoneCatalogue()[CohesiveZoneType::Name()] = this;
  }

#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase>
#else
  CohesiveZoneBase*
#endif
  InitializeCohesiveZone(TICPP::HierarchicalDataNode* hdn)
  {

#if USECPP11==1
    std::unique_ptr<CohesiveZoneBase> tmp(new CohesiveZoneType());
#else
    CohesiveZoneBase* tmp = new CohesiveZoneType();
#endif

    if (hdn)
    {
      tmp->resize(0, 1);
      tmp->ReadXML(*hdn);
    }
    return tmp;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_COHESIVEZONE( ClassName ) namespace { CohesiveZoneRegistrator<ClassName> reg_ ## ClassName; }

#endif
