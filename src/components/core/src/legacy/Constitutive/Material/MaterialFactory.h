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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file MaterialFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef MATERIALFACTORY_H_
#define MATERIALFACTORY_H_

#include "MaterialBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// Material Factory
//
// Consists of the following parts:
//   * The function to generate new material pointers: "NewMaterial"
//   * A base class to derive the functions to generate material pointers: "MaterialInitializer"
//   * A String-to-Material-Intializer map hidden behind the GetMaterialCatalogue function
//   * A template to create material initializers: "MaterialRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_MATERIAL"
// 
// Most materials will only need to use one or two of the parts:
//   * To register a new material in the factory: REGISTER_MATERIAL( MaterialClassName )
//   * To load a material pointer from the factory:       MaterialBase* aMaterialPtr = MaterialFactory::NewMaterial(materialString, args );

/// Base class to generate new Material pointers
class MaterialInitializer
{
public:
  virtual
#if USECPP11==1
  std::unique_ptr<MaterialBase>
#else
  MaterialBase*
#endif
  InitializeMaterial(TICPP::HierarchicalDataNode* hdn) = 0;

  virtual ~MaterialInitializer()
  {
  }
};

typedef std::map<std::string, MaterialInitializer*> MaterialCatalogueType;

class MaterialFactory
{
public:
  /// The Material Factory.
  static
#if USECPP11==1
  std::unique_ptr<MaterialBase>
#else
  MaterialBase*
#endif
 NewMaterial(const std::string& materialName,
                                   TICPP::HierarchicalDataNode* hdn = 0);

  /// Interface to the Material name -> Material initializer map
  static MaterialCatalogueType& GetMaterialCatalogue();

  /// Return a list of supported material names
  static void GetMaterialNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from MaterialInitializer
/// But note that not all materials need to be registered in this way (eg. Geodyn material).
template<class MaterialType>
class MaterialRegistrator: public MaterialInitializer
{
public:
  MaterialRegistrator(void)
  {
    MaterialFactory::GetMaterialCatalogue()[MaterialType::Name()] = this;
  }

#if USECPP11==1
  std::unique_ptr<MaterialBase>
#else
  MaterialBase*
#endif
  InitializeMaterial(TICPP::HierarchicalDataNode* hdn)
  {

#if USECPP11==1
    std::unique_ptr<MaterialBase> tmp = std::unique_ptr<MaterialBase>(new MaterialType());
#else
    MaterialBase* tmp = new MaterialType();
#endif

    if(hdn)
    {
      tmp->resize(0,1);
      tmp->ReadXML(*hdn);
    }
    return tmp;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_MATERIAL( ClassName ) namespace{ MaterialRegistrator<ClassName> reg_##ClassName; }

#endif
