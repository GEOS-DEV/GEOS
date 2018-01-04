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
 * @file FunctionManager.h
 * @author walsh24
 * @date May 24, 2011
 */

#ifndef FUNCTIONMANAGER_H_
#define FUNCTIONMANAGER_H_

//#include "Common/typedefs.h"

#include <map>
#include "codingUtilities/Functions.hpp"

#ifdef USE_ATK
#include <slic/slic.hpp>
#endif

class FunctionManager
{
public:

  static FunctionManager& Instance()
  {
    static FunctionManager theFunctionManager;

    return theFunctionManager;
  }

  inline std::map< std::string, Function* >& Functions()
  { return m_functions; }

  inline Function& GetFunction( const std::string& functionName ) const
  {
    std::map<std::string,Function* >::const_iterator function = m_functions.find( functionName );
    if( function == m_functions.end() )
    {
#ifdef USE_ATK
      SLIC_ERROR("Error FunctionManager: Function name `" + functionName + "' not found\n");
#endif
    }

    return *(function->second);
  }


  inline void AddFunction( const std::string& functionName, Function* aFunction ){
    m_functions.insert(std::make_pair(functionName, aFunction));
  }

  inline realT Eval( const std::string& functionName, const realT& key ) const
  {
    Function& f = GetFunction(functionName);
    return f(key);
  }

private:
  std::map<std::string,Function* > m_functions;

  FunctionManager():
    m_functions()
  {}

  ~FunctionManager() {
    for( std::map<std::string,Function*>::iterator it_func = m_functions.begin() ;
         it_func != m_functions.end() ; ++it_func )
    {
      delete it_func->second;
    }
  }

  FunctionManager( const FunctionManager& );
  FunctionManager& operator=( const FunctionManager& );
};

#endif /* FUNCTIONMANAGER_H_ */
