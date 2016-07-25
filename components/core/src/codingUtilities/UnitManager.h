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
 * @file UnitManager.h
 * @author walsh24
 * @date July 27, 2011
 */

#ifndef UNITMANAGER_H_
#define UNITMANAGER_H_

#include "common/DataTypes.hpp"
//#include "Utilities/Utilities.h"
//#include "Utilities/Functions.h"
#include "common/Utilities.hpp"
#include <string>
#include <cmath>
#include <map>

namespace geosx
{
namespace TICPP { class HierarchicalDataNode; }

/// Units in the GPUnits namespace are expressed in m,g,s,K,mol
/// For unit conversion use the UnitManager class. 
namespace GPUnits
{
  enum basic_dimensions {
     UNIT_LENGTH, UNIT_MASS, UNIT_TIME, UNIT_TEMPERATURE, UNIT_MOLE, NUM_UNITS
  };
  
  // SI prefixes
  const real64 giga = 1e9;
  const real64 mega = 1e6;
  const real64 kilo = 1e3;
  const real64 deca = 10.0;
  const real64 deci = 0.1;
  const real64 centi = 1e-2;
  const real64 milli = 1e-3;
  const real64 micro = 1e-6;
  const real64 nano = 1e-9;

  // length units
  const real64 m = 1.0;
  const real64 cm = centi*m;
  const real64 mm = milli*m;
  const real64 km = kilo*m;
  const real64 um = micro*m; // micrometers
  const real64 nm = nano*m;
  const real64 in = 2.54e-2*m;
  const real64 ft = 0.3048*m; //12.0*in
  const real64 yd = 36.0*in;
  const real64 mi = 1760.0*yd;

  // area/permeability units
  const real64 darcy = 9.869233e-13*m*m;
  
  
  // mass units
  const real64 g = 1.0;
  const real64 mg = milli*g;
  const real64 kg = kilo*g;
  const real64 lb = 453.59237*g;
  const real64 slug = 14.593903*kg; // unit of mass used in lbf calculations
  
  // time units
  const real64 s = 1.0;
  const real64 ms = milli*s;
  const real64 us = micro*s; // microseconds
  const real64 min = 60*s;
  const real64 hour = 60*min;
  const real64 day = 24*hour;
  const real64 year = 365.25*day;// Julian year
  const real64 shake = 1e-8*s; // shake = 10 nano seconds;
  
  // temperature units
  const real64 K = 1.0;
  
  // amount units
  const real64 mol = 1.0;
  const real64 mmol = milli*mol;

  struct UnitData{
	  real64 value;
	  std::string name;
  };
 
}

/**
 * @author walsh24
 * 
 * @brief Singleton class to control unit conversion and set default units. 
 * 
 **/
class UnitManager{
	
  public:
    static UnitManager& Instance()
    {
      static UnitManager theUnitManager;
      return theUnitManager;
    };
    void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;
    
    /// Convert string expression to base units
    real64 Convert(const std::string& quantityString);
    /// Convert quantity and string units to base units
    real64 Convert(real64 value, std::string unitString);
    /// Convert quantity with explicitly defined units to base units (faster than string methods)
    real64 Convert(real64 value, real64 lengthUnit, int lengthdim,
                               real64 massUnit= 1.0, int massdim = 0,
                               real64 timeUnit= 1.0, int timedim = 0,
                               real64 tempUnit= 1.0, int tempdim = 0,
                               real64 moleUnit= 1.0, int moledim = 0);
    
    
    /// Convert quantity from base units to units given in a string
    real64 ConvertTo(const std::string& unitString, real64 value);
    /// Convert quantity from base units to explicitly defined units  (faster than string methods)
    real64 ConvertTo(real64 lengthUnit, int lengthdim,
                    real64 massUnit, int massdim,
                    real64 timeUnit, int timedim,
                    real64 tempUnit, int tempdim,
                    real64 moleUnit, int moledim, real64 value);
    // overloaded functions to allow ConvertTo to mirror the form of Convert
    real64 ConvertTo(real64 lengthUnit, int lengthdim, real64 massUnit, int massdim, real64 timeUnit, int timedim,
                    real64 tempUnit, int tempdim, real64 value)
    { 
      return ConvertTo(lengthUnit, lengthdim, massUnit, massdim, timeUnit, timedim,tempUnit, tempdim, 1.0, 0, value);
    };
    real64 ConvertTo(real64 lengthUnit, int lengthdim, real64 massUnit, int massdim, real64 timeUnit, int timedim, real64 value)
    { 
      return ConvertTo(lengthUnit, lengthdim, massUnit, massdim, timeUnit, timedim,1.0, 0, 1.0, 0, value);
    };
    real64 ConvertTo(real64 lengthUnit, int lengthdim, real64 massUnit, int massdim, real64 value)
    { 
      return ConvertTo(lengthUnit, lengthdim, massUnit, massdim, 1.0, 0,1.0, 0, 1.0, 0, value);
    };
    real64 ConvertTo(real64 lengthUnit, int lengthdim, real64 value)
    { 
      return ConvertTo(lengthUnit, lengthdim, 1.0, 0, 1.0, 0,1.0, 0, 1.0, 0, value);
    };


    void ReportAllUnits(void){
      std::map<std::string,GPUnits::UnitData>::iterator itr =  m_unitValueMap.begin();
	  std::map<std::string,GPUnits::UnitData>::iterator iend = m_unitValueMap.end() ;
      std::cout << std::endl;
      std::cout << "Symbol" << "\t \t" << "Description"<< "\n";
      std::cout << "------" << "\t \t" << "-----------"<< "\n";
      for(;itr != iend; ++itr){
        std::cout << itr->first << "\t \t" << itr->second.name << "\n";
      }
      std::cout << std::endl;
    };
 
  private:


 
    /**
     *  @brief Constructor
     * 
     *  NB. The units are set to SI units by default, 
     *  Within the class, the base units are expressed in m,g,s,K,mol.
     * 
     **/
    UnitManager()//:
//      m_fParser()
    {
      using namespace GPUnits;
    	// SI units - m kg s K mol
      m_baseUnits[UNIT_LENGTH]= m;
      m_baseUnits[UNIT_MASS]= kg;
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
      UpdateParser();
    };

    ~UnitManager() {};
    UnitManager( const UnitManager& );
    UnitManager& operator=( const UnitManager& );

    real64 Units( real64 length= 0.0,
                 real64 mass= 0.0,
                 real64 time= 0.0,
                 real64 temp  = 0.0,
                 real64 mole =0.0)
    {
      (void)temp;
       real64 rv = 1.0;
       if( !asctoolkit::utilities::isNearlyEqual( length, 0.0 ) ) rv /= pow(m_baseUnits[0],length);
       if( !asctoolkit::utilities::isNearlyEqual( mass, 0.0 ) )   rv /= pow(m_baseUnits[1],mass);
       if( !asctoolkit::utilities::isNearlyEqual( time, 0.0 ) )   rv /= pow(m_baseUnits[2],time);
       if( !asctoolkit::utilities::isNearlyEqual( mole, 0.0 ) )   rv /= pow(m_baseUnits[3],mole);
       return rv;
    };
    void UpdateParser(void);
    
    /// Add a new unit symbol to the unit manager
    void DefineUnit(const std::string& name,const std::string& value,const std::string& description="");
    void DefineUnit(const std::string& name,const real64& value,const std::string& description="");


    real64 m_baseUnits[GPUnits::NUM_UNITS]; // problem default units given in m,g,s,K,mol.
//    FunctionParser m_fParser;
    std::map<std::string,GPUnits::UnitData> m_unitValueMap;
};

//////////////////////////////////////////////////////////////////////////////////////////

/**
 * @author walsh24
 * @brief Converts a quantity in the default units to another set of units.
 * 
 * eg.
 *  UnitManager& theUnitManager = UnitManager::Instance();
 *  real64 gbase = theUnitManager.Convert("9.8 m/s^2");  // in base units
 *  real64 g = theUnitManager.ConvertTo("km/s^2",gbase);  // g = 0.0098
 * 
 **/
 inline
real64 UnitManager::ConvertTo(const std::string& unitString,real64 value){
	real64 invVal = Convert(1.0, unitString);
	return value/invVal;
}

/**
 * @author walsh24
 * @brief Converts a quantity in the default units to another set of units.
 * 
 * eg.
 *  UnitManager& theUnitManager = UnitManager::Instance();
 *  real64 gbase = theUnitManager.Convert("9.8 m/s^2");  // in base units
 *  using namespace GPUnits;
 *  real64 g = theUnitManager.ConvertTo(km,1,s,-2,gbase);  // g = 0.0098
 * 
 **/
 inline
real64 UnitManager::ConvertTo(real64 lengthUnit, int lengthdim,
                             real64 massUnit, int massdim, real64 timeUnit, int timedim,
                             real64 tempUnit, int tempdim, real64 moleUnit, int moledim, real64 value){
   
   real64 invVal = Convert(1.0, lengthUnit, lengthdim, massUnit, massdim, timeUnit, timedim,
                              tempUnit, tempdim, moleUnit, moledim);
   return value/invVal;                     	
}


}


#endif /* UNITMANAGER_H_ */
