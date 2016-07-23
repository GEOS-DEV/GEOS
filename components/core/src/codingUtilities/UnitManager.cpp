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
 * @file UnitManager.cpp
 * @author walsh24
 * @date July 27, 2011
 */

#include "UnitManager.h"
#include "codingUtilities/StringUtilities.hpp"
#include "fileIO/ticpp/TinyXMLParser.h"
#include <cstring>
#include <algorithm>

namespace geosx
{

/**
 * @author walsh24
 * 
 * Units may be set in terms of predefined default systems of units.
 *   <Units default_units="SI_Units" />      * SI units: m, kg, s, K, mol
 *   <Units default_units="cgs_Units" />     * cm, g, s, K, mol
 * 
 * Units may also be set individually - either in terms of m,g,s,K,mol
 *   <Units length="1.0"            * Default length unit in meters
 *          mass= "1.0"             * Default mass unit in grams
 *          time= "1.0"             * Default time unit in seconds
 *          temperature= "1.0"      * Default temperature unit in Kelvin
 *          mole= "1.0" />          * Default chemical unit in mole
 * 
 * or defined explicitly in terms of known quantities
 *   <Units length="km" mass="mg" time="year" temperature= "K" mole= "1e-3 mol" /> 
 * 
 * The methods may also be combined to override the default units eg. 
 *   <Units default_units="SI_Units" time="day">  * m, kg, day, K, mol  
 * 
 * 
 **/
void UnitManager::ReadXML( TICPP::HierarchicalDataNode* hdn ){

  using namespace GPUnits;
  
  real64 length =0.0;
  real64 mass =0.0;
  real64 time =0.0;
  real64 temp =0.0;
  real64 mole =0.0;
  
  std::string lengthStr = hdn->GetAttributeStringOrDefault("length","");
  std::string massStr = hdn->GetAttributeStringOrDefault("mass","");
  std::string timeStr = hdn->GetAttributeStringOrDefault("time","");
  std::string tempStr = hdn->GetAttributeStringOrDefault("temperature","");
  std::string moleStr = hdn->GetAttributeStringOrDefault("mole","");
  
  // reset base units to 1 so we can parse unit length mass etc. strings
  bool runParser = false;
  if(!lengthStr.empty()){ m_baseUnits[UNIT_LENGTH]= 1; runParser = true;} 
  if(!massStr.empty()){ m_baseUnits[UNIT_MASS]= 1; runParser = true;} 
  if(!timeStr.empty()){ m_baseUnits[UNIT_TIME]= 1; runParser = true;} 
  if(!tempStr.empty()){ m_baseUnits[UNIT_TEMPERATURE]= 1; runParser = true;} 
  if(!moleStr.empty()){ m_baseUnits[UNIT_MOLE]= 1; runParser = true;} 
  
  if(runParser){
    UpdateParser();
    if(!lengthStr.empty()) length = Convert(lengthStr);
    if(!massStr.empty())   mass = Convert(massStr);
    if(!timeStr.empty())   time = Convert(timeStr);
    if(!tempStr.empty())   temp = Convert(tempStr);
    if(!moleStr.empty())   mole = Convert(moleStr);  
  }

  std::string defaultUnitStr = hdn->GetAttributeStringOrDefault("default_units","");
  
  if(!defaultUnitStr.empty()){
  	if( stringutilities::streq(defaultUnitStr,"SI_Units") ){
  	  // m kg s K mol
      m_baseUnits[UNIT_LENGTH]= m;
      m_baseUnits[UNIT_MASS]= kg;
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"SI_Units_mm") ){
      // millimeters, metric tons, seconds, Kelvin, mol
      m_baseUnits[UNIT_LENGTH]= mm;
      m_baseUnits[UNIT_MASS]= 1e3*kg;
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"Meters_Geodyn_Units") ){
        // meters, metric tons, milliseconds, Kelvin, mol
        m_baseUnits[UNIT_LENGTH]= m;
        m_baseUnits[UNIT_MASS]= 1e3*kg;
        m_baseUnits[UNIT_TIME]= ms;
        m_baseUnits[UNIT_TEMPERATURE]= K;
        m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"Geodyn_Units") || stringutilities::streq(defaultUnitStr,"Microscale_Units") ){
      // millimeters, milligrams, microseconds, Kelvin, mol
      m_baseUnits[UNIT_LENGTH]= mm;
      m_baseUnits[UNIT_MASS]= mg;
      m_baseUnits[UNIT_TIME]= us;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"cgs_Units") ){
      // centimeters, grams, seconds, Kelvin, mol
      m_baseUnits[UNIT_LENGTH]= cm;
      m_baseUnits[UNIT_MASS]= g;
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"BDiv_Units") ){
      // centimeters, grams, seconds, Kelvin, mol
      m_baseUnits[UNIT_LENGTH]= cm;
      m_baseUnits[UNIT_MASS]= g;
      m_baseUnits[UNIT_TIME]= us;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"US_Units_ft") ){
      // ft, slug, seconds, Kelvin, mol
      m_baseUnits[UNIT_LENGTH]= ft;
      m_baseUnits[UNIT_MASS]= slug;
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else if( stringutilities::streq(defaultUnitStr,"US_Units_in") ){
      // in, lbf*s^2/in, seconds, Kelvin, mol
      m_baseUnits[UNIT_LENGTH]= in;
      m_baseUnits[UNIT_MASS]= 175.12684*kg; // lbf*s^2/in -> gives lbf for force units
      m_baseUnits[UNIT_TIME]= s;
      m_baseUnits[UNIT_TEMPERATURE]= K;
      m_baseUnits[UNIT_MOLE]= mol;
    } else {
  	  SLIC_ERROR("UnitManager: did not recognize " + defaultUnitStr + " default units");
    }
  }  
  

  // change base units (override default strings)	
  if(length > 0.0) m_baseUnits[UNIT_LENGTH] = length;
  if(mass > 0.0) m_baseUnits[UNIT_MASS] = mass;
  if(time > 0.0) m_baseUnits[UNIT_TIME] = time;
  if(temp > 0.0) m_baseUnits[UNIT_TEMPERATURE] = temp;
  if(mole > 0.0) m_baseUnits[UNIT_MOLE] = mole;
  
  UpdateParser();
  
  
  // User-defined units
  //
  // eg.
  //   <Units length="mm" mass="g" time="s" temperature="K" mole="mol" >
  //     <AddUnit name="ha" value="10000 m^2" />
  //   </Units>
  //
  // User-defined units cannot be used to set base units.
  for (TICPP::HierarchicalDataNode* childNode = hdn->Next(true); childNode; childNode = hdn->Next())
  {

    if (stringutilities::streq(childNode->Heading(), "AddUnit"))
    {
      std::string nameStr = childNode->GetAttributeString("name");
      std::string valueStr = childNode->GetAttributeString("value");
      DefineUnit(nameStr,valueStr);
    }

  }
}

/**
 * @author walsh24
 * @brief Updates unit definitions following changes to base units.
 * 
 * This method is called after the base units have been changed. 
 * Units supported by the parser are defined in this function. 
 * 
 **/
void UnitManager::UpdateParser(void){

  using namespace GPUnits;
  // Length
  DefineUnit("m",Units(1),"meter"); // meters
  DefineUnit("cm",cm*Units(1),"centimeter");
  DefineUnit("mm",mm*Units(1),"millimeter");
  DefineUnit("um",um*Units(1),"micron");  // micrometers
  DefineUnit("nm",nm*Units(1),"nanometer");
  DefineUnit("km",km*Units(1),"kilometer");
  
  DefineUnit("in",in*Units(1),"inch"); // inches
  DefineUnit("ft",ft*Units(1),"foot"); // feet
  DefineUnit("yd",yd*Units(1),"yard"); // yards
  DefineUnit("mi",mi*Units(1),"mile"); // miles

  // Area (Permeability)
  DefineUnit("darcy",darcy*Units(2),"darcy");
  DefineUnit("mD",milli*darcy*Units(2),"millidarcy");
  DefineUnit("uD",micro*darcy*Units(2),"microdarcy");
  

  // Volume
  DefineUnit("L",1e-3*Units(3),"liter"); // liter
  DefineUnit("l",1e-3*Units(3),"liter");
  DefineUnit("ml",1e-6*Units(3),"milliliter");
  DefineUnit("mL",1e-6*Units(3),"milliliter");
  DefineUnit("cc",1e-6*Units(3),"cubic centimeter"); // cubic centimeter
  DefineUnit("gal",3.785411784*1e-3*Units(3),"US liquid gallon"); // US liquid gallon (gal = 3.785 l)
  DefineUnit("bbl",158.987294928*1e-3*Units(3),"oil barrel"); // oil barrel (bbl = 158.987 l)

  // Mass
  DefineUnit("g",Units(0,1),"gram"); // gram
  DefineUnit("mg",mg*Units(0,1),"milligram");
  DefineUnit("kg",kg*Units(0,1),"kilogram");
  
  DefineUnit("lb",lb*Units(0,1),"pound (mass)"); // pound
  DefineUnit("slug",slug*Units(0,1),"slug"); // slug

  // Time
  DefineUnit("s",Units(0,0,1),"second");
  DefineUnit("ms",ms*Units(0,0,1),"millisecond");
  DefineUnit("us",us*Units(0,0,1),"microsecond");  // micro seconds
  DefineUnit("minute",GPUnits::min*Units(0,0,1),"minute");  // "min" not used to avoid conflict with function parser minimum
  DefineUnit("hour",hour*Units(0,0,1),"hour");
  DefineUnit("day",day*Units(0,0,1),"day");
  DefineUnit("year",year*Units(0,0,1),"year");
  
  // Temperature
  DefineUnit("K",Units(0,0,0,1),"Kelvin");
  

  // Amount
  DefineUnit("mol",Units(0,0,0,0,1),"mole");
  DefineUnit("mmol",mmol*Units(0,0,0,0,1),"millimole");

  // Force
  DefineUnit("N",kg*Units(1,1,-2),"Newton"); // Newton = kg.m/s^2
  DefineUnit("dyn",centi*Units(1,1,-2),"dyne"); // dyne = g.cm/s^2
  DefineUnit("lbf",4.44822162*kg*Units(1,1,-2),"pound (force)"); // pound force = 4.44822162*N
  
  // Energy
  DefineUnit("J",kg*Units(2,1,-2),"Joule"); // Joule = kg.m^2/s^2
  DefineUnit("kJ",kilo*kg*Units(2,1,-2),"kilojoule"); // kJ = 1e3*kg.m^2/s^2
  
  DefineUnit("cal",4.184*kg*Units(2,1,-2),"Thermochemical calorie"); // Thermochemical calorie = 4.184 J
  DefineUnit("kcal",4.184*kilo*kg*Units(2,1,-2),"kilocalorie"); //
  
  DefineUnit("eV",1.602176487e-19*kg*Units(2,1,-2),"electron volt"); //Electron volt = 1.602176487e-19 J // nb to use eV in the parser leave a space between the number and the units, eg. "12 eV"
 
  // Work
  DefineUnit("W",kg*Units(2,1,-3),"Watt"); // Watt = kg.m^2/s^3
  DefineUnit("kW",kilo*kg*Units(2,1,-3),"kilowatt"); // kW = 1e3*kg.m^2/s^3
 
  // Pressure
  DefineUnit("Pa",kg*Units(-1,1,-2),"Pascal"); // Pascal = kg/(m.s^2)
  DefineUnit("kPa",kilo*kg*Units(-1,1,-2),"kiloPascal"); // kPa = 10^3 Pa
  DefineUnit("MPa",mega*kg*Units(-1,1,-2),"Megapascal"); // MPa = 10^6 Pa
  DefineUnit("GPa",giga*kg*Units(-1,1,-2),"Gigapascal"); // GPa = 10^9 Pa
  DefineUnit("bar",1e5*kg*Units(-1,1,-2),"bar"); // bar = 10^5 Pa
  DefineUnit("atm",1.01325e5*kg*Units(-1,1,-2),"atmosphere"); // atm = 1.01325 x10^5 Pa
  DefineUnit("Torr",133.322*kg*Units(-1,1,-2),"Torr"); // Torr = 133.322 Pa
  DefineUnit("psi",6.895e3*kg*Units(-1,1,-2),"pounds per square inch"); // psi = 6.895 x10^3 Pa
  
  // Density
  DefineUnit("ppg",119.826427*kg*Units(-3,1),"pounds per US liquid gallon"); // pounds per US liquid gallon = 119.826427 kg / m3

  // Dynamic Viscosity
  DefineUnit("cP",1e-3*kg*Units(-1,1,-1),"centipoise"); // centipoise = 1e-3 Pa.s


}

/**
 * @author walsh24
 * @brief Add a user defined unit to the unit manager.
 * 
 **/
void UnitManager::DefineUnit(const std::string& symbolString,const std::string& valueStr, const std::string& name){
  real64 value = Convert(valueStr);
  //m_fParser.AddConstant(symbolString.c_str(),value);
  this->DefineUnit(symbolString, value, name);
}

void UnitManager::DefineUnit(const std::string& symbolString, const real64& value, const std::string& name){
  GPUnits::UnitData udStruct;
  udStruct.name = name;
  udStruct.value = value;
  m_unitValueMap[symbolString] = udStruct;
//  m_fParser.AddConstant(symbolString.c_str(),value);
}

/**
 * @author walsh24
 * @brief Converts a string to the preset default units.
 * 
 * The following are valid (and equivalent) strings:
 * "9.8m.s^-2"
 * "0.98e1m/s^2"
 * "9.8m/(s.s)"
 * "9.8 m/(s.s)"
 * "0.0098 km/(s*s)"
 * "9.8" (Assuming the default units are m and s)
 *
 *
 * NB. if the unit symbol begins with "e" or "E"
 * a space must be left between the number and the 
 * units to parse the string properly eg. "1.23 eV"
 * otherwise the e will be interpreted as an exponent. 
 **/
real64 UnitManager::Convert(const std::string& unitString){

  // find first non-numeric
  size_t  indx = unitString.find_first_not_of(" \t"); // skip leading whitespace
  indx = unitString.find_first_not_of("0123456789.eE-+",indx);  // nb cannot include " " - needed to parse units starting with "e" or "E"
  real64 number = (indx > 0)? atof(unitString.substr(0,indx).c_str()) :  1.0;

  if(indx == std::string::npos) return number; // dimensionless or in default units

  return Convert(number, unitString.substr(indx));
}

/**
 * @author walsh24
 * @brief Converts a quantity and units to the preset default units.
 * 
 * eg.
 *  UnitManager& theUnitManager = UnitManager::Instance();
 *  real64 g = theUnitManager.Convert(9.8,"m/s^2");
 * 
 **/
//real64 UnitManager::Convert(real64 number, std::string unitString){
real64 UnitManager::Convert(real64 , std::string ){
#pragma message ("reenable this function")
//
//
//
//  // replace "." with "*"
//  std::replace(unitString.begin(), unitString.end(), '.', '*');
//
//  int err = m_fParser.Parse(unitString.c_str(), "");
//  if(err >=0){
//  	throw GPException("Error detected in UnitManager expression: '"+ unitString+ "' at character " + toString(err+1) + "." );
//  }
//
//  const double Dummy=0;
//  real64 scale =  m_fParser.Eval(&Dummy);
//
//  return number*scale;
  return 0.0;
}
/**
 * @author walsh24
 * @brief Converts a quantity to the preset default units.
 * 
 * This method is faster than the string parser.
 * 
 * eg.
 *  UnitManager& theUnitManager = UnitManager::Instance();
 *  using namespace GPUnits; // simplifies use of GPUnits::m, GPUnits::s etc. 
 *  real64 g = theUnitManager.Convert(9.8,m,1,s,-2);
 * 
 **/
real64 UnitManager::Convert( real64 /*value*/ ,
                            real64 lengthUnit,
                            int lengthdim,
                            real64 massUnit,
                            int massdim,
                            real64 timeUnit,
                            int timedim,
                            real64 tempUnit,
                            int tempdim,
                            real64 moleUnit,
                            int moledim)
{
   real64 rv = Units(lengthdim,massdim,timedim,tempdim,moledim);
   if(lengthdim != 0) rv *= pow(lengthUnit,lengthdim);   
   if(massdim != 0) rv *= pow(massUnit,massdim);   
   if(timedim != 0) rv *= pow(timeUnit,timedim);   
   if(tempdim != 0) rv *= pow(tempUnit,tempdim);   
   if(moledim != 0) rv *= pow(moleUnit,moledim);
   
   return rv;                         	
}
                               
}


    
    
