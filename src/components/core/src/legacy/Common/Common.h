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
 * @file ProblemManagerT.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef _COMMON_H_
#define _COMMON_H_



#include <vector>
#include <algorithm>
#include <map>


#include <memory>

#if __ICC
  #include <cmath>
#else
  #include <cmath>
//  #include <omp.h>
#endif


//#include "../ArrayT/bufvector.h"
#include "typedefs.h"
#include "GPException.h"
#include "GlobalIndexManager.h"

const int one = 1;

#define USECPP11 0

#if USECPP11==1
#else
#endif

#if USECPP11==1
#include <array>
#else
#include <stdio.h>
#endif


/**
 * @author Randolph Settgast
 * @brief base class for the field structures
 * This class serves as the base class for all the field structures which hold attribute information
 * for a field.
 *
 */
class FieldBase
{
public:
  /**
   * @author Randolph Settgast
   * @param key the key of the field
   * @param name the name of the field
   * @param WriteToRestart flag to see if the field gets written to restart files
   * @param WriteToPlot flag to see if the field gets written to plot files
   */
  FieldBase( const size_t key, const std::string& name, const bool WriteToRestart, const bool WriteToPlot )
  :m_key(key), m_name(name), m_WriteToRestart(WriteToRestart), m_WriteToPlot(WriteToPlot)
  {}

  virtual ~FieldBase(){}

public:
  /// the field key
  size_t m_key;

  /// the field name
  std::string m_name;

  /// flag to see if the field gets written to restart files
  bool m_WriteToRestart;

  /// flag to see if the field gets written to plot files
  bool m_WriteToPlot;

private:
  /// non-callable default constructor
  FieldBase();

  /// non-callable copy constructor
  FieldBase(const FieldBase&);

};

namespace FieldInfo
{
  enum FieldEnum
  {
    noKey,
    isDomainBoundary,
    isExternal,
    processColor,
    ownedByRank,
    ghostRank,
    referencePosition,
    currentPosition,
    displacement,
    incrementalDisplacement,
    relativePosition,
    velocity,
    acceleration,
    force,
    contactForce,
    hgforce,
    mass,
    rotationAxis,
    rotationMagnitude,
    rotationalAxisIncrement,
    rotationalMagnitudeIncrement,
    rotationalVelocity,
    rotationalAcceleration,
    moment,
    rotationalInertia,
    pressure,
    fluidPressure,
    deviatorStress,
    totalStress,
    volume,
    density,
    fluidDensity,
    massFlux,
    demIndex,
    seismicBeginTime,
    seismicRiseTime,
    slip,
    area,
    seismicMagnitude,
    seismicMoment,
    damageIndicator,
    numFieldEnums
  };

  enum FieldType
  {
    integerField = 0,
    localIndexField = 1,
    globalIndexField = 2,
    realField = 3,
    R1TensorField = 4,
    R2TensorField = 5,
    R2SymTensorField = 6,
    integerParameter = 7,
    realParameter = 8,
    numFieldTypes = 9
  };

  extern std::vector<FieldBase*> AttributesByKey;
  extern std::map<std::string, FieldBase*> AttributesByName;

  inline size_t FieldSize(FieldType ft)
  {
    /*size_t rval;
     if( i == integerField )
     rval = 1;
     else if( i == realField )
     rval = 1;
     else if( i == R1TensorField )
     rval = 3;
     else if( i == R2TensorField )
     rval = 9;
     else if( i == R2SymTensorField )
     rval = 6;
     else
     throw GPException("Common.h::FieldSize() invalid\n");*/
    switch (ft)
    {
      case integerField:
        return 1;
        break;
      case localIndexField:
        return 1;
        break;
      case globalIndexField:
        return 1;
        break;
      case realField:
        return 1;
        break;
      case R1TensorField:
        return 3;
        break;
      case R2TensorField:
        return 9;
        break;
      case R2SymTensorField:
        return 6;
        break;
      case integerParameter:
        return 1;
        break;
      case realParameter:
        return 1;
        break;
      case numFieldTypes:
      default:
        throw GPException("Common.h::FieldSize() invalid\n");
    }
    /*Should not get here*/
    throw GPException("Common.h::FieldSize() invalid\n");
    return 0;
  }

  // Type strings
  static const std::string IntegerStr = "Integer";
  static const std::string LocalIndexStr = "LocalIndex";
  static const std::string GlobalIndexStr = "GlobalIndex";
  static const std::string RealStr = "Scalar";
  static const std::string R1TensorStr = "Vector";
  static const std::string R2TensorStr = "Tensor";
  static const std::string R2SymTensorStr = "SymmetricTensor";
  static const std::string IntegerParameterStr = "IntegerParameter";
  static const std::string RealParameterStr = "ScalarParameter";

  inline const std::string& FieldTypeName(FieldType ft)
  {
    switch (ft)
    {
      case integerField:
        return IntegerStr;
        break;
      case localIndexField:
        return LocalIndexStr;
        break;
      case globalIndexField:
        return GlobalIndexStr;
        break;
      case realField:
        return RealStr;
        break;
      case R1TensorField:
        return R1TensorStr;
        break;
      case R2TensorField:
        return R2TensorStr;
        break;
      case R2SymTensorField:
        return R2SymTensorStr;
        break;
      case integerParameter:
        return IntegerParameterStr;
        break;
      case realParameter:
        return RealParameterStr;
        break;
      case numFieldTypes:
      default:
        throw GPException("FieldTypeName: Unrecognized field type");
    }
    /*Should not get here*/
    throw GPException("FieldTypeName: Unrecognized field type");
    return RealStr;
  }

  void AllocateAttributes();
  void DeleteAttributes();

} // Field Info

typedef unsigned int FieldKey;
typedef FieldInfo::FieldType FieldType;




template <FieldKey FIELDKEY>
struct FieldStructWrapper
{};




//*************************************************************************************************
struct isDomainBoundary : public FieldBase
{
  typedef int Type;
  static size_t key() { return FieldInfo::isDomainBoundary; }
  static const char* Name() { return "isDomainBoundary"; }

  isDomainBoundary():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::isDomainBoundary>
{  typedef isDomainBoundary FieldStruct; };


//*************************************************************************************************
struct isExternal : public FieldBase
{
  typedef int Type;
  static size_t key() { return FieldInfo::isExternal; }
  static const char* Name() { return "isExternal"; }

  isExternal():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::isExternal>
{  typedef isExternal FieldStruct; };

//*************************************************************************************************
struct processColor : public FieldBase
{
  typedef int Type;
  static size_t key() { return FieldInfo::processColor; }
  static const char* Name() { return "processGroup"; }

  processColor():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::processColor>
{  typedef processColor FieldStruct; };


//*************************************************************************************************
struct ghostRank : public FieldBase
{
  typedef int Type;
  static size_t key() { return FieldInfo::ghostRank; }
  static const char* Name() { return "ghostRank"; }

  ghostRank():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::ghostRank>
{  typedef ghostRank FieldStruct; };


//*************************************************************************************************
struct ownership : public FieldBase
{
  typedef int Type;
  static size_t key() { return FieldInfo::ownedByRank; }
  static const char* Name() { return "ownership"; }

  ownership():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::ownedByRank>
{  typedef ownership FieldStruct; };


//*************************************************************************************************
struct demIndex : public FieldBase
{
  typedef localIndex Type;
  static size_t key() { return FieldInfo::demIndex; }
  static const char* Name() { return "demIndex"; }

  demIndex():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::demIndex>
{  typedef demIndex FieldStruct; };











//*************************************************************************************************
struct referencePosition : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::referencePosition; }
  static const char* Name() { return "ReferencePosition"; }

  referencePosition():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::referencePosition>
{  typedef referencePosition FieldStruct; };



//*************************************************************************************************
struct currentPosition : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::currentPosition; }
  static const char* Name() { return "CurrentPosition"; }

  currentPosition():
    FieldBase( key(), Name(), true, true )
  {}

};
template<> struct FieldStructWrapper<FieldInfo::currentPosition>
{  typedef currentPosition FieldStruct; };



//*************************************************************************************************
struct displacement : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::displacement; }
  static const char* Name() { return "Displacement"; }

  displacement():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::displacement>
{  typedef displacement FieldStruct; };

//*************************************************************************************************
struct incrementalDisplacement : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::incrementalDisplacement; }
  static const char* Name() { return "IncrementalDisplacement"; }

  incrementalDisplacement():
    FieldBase( key(), Name(), true, false )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::incrementalDisplacement>
{  typedef incrementalDisplacement FieldStruct; };

//*************************************************************************************************
struct relativePosition : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::relativePosition; }
  static const char* Name() { return "RelativePosition"; }

  relativePosition():
    FieldBase( key(), Name(), true, false )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::relativePosition>
{  typedef relativePosition FieldStruct; };

//*************************************************************************************************
struct velocity : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::velocity; }
  static const char* Name() { return "Velocity"; }
  velocity():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::velocity>
{  typedef velocity FieldStruct; };

//*************************************************************************************************
struct acceleration : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::acceleration; }
  static const char* Name() { return "Acceleration"; }
  acceleration():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::acceleration>
{  typedef acceleration FieldStruct; };

//*************************************************************************************************
struct force : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::force; }
  static const char* Name() { return "Force"; }
  force():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::force>
{  typedef force FieldStruct; };

//*************************************************************************************************
struct hgforce : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::hgforce; }
  static const char* Name() { return "hgforce"; }
  hgforce():
    FieldBase( key(), Name(), false, false )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::hgforce>
{  typedef hgforce FieldStruct; };

//*************************************************************************************************
struct contactForce : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::contactForce; }
  static const char* Name() { return "contactForce"; }
  contactForce():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::contactForce>
{  typedef contactForce FieldStruct; };







//*************************************************************************************************
struct rotationAxis : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::rotationAxis; }
  static const char* Name() { return "RotationAxis"; }

  rotationAxis():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationAxis>
{  typedef rotationAxis FieldStruct; };


//*************************************************************************************************
struct rotationMagnitude : public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::rotationMagnitude; }
  static const char* Name() { return "RotationMagnitude"; }

  rotationMagnitude():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationMagnitude>
{  typedef rotationMagnitude FieldStruct; };

//*************************************************************************************************
struct rotationalAxisIncrement : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::rotationalAxisIncrement; }
  static const char* Name() { return "RotationAxisIncrement"; }

  rotationalAxisIncrement():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationalAxisIncrement>
{  typedef rotationalAxisIncrement FieldStruct; };


//*************************************************************************************************
struct rotationalMagnitudeIncrement : public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::rotationalMagnitudeIncrement; }
  static const char* Name() { return "RotationIncrement"; }

  rotationalMagnitudeIncrement():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationalMagnitudeIncrement>
{  typedef rotationalMagnitudeIncrement FieldStruct; };




//*************************************************************************************************
struct rotationalVelocity : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::rotationalVelocity; }
  static const char* Name() { return "RotationalVelocity"; }

  rotationalVelocity():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationalVelocity>
{  typedef rotationalVelocity FieldStruct; };


//*************************************************************************************************
struct rotationalAcceleration : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::rotationalAcceleration; }
  static const char* Name() { return "RotationalAcceleration"; }

  rotationalAcceleration():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationalAcceleration>
{  typedef rotationalAcceleration FieldStruct; };


//*************************************************************************************************
struct moment : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::moment; }
  static const char* Name() { return "Moment"; }

  moment():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::moment>
{  typedef moment FieldStruct; };

//*************************************************************************************************
struct rotationalInertia : public FieldBase
{
  typedef R1Tensor Type;
  static size_t key() { return FieldInfo::rotationalInertia; }
  static const char* Name() { return "RotationalInertia"; }

  rotationalInertia():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::rotationalInertia>
{  typedef rotationalInertia FieldStruct; };


//*************************************************************************************************
struct volume : public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::volume; }
  static const char* Name() { return "Volume"; }
  volume():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::volume>
{  typedef volume FieldStruct; };


//*************************************************************************************************
struct mass : public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::mass; }
  static const char* Name() { return "Mass"; }
  mass():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::mass>
{  typedef mass FieldStruct; };


//*************************************************************************************************
struct Pressure: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::pressure; }
  static const char* Name() { return "Pressure"; }
  Pressure():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::pressure>
{  typedef Pressure FieldStruct; };


//*************************************************************************************************
struct FluidPressure: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::fluidPressure; }
  static const char* Name() { return "FluidPressure"; }
  FluidPressure():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::fluidPressure>
{  typedef FluidPressure FieldStruct; };


//*************************************************************************************************
struct Density: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::density; }
  static const char* Name() { return "Density"; }
  Density():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::density>
{  typedef Density FieldStruct; };


//*************************************************************************************************
struct FluidDensity: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::fluidDensity; }
  static const char* Name() { return "FluidDensity"; }
  FluidDensity():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::fluidDensity>
{  typedef FluidDensity FieldStruct; };



//*************************************************************************************************
struct massFlux: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::massFlux; }
  static const char* Name() { return "massFlux"; }
  massFlux():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::massFlux>
{  typedef massFlux FieldStruct; };


//*************************************************************************************************
struct damageIndicator: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::damageIndicator; }
  static const char* Name() { return "DamageIndicator_ip00"; }
  damageIndicator():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::damageIndicator>
{  typedef damageIndicator FieldStruct; };


//*************************************************************************************************
struct DeviatorStress: public FieldBase
{
  typedef R2SymTensor Type;
  static size_t key() { return FieldInfo::deviatorStress; }
  static const char* Name() { return "DeviatorStress"; }
  DeviatorStress():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::deviatorStress>
{  typedef DeviatorStress FieldStruct; };


//*************************************************************************************************
struct SeismicBeginTime: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::seismicBeginTime; }
  static const char* Name() { return "SeismicBeginTime"; }
  SeismicBeginTime():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::seismicBeginTime>
{  typedef SeismicBeginTime FieldStruct; };


//*************************************************************************************************
struct SeismicRiseTime: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::seismicRiseTime; }
  static const char* Name() { return "SeismicRiseTime"; }
  SeismicRiseTime():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::seismicRiseTime>
{  typedef SeismicRiseTime FieldStruct; };


//*************************************************************************************************
struct Slip: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::slip; }
  static const char* Name() { return "Slip"; }
  Slip():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::slip>
{  typedef Slip FieldStruct; };


//*************************************************************************************************
struct Area: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::area; }
  static const char* Name() { return "Area"; }
  Area():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::area>
{  typedef Area FieldStruct; };


//*************************************************************************************************
struct SeismicMagnitude: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::seismicMagnitude; }
  static const char* Name() { return "SeismicMagnitude"; }
  SeismicMagnitude():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::seismicMagnitude>
{  typedef SeismicMagnitude FieldStruct; };


//*************************************************************************************************
struct SeismicMoment: public FieldBase
{
  typedef realT Type;
  static size_t key() { return FieldInfo::seismicMoment; }
  static const char* Name() { return "SeismicMoment"; }
  SeismicMoment():
    FieldBase( key(), Name(), true, true )
  {}
};
template<> struct FieldStructWrapper<FieldInfo::seismicMoment>
{  typedef SeismicMoment FieldStruct; };



//*************************************************************************************************
template <FieldKey FIELDKEY>
struct Field
{
  typedef typename FieldStructWrapper<FIELDKEY>::FieldStruct::Type Type;

  static const char* Name() { return FieldStructWrapper<FIELDKEY>::FieldStruct::Name(); }
  static size_t key() { return FieldStructWrapper<FIELDKEY>::FieldStruct::key(); }
};

template< typename TYPE >
inline realT& FieldComponent(  TYPE& field, const int i=0 )
{
  return field[i];
}

template<> inline realT& FieldComponent<realT>( realT& field, const int ) {return field;}


#endif /* _CONSTANTS_H_ */
