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

#include "Common.h"



std::vector<FieldBase*> FieldInfo::AttributesByKey(numFieldEnums);
std::map<std::string, FieldBase*> FieldInfo::AttributesByName;

// runtime callable functions
void FieldInfo::AllocateAttributes()
{
  AttributesByKey[isDomainBoundary]        = new FieldStructWrapper<isDomainBoundary>::FieldStruct;
  AttributesByKey[isExternal]              = new FieldStructWrapper<isExternal>::FieldStruct;
  AttributesByKey[processColor]              = new FieldStructWrapper<processColor>::FieldStruct;
  AttributesByKey[ghostRank]                 = new FieldStructWrapper<ghostRank>::FieldStruct;
  AttributesByKey[ownedByRank]                 = new FieldStructWrapper<ownedByRank>::FieldStruct;
  AttributesByKey[demIndex]                 = new FieldStructWrapper<demIndex>::FieldStruct;

  AttributesByKey[referencePosition]        = new FieldStructWrapper<referencePosition>::FieldStruct;
  AttributesByKey[currentPosition]          = new FieldStructWrapper<currentPosition>::FieldStruct;
  AttributesByKey[displacement]             = new FieldStructWrapper<displacement>::FieldStruct;
  AttributesByKey[incrementalDisplacement]  = new FieldStructWrapper<incrementalDisplacement>::FieldStruct;
  AttributesByKey[relativePosition]         = new FieldStructWrapper<relativePosition>::FieldStruct;
  AttributesByKey[velocity]                 = new FieldStructWrapper<velocity>::FieldStruct;
  AttributesByKey[acceleration]             = new FieldStructWrapper<acceleration>::FieldStruct;
  AttributesByKey[force]                    = new FieldStructWrapper<force>::FieldStruct;
  AttributesByKey[hgforce]                  = new FieldStructWrapper<hgforce>::FieldStruct;
  AttributesByKey[contactForce]             = new FieldStructWrapper<contactForce>::FieldStruct;
  AttributesByKey[mass]                     = new FieldStructWrapper<mass>::FieldStruct;
  AttributesByKey[rotationAxis]             = new FieldStructWrapper<rotationAxis>::FieldStruct;
  AttributesByKey[rotationMagnitude]             = new FieldStructWrapper<rotationMagnitude>::FieldStruct;
  AttributesByKey[rotationalAxisIncrement]       = new FieldStructWrapper<rotationalAxisIncrement>::FieldStruct;
  AttributesByKey[rotationalMagnitudeIncrement]    = new FieldStructWrapper<rotationalMagnitudeIncrement>::FieldStruct;
  AttributesByKey[rotationalVelocity]             = new FieldStructWrapper<rotationalVelocity>::FieldStruct;
  AttributesByKey[rotationalAcceleration]             = new FieldStructWrapper<rotationalAcceleration>::FieldStruct;
  AttributesByKey[moment]                        = new FieldStructWrapper<moment>::FieldStruct;
  AttributesByKey[rotationalInertia]             = new FieldStructWrapper<rotationalInertia>::FieldStruct;

  AttributesByKey[pressure]                 = new FieldStructWrapper<pressure>::FieldStruct;
  AttributesByKey[fluidPressure]                 = new FieldStructWrapper<fluidPressure>::FieldStruct;
  AttributesByKey[deviatorStress]           = new FieldStructWrapper<deviatorStress>::FieldStruct;
  AttributesByKey[volume]                   = new FieldStructWrapper<volume>::FieldStruct;
  AttributesByKey[density]                  = new FieldStructWrapper<density>::FieldStruct;
  AttributesByKey[fluidDensity]                  = new FieldStructWrapper<fluidDensity>::FieldStruct;
  AttributesByKey[massFlux]                 = new FieldStructWrapper<massFlux>::FieldStruct;
  AttributesByKey[damageIndicator]          = new FieldStructWrapper<damageIndicator>::FieldStruct;

  AttributesByKey[seismicBeginTime]         = new FieldStructWrapper<seismicBeginTime>::FieldStruct;
  AttributesByKey[seismicRiseTime]          = new FieldStructWrapper<seismicRiseTime>::FieldStruct;
  AttributesByKey[slip]                     = new FieldStructWrapper<slip>::FieldStruct;
  AttributesByKey[area]                     = new FieldStructWrapper<area>::FieldStruct;
  AttributesByKey[seismicMagnitude]         = new FieldStructWrapper<seismicMagnitude>::FieldStruct;
  AttributesByKey[seismicMoment]            = new FieldStructWrapper<seismicMoment>::FieldStruct;


  AttributesByName[Field<isDomainBoundary>::Name()]        = new FieldStructWrapper<isDomainBoundary>::FieldStruct;
  AttributesByName[Field<isExternal>::Name()]        = new FieldStructWrapper<isExternal>::FieldStruct;
  AttributesByName[Field<processColor>::Name()]        = new FieldStructWrapper<processColor>::FieldStruct;
  AttributesByName[Field<ghostRank>::Name()]        = new FieldStructWrapper<ghostRank>::FieldStruct;
  AttributesByName[Field<ownedByRank>::Name()]        = new FieldStructWrapper<ownedByRank>::FieldStruct;
  AttributesByName[Field<demIndex>::Name()]                 = new FieldStructWrapper<demIndex>::FieldStruct;



  AttributesByName[Field<referencePosition>::Name()]        = new FieldStructWrapper<referencePosition>::FieldStruct;
  AttributesByName[Field<currentPosition>::Name()]          = new FieldStructWrapper<currentPosition>::FieldStruct;
  AttributesByName[Field<displacement>::Name()]             = new FieldStructWrapper<displacement>::FieldStruct;
  AttributesByName[Field<incrementalDisplacement>::Name()]  = new FieldStructWrapper<incrementalDisplacement>::FieldStruct;
  AttributesByName[Field<relativePosition>::Name()]         = new FieldStructWrapper<relativePosition>::FieldStruct;
  AttributesByName[Field<velocity>::Name()]                 = new FieldStructWrapper<velocity>::FieldStruct;
  AttributesByName[Field<acceleration>::Name()]             = new FieldStructWrapper<acceleration>::FieldStruct;
  AttributesByName[Field<force>::Name()]                    = new FieldStructWrapper<force>::FieldStruct;
  AttributesByName[Field<hgforce>::Name()]                    = new FieldStructWrapper<hgforce>::FieldStruct;
  AttributesByName[Field<contactForce>::Name()]                    = new FieldStructWrapper<contactForce>::FieldStruct;
  AttributesByName[Field<mass>::Name()]                     = new FieldStructWrapper<mass>::FieldStruct;
  AttributesByName[Field<rotationAxis>::Name()]             = new FieldStructWrapper<rotationAxis>::FieldStruct;
  AttributesByName[Field<rotationMagnitude>::Name()]        = new FieldStructWrapper<rotationMagnitude>::FieldStruct;
  AttributesByName[Field<rotationalAxisIncrement>::Name()]             = new FieldStructWrapper<rotationalAxisIncrement>::FieldStruct;
  AttributesByName[Field<rotationalMagnitudeIncrement>::Name()]        = new FieldStructWrapper<rotationalMagnitudeIncrement>::FieldStruct;
  AttributesByName[Field<rotationalVelocity>::Name()]        = new FieldStructWrapper<rotationalVelocity>::FieldStruct;
  AttributesByName[Field<rotationalAcceleration>::Name()]        = new FieldStructWrapper<rotationalAcceleration>::FieldStruct;
  AttributesByName[Field<moment>::Name()]        = new FieldStructWrapper<moment>::FieldStruct;
  AttributesByName[Field<rotationalInertia>::Name()]        = new FieldStructWrapper<rotationalInertia>::FieldStruct;



  AttributesByName[Field<pressure>::Name()]                 = new FieldStructWrapper<pressure>::FieldStruct;
  AttributesByName[Field<fluidPressure>::Name()]            = new FieldStructWrapper<fluidPressure>::FieldStruct;
  AttributesByName[Field<deviatorStress>::Name()]           = new FieldStructWrapper<deviatorStress>::FieldStruct;
  AttributesByName[Field<volume>::Name()]                   = new FieldStructWrapper<volume>::FieldStruct;
  AttributesByName[Field<density>::Name()]                  = new FieldStructWrapper<density>::FieldStruct;
  AttributesByName[Field<fluidDensity>::Name()]             = new FieldStructWrapper<fluidDensity>::FieldStruct;
  AttributesByName[Field<massFlux>::Name()]                 = new FieldStructWrapper<massFlux>::FieldStruct;
  AttributesByName[Field<damageIndicator>::Name()]          = new FieldStructWrapper<damageIndicator>::FieldStruct;

  AttributesByName[Field<seismicBeginTime>::Name()]         = new FieldStructWrapper<seismicBeginTime>::FieldStruct;
  AttributesByName[Field<seismicRiseTime>::Name()]          = new FieldStructWrapper<seismicRiseTime>::FieldStruct;
  AttributesByName[Field<slip>::Name()]                     = new FieldStructWrapper<slip>::FieldStruct;
  AttributesByName[Field<area>::Name()]                     = new FieldStructWrapper<area>::FieldStruct;
  AttributesByName[Field<seismicMagnitude>::Name()]         = new FieldStructWrapper<seismicMagnitude>::FieldStruct;
  AttributesByName[Field<seismicMoment>::Name()]            = new FieldStructWrapper<seismicMoment>::FieldStruct;
}

void FieldInfo::DeleteAttributes()
{
  for( std::vector<FieldBase*>::iterator i=AttributesByKey.begin() ; i!=AttributesByKey.end() ; ++i )
    delete (*i);

  for( std::map<std::string, FieldBase*>::iterator i=AttributesByName.begin() ; i!=AttributesByName.end() ; ++i )
    delete (i->second);
}
