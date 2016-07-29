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
/*
 * ObjectManager.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "ObjectManager.h"

/**
 * @author R.R. Settgast
 * @return
 */
ObjectManager::ObjectManager():
m_objects(),
m_numObjects(0),
m_objectMemory(),
m_requestedFields(),
m_fieldRegistry()
{}


/**
 * @author R.R. Settgast
 * @return
 */
ObjectManager::ObjectManager( const ObjectManager& init )
{

}


/**
 * @author R.R. Settgast
 * @return
 */
ObjectManager::~ObjectManager()
{
}


void ObjectManager::RegisterFields( std::vector< std::pair< FieldKey, FieldType > > requestedFields )
{
  for(std::vector< std::pair< FieldKey, FieldType > >::const_iterator field = requestedFields.begin() ;
      field != requestedFields.end() ; ++field )
  {
    m_requestedFields.insert( *field );
  }

}

void ObjectManager::OrganizeFields(  )
{
  m_fieldRegistry.SetRegistry(m_requestedFields);
}


void ObjectManager::AllocateObjectFields()
{
  m_objectMemory.resize( m_fieldRegistry.m_objectSize * m_numObjects );
}

void ObjectManager::resize( const size_t size )
{
  m_numObjects = size;
  m_objects.resize( size );


  this->AllocateObjectFields();

  realT* p_object = &(m_objectMemory[0]);
  for( std::vector<Object*>::iterator object=m_objects.begin() ; object!=m_objects.end() ; ++object )
  {
    (*object) = (Object*)p_object;
    (*object)->m_fieldRegistry = &(this->m_fieldRegistry);
    p_object += m_fieldRegistry.m_objectSize ;
  }

}
