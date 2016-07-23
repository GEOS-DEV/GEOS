/*
 * ConstitutivePropertiesTable.cpp
 *
 *  Created on: Mar 12, 2013
 *      Author: johnson346
 */

#include "ConstitutivePropertiesTable.h"
#include "ObjectManagers/TableManager.h"
#include "DataStructures/Tables/Table.h"
#include "DataStructures/VectorFields/ElementRegionT.h"

ConstitutivePropertiesTable::ConstitutivePropertiesTable() :
m_fieldNames(), m_tableNames(), m_scalars(), m_managerFields(),
m_feNodeManager(0), m_feElementManagerIndices()
{
}

ConstitutivePropertiesTable::~ConstitutivePropertiesTable()
{
}

void
ConstitutivePropertiesTable::Add(const std::string& fieldName,
                                 const std::string& tableName,
                                 ObjectDataStructureBaseT& ods,
                                 NodeManagerT* nm)
{
  const localIndex last = m_fieldNames.size();
  if(nm)
  {
    m_feNodeManager = nm;
    m_feElementManagerIndices.insert(last);
  }
  m_managerFields[&ods].push_back(last);
  m_fieldNames.push_back(fieldName);
  m_tableNames.push_back(tableName);
  m_scalars.push_back(0.0);
}

void
ConstitutivePropertiesTable::Add(const std::string& fieldName,
                                 const realT scalar,
                                 ObjectDataStructureBaseT& ods,
                                 NodeManagerT* nm)
{
  const localIndex last = m_fieldNames.size();
  if(nm)
  {
    m_feNodeManager = nm;
    m_feElementManagerIndices.insert(last);
  }
  m_managerFields[&ods].push_back(last);
  m_fieldNames.push_back(fieldName);
  m_tableNames.push_back("");
  m_scalars.push_back(scalar);
}

void ConstitutivePropertiesTable::Apply(const realT time)
{
  const TableManager& tableManager = TableManager::Instance();

  //foreach name and manager pair
  for(std::map<ObjectDataStructureBaseT*, lArray1d>::iterator im = m_managerFields.begin(); im != m_managerFields.end(); ++im)
  {
    Array1dT<rArray1d> fields;
    sArray1d names;
    ObjectDataStructureBaseT& ods = *(im->first);
    for(lArray1d::const_iterator it = im->second.begin(); it != im->second.end(); ++it)
    {
      names.push_back(m_fieldNames[*it]);
      const std::string& tableName = m_tableNames[*it];

      //fill the temporary array
      rArray1d field;
      if(!tableName.empty())
      {
        field.resize(ods.DataLengths());

        //get the table to query
        const Table3D* t3dp = stlMapLookupPointer(tableManager.Tables<3>(), tableName);
        if(!t3dp)
          throw GPException(
              "ConstitutivePropertiesTable::Apply : Cannot find requested table in the table manager: " + tableName);

        if (m_feElementManagerIndices.find(*it) != m_feElementManagerIndices.end() )
        {
          for (localIndex i = 0; i < ods.DataLengths(); ++i)
          {
            R1Tensor cpos(((ElementRegionT&) ods).GetElementCenter(i, *m_feNodeManager));
            field[i] = t3dp->Lookup(cpos); // , TableInterpolation::linear
          }
        }
        else
        {
          localIndex i = 0;
          const Array1dT<R1Tensor>& pos = ods.GetFieldData<FieldInfo::referencePosition>();
          for (Array1dT<R1Tensor>::const_iterator it4 = pos.begin(); it4 != pos.end(); ++it4, ++i)
            field[i] = t3dp->Lookup(*it4); // , TableInterpolation::linear
        }
      }
      else
      {
        //fill "field" from a scalar
        field.resize(ods.DataLengths(), m_scalars[*it]);
      }
      fields.push_back(field);
    }

    //copy the field entries into the appropriate object members
    ods.DeserializeObjectFields(names, fields);
  }
}

