/*
 * ConstitutivePropertiesTable.h
 *
 *  Created on: Mar 12, 2013
 *      Author: johnson346
 */

#ifndef CONSTITUTIVEPROPERTIESTABLE_H_
#define CONSTITUTIVEPROPERTIESTABLE_H_

#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "DataStructures/Tables/Table.h"
#include "ConstitutiveBase.h"

class ConstitutivePropertiesTable
{
public:
  ConstitutivePropertiesTable();
  virtual ~ConstitutivePropertiesTable();

  void Add(const std::string& fieldName,
           const std::string& tableName,
           ObjectDataStructureBaseT& ods,
           NodeManagerT* nm);

  void Add(const std::string& fieldName,
           const realT scalar,
           ObjectDataStructureBaseT& ods,
           NodeManagerT* nm);

  void Apply(const realT time = 0.0);

  sArray1d m_fieldNames;
  sArray1d m_tableNames;
  rArray1d m_scalars;
  std::map<ObjectDataStructureBaseT*, lArray1d> m_managerFields;

  //Have to do something special for the finite element manager
  NodeManagerT* m_feNodeManager;
  lSet m_feElementManagerIndices;
};

#endif /* CONSTITUTIVEPROPERTIESTABLE_H_ */
