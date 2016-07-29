/*
 * JointPopulator.h
 *
 *  Created on: Nov 19, 2012
 *      Author: johnson346
 */

#ifndef JOINTPOPULATOR_H_
#define JOINTPOPULATOR_H_

#include "JointSetT.h"
#include "ObjectManagers/FaceManagerT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "DataStructures/VectorFields/ElementRegionT.h"

class JointPopulator
{
public:
  JointPopulator();
  virtual ~JointPopulator();

  virtual void
  ReadXML(TICPP::HierarchicalDataNode* hdn);

  void Initialize();
  void AddJointSet(const realT weight, const JointSetT& jointSet);

  void Populate(const ElementRegionT& elementRegion,
                const Array1dT<R1Tensor>& nodesRef,
                const Array1dT<R1Tensor>& nodesDisp,
                const R1Tensor& min,
                const R1Tensor& max,
                Array1dT<rArray1d>& frequencies);

  bool Populate(const std::map< std::string, ElementRegionT >& m_ElementRegions,
                const R1Tensor& min,
                const R1Tensor& max,
                FaceManagerT& faceManager,
                NodeManagerT& nodeManager);

  bool Populate(const std::map< std::string, ElementRegionT >& elementRegions,
                FaceManagerT& faceManager,
                NodeManagerT& nodeManager);

  bool Populate(Array1dT<R1Tensor>& positions,
                Array1dT<R1Tensor>& normals,
                Array1dT<R1Tensor>& strikes,
                Array1dT<R1Tensor>& dips) const;

  static bool Populate(const std::map< std::string, ElementRegionT >& elementRegions,
                       const sArray1d& elementRegionNames,
                       const std::string& nodeSetName,
                       Array1dT<R1Tensor>& positions,
                       Array1dT<R1Tensor>& normals,
                       Array1dT<R1Tensor>& strikes,
                       Array1dT<R1Tensor>& dips,
                       FaceManagerT& faceManager,
                       NodeManagerT& nodeManager);

  bool Next(R1Tensor& strikeVector,
            R1Tensor& dipVector,
            R1Tensor& normalVector,
            realT& strikeLength,
            realT& dipLength);

  inline R1Tensor North() const
  {
    if(m_jointSets.size() > 0)
    {
      return m_jointSets.back().m_north;
    }
    else
    {
      R1Tensor ret(0);
      ret(1) = 1.0;
      return ret;
    }
  }

  inline R1Tensor Up() const
  {
    if(m_jointSets.size() > 0)
    {
      return m_jointSets.back().m_up;
    }
    else
    {
      R1Tensor ret(0);
      ret(2) = 1.0;
      return ret;
    }
  }

  inline bool Active() const { return m_isActive; }

private:
  std::vector<JointSetT> m_jointSets;
  rArray1d m_jointWeights;
  bool m_isActive;
public:
  std::string m_elementRegionName, m_nodeSetName, m_fileName;
};

#endif /* JOINTPOPULATOR_H_ */
