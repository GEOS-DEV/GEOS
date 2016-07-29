/*
 * JointPopulator.cpp
 *
 *  Created on: Nov 19, 2012
 *      Author: johnson346
 */

#include "JointPopulator.h"
#include "ObjectManagers/PhysicalDomainT.h"

#if GPAC_MPI
#include <mpi.h>
#endif

JointPopulator::JointPopulator() : m_jointSets(), m_jointWeights(), m_isActive(false), m_elementRegionName(""), m_nodeSetName(""), m_fileName("")
{
}

JointPopulator::~JointPopulator()
{
}

void
JointPopulator::Initialize()
{
  realT weightTotal = 0.0;
  for(rArray1d::const_iterator it = m_jointWeights.begin(); it != m_jointWeights.end(); ++it)
    weightTotal += *it;
  if(!isEqual(weightTotal, 1.0) && !isZero(weightTotal))
  {
    for(rArray1d::iterator it = m_jointWeights.begin(); it != m_jointWeights.end(); ++it)
      *it /= weightTotal;
  }
}

void
JointPopulator::ReadXML(TICPP::HierarchicalDataNode* hdn)
{
  R1Tensor north(0.0), up(0.0);

  ///element region name
  m_elementRegionName = hdn->GetAttributeStringOrDefault("elementRegionName", "");
  m_nodeSetName = hdn->GetAttributeStringOrDefault("nodeSetName", "default_joint_nodeset");
  m_fileName = hdn->GetAttributeString("file");
  m_isActive = true;

  //Orientation
  {
    ///Up vector
    R1Tensor upDef(0.0);
    upDef(2) = 1.0;
    up = hdn->GetAttributeTensorOrDefault("upVector", upDef);
    if (isEqual(Dot(up, up), 0))
      up = upDef;
    else
      up.Normalize();

    ///North vector
    R1Tensor northDef(0.0);
    northDef(1) = 1.0;
    north = hdn->GetAttributeTensorOrDefault("northVector", northDef);
    if (isEqual(Dot(north, north), 0))
      north = northDef;
    else
      north.Normalize();
  }

  for(TICPP::HierarchicalDataNode* jointNode = hdn->Next(true); jointNode; jointNode = hdn->Next())
  {
    m_jointSets.resize(m_jointSets.size()+1);
    m_jointSets.back().m_north = north;
    m_jointSets.back().m_up = up;
    m_jointSets.back().ReadXML(jointNode);
    m_jointWeights.push_back(jointNode->GetAttributeOrDefault("weight", 1.0));
  }
  Initialize();
}

void
JointPopulator::AddJointSet(const realT weight, const JointSetT& jointSet)
{
  m_jointSets.push_back(jointSet);
  m_jointWeights.push_back(weight);
}

bool JointPopulator::Next(R1Tensor& strikeVector,
                          R1Tensor& dipVector,
                          R1Tensor& normalVector,
                          realT& strikeLength,
                          realT& dipLength)
{
  if(m_jointWeights.size()==0)
    return false;

  //find the joint set to sample
  localIndex i = -1;
  {
    const realT next = StatisticalDistributionBaseT::UniformSample(0.0, 1.0);
    {
      realT cumWgt = 0.0;
      for(i = 0; i < m_jointWeights.size(); ++i)
      {
        cumWgt += m_jointWeights[i];
        if(cumWgt >= next)
          break;
      }
    }
  }

  //sample the joint set
  JointSetT& js = m_jointSets[i];
  js.NextStrikeDipNormal(strikeVector,
                         dipVector,
                         normalVector);
  strikeLength = js.NextStrikeLengthPowerLaw(); //FIXME!!

  dipLength = js.NextAspectRatioGaussian() * strikeLength;

  return true;
}

///Creates a joint distribution
/**
 * This function creates a realization of joint hypocenters and bins them to
 * elements within the specified element region; the result is a list of frequencies
 * of joint hypocenters for each element
 */
void
JointPopulator::Populate(const ElementRegionT& elementRegion,
                         const Array1dT<R1Tensor>& nodesRef,
                         const Array1dT<R1Tensor>& nodesDisp,
                         const R1Tensor& min,
                         const R1Tensor& max,
                         Array1dT<rArray1d>& frequencies)
{
  const gArray1d& localToGlobal = elementRegion.m_localToGlobalMap;

  //GET EVERY ELEMENT'S CENTROID
  Array1dT<R1Tensor> centroids(localToGlobal.size(), static_cast<R1Tensor>(0.0) );
  {
    const FixedOneToManyRelation& elementToNodes = elementRegion.m_toNodesRelation;
    for (localIndex i = 0; i < elementToNodes.Dimension(0); ++i)
    {
      R1Tensor& center = centroids[i];
      for (localIndex j = 0; j < elementToNodes.Dimension(1); ++j)
      {
        center += nodesRef[elementToNodes(i, j)];
        center += nodesDisp[elementToNodes(i, j)];
      }
      center *= 1.0 / elementToNodes.Dimension(1);
    }
  }

  //DO THE SPATIAL DISTRIBUTION OF PROPERTIES
  frequencies.resize(m_jointSets.size());
  for(localIndex i = 0; i < m_jointSets.size(); i++)
  {
    frequencies[i].resize(centroids.size(), 0.0);
    m_jointSets[i].SampleFrequenciesFractal(centroids, localToGlobal, min, max, frequencies[i]);
  }
  //frequencies is now filled with the joint counts for all elements in the element region on the local process
}

///Creates a joint distribution
/**
 * This function snaps a statistical realization of joints to the closest faces and then defines
 * the corresponding node set
 */
bool
JointPopulator::Populate(const std::map< std::string, ElementRegionT >& elementRegions,
                         const R1Tensor& min,
                         const R1Tensor& max,
                         FaceManagerT& faceManager,
                         NodeManagerT& nodeManager)
{
  if(!m_isActive || m_nodeSetName.length() == 0)
    return false;

  sArray1d elementRegionNames;
  if(m_elementRegionName.length() == 0)
  {
    elementRegionNames.reserve(elementRegions.size());
    for(std::map< std::string, ElementRegionT >::const_iterator it = elementRegions.begin(); it != elementRegions.end(); ++it)
      elementRegionNames.push_back(it->first);
  }
  else
  {
    elementRegionNames.reserve(1);
    elementRegionNames.push_back(m_elementRegionName);
  }

  //get position references
  const Array1dT<R1Tensor>& ref = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& disp = nodeManager.GetFieldData<FieldInfo::displacement>();

  //first, generate the spatial distribution of joints and get normals/dimensions
  Array1dT<R1Tensor> positions, normals, strikes, dips;
  for(localIndex ijs = 0; ijs < m_jointSets.size(); ++ijs)
  {
    JointSetT& js = m_jointSets[ijs];
    js.SamplePositionsFractal(ref, disp, min, max,
                              positions, normals, strikes,
                              dips, m_jointWeights[ijs]);
  }

  std::cout<<"--->I made "<<positions.size()<<" joints"<<std::endl;
  return Populate(elementRegions, elementRegionNames, m_nodeSetName, positions, normals, strikes, dips, faceManager, nodeManager);
}

bool JointPopulator::Populate(Array1dT<R1Tensor>& positions,
                              Array1dT<R1Tensor>& normals,
                              Array1dT<R1Tensor>& strikes,
                              Array1dT<R1Tensor>& dips) const
{
  std::ifstream geometry;
  geometry.open(m_fileName.c_str());
  if(!geometry.is_open())
    return false;
  while(!geometry.eof())
  {
    //realT xx;
    R1Tensor x, n;
    geometry >> x(0) >> x(1) >> x(2) >> n(0) >> n(1) >> n(2);
    n.Normalize();
    if(isZero(Dot(n,n)))
      break;

    R1Tensor dip, strike;
    {
      realT mstrike, mdip;
      geometry >> mstrike >> mdip;

      GeometryUtilities::Strike(n, Up(), strike);
      dip.Cross(strike, n);
      dip.Normalize();
      dip *= mdip;
      strike *= mstrike;
    }


    positions.push_back(x);
    normals.push_back(n);
    strikes.push_back(strike);
    dips.push_back(dip);
  }
  geometry.close();
  return true;
}

bool JointPopulator::Populate(const std::map<std::string, ElementRegionT>& elementRegions,
                              FaceManagerT& faceManager,
                              NodeManagerT& nodeManager)
{
  if(!m_isActive || m_nodeSetName.length() == 0 || m_fileName.empty())
    return false;

  sArray1d elementRegionNames;
  if(m_elementRegionName.length() == 0)
  {
    elementRegionNames.reserve(elementRegions.size());
    for(std::map< std::string, ElementRegionT >::const_iterator it = elementRegions.begin(); it != elementRegions.end(); ++it)
      elementRegionNames.push_back(it->first);
  }
  else
  {
    elementRegionNames.reserve(1);
    elementRegionNames.push_back(m_elementRegionName);
  }

  //get position references
  //const Array1dT<R1Tensor>& ref = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  //const Array1dT<R1Tensor>& disp = nodeManager.GetFieldData<FieldInfo::displacement>();

  //first, generate the spatial distribution of joints and get normals/dimensions
  Array1dT<R1Tensor> positions, normals, strikes, dips;
  Populate(positions, normals, strikes, dips);
  return Populate(elementRegions, elementRegionNames, m_nodeSetName,
                  positions, normals, strikes, dips,
                  faceManager, nodeManager);
}


///Creates a joint distribution
/**
 * This function snaps a set of joint patches to the closest faces and then defines
 * the corresponding node set
 */
bool
JointPopulator::Populate(const std::map< std::string, ElementRegionT >& elementRegions,
                         const sArray1d& elementRegionNames,
                         const std::string& nodeSetName,
                         Array1dT<R1Tensor>& positions,
                         Array1dT<R1Tensor>& normals,
                         Array1dT<R1Tensor>& strikes,
                         Array1dT<R1Tensor>& dips,
                         FaceManagerT& faceManager,
                         NodeManagerT& nodeManager)
{
  //get position references
  const Array1dT<R1Tensor>& ref = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& disp = nodeManager.GetFieldData<FieldInfo::displacement>();

  lSet nodes;
  for(sArray1d::const_iterator it = elementRegionNames.begin(); it != elementRegionNames.end(); ++it)
  {
    std::map< std::string, ElementRegionT >::const_iterator iter = elementRegions.find(*it);
    if(iter == elementRegions.end())
      return false;

    const ElementRegionT& elemRegion = iter->second;
    const localIndex nnodes = elemRegion.m_toNodesRelation.Dimension(1);
    if(nnodes == 0)
      return false;
    const localIndex nfaces = elemRegion.m_toFacesRelation.Dimension(1);
    if(nfaces == 0)
      return false;

    //second, cut through elements and find possible nodes to sever
    lSet nodeCandidates;
    for (localIndex iel = 0; iel < elemRegion.m_numElems; ++iel)
    {
      //get the element center
      R1Tensor xc = elemRegion.GetElementCenter( iel, nodeManager );

      for(localIndex ijoint = 0; ijoint < positions.size(); ++ijoint)
      {
        //skip the element if the projection of the center does not lie on the joint patch;
        //this (kind of) assumes that if a joint cuts more than halfway through any element
        //it cuts all the way through; otherwise, it does not cut at all
        R1Tensor dx(positions[ijoint]);
        dx -= xc;
        if(fabs(Dot(dx,strikes[ijoint]) / Dot(strikes[ijoint],strikes[ijoint])) > 1)
          continue;
        else if(fabs(Dot(dx,dips[ijoint]) / Dot(dips[ijoint],dips[ijoint])) > 1)
          continue;

        lSet candidates;
        for(localIndex i = 0; i < nnodes; i++)
        {
          dx = ref[elemRegion.m_toNodesRelation(iel, i)];
          dx += disp[elemRegion.m_toNodesRelation(iel, i)];
          dx -= positions[ijoint];
          if(Dot(dx, normals[ijoint]) < 0)
            candidates.insert(elemRegion.m_toNodesRelation(iel, i));
        }
        if(candidates.size() != nnodes)
          nodeCandidates.insert(candidates.begin(), candidates.end());
      }
    }

    //make sure all nodes on the face are ready to split before a node is added
    for (localIndex iel = 0; iel < elemRegion.m_numElems; ++iel)
    {
      for(localIndex i = 0; i < nfaces; i++)
      {
        const localIndex iface = elemRegion.m_toFacesRelation(iel, i);
        const lArray1d& inodes = faceManager.m_toNodesRelation[iface];
        bool ok = true;
        for(lArray1d::const_iterator itn = inodes.begin(); itn != inodes.end(); ++itn)
        {
          if(nodeCandidates.find(*itn) == nodeCandidates.end())
          {
            ok = false;
            break;
          }
        }
        if(ok)
          nodes.insert(inodes.begin(), inodes.end());
      }
    }

    //add nodes into set
    nodeManager.m_Sets[nodeSetName].insert(nodes.begin(), nodes.end());
  }

  faceManager.ConstructSetFromSetAndMap( nodeManager.m_Sets[nodeSetName],
                                         faceManager.m_toNodesRelation,
                                         nodeSetName );



  std::cout<<"--->I made "<<nodes.size()<<" nodes altogether"<<std::endl;
  return nodes.size() > 1;
}
