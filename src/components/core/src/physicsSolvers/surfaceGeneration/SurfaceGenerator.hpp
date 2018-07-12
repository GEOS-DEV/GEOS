/*
 * SurfaceGenerator.hpp
 *
 *  Created on: Jul 3, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_
#define SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

class SpatialPartition;
class ModifiedObjectLists;

class NodeManager;
class EdgeManager;
class FaceManager;
class ExternalFaceManager;
class ElementRegionManager;
class ElementRegion;


class SurfaceGenerator : public SolverBase
{
public:
  SurfaceGenerator( const std::string& name,
                    ManagedGroup * const parent );
  ~SurfaceGenerator();


  static string CatalogName() { return "SurfaceGenerator"; }


  virtual void TimeStep( real64 const & time_n,
                         real64 const & dt,
                         int const cycleNumber,
                         dataRepository::ManagedGroup * domain ) override;

  int SeparationDriver( NodeManager& nodeManager,
                        EdgeManager& edgeManager,
                        FaceManager& faceManager,
                        ExternalFaceManager& externalFaceManager,
                        ElementRegionManager& elementManager,
                        SpatialPartition& partition,
                        const bool prefrac,
                        const realT time );
private:
  integer m_failCriterion;
};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_ */
