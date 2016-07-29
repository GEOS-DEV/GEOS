/*
 * DistributedBifroest.h
 *
 *  Created on: Aug 15, 2012
 *      Author: johnson346
 */

#ifndef DISTRIBUTEDBIFROEST_H_
#define DISTRIBUTEDBIFROEST_H_

#include "BifroestBase.h"

class RankCurrentContactsTuple
{
public:
  int rank;
  int status;
  RankCurrentContactsTuple() : rank(0), status(0) {};
  ~RankCurrentContactsTuple() {};
};

class DistributedBifroest: public BifroestBase
{
public:
  DistributedBifroest();
  virtual ~DistributedBifroest();
  virtual void CoordinateShares();
protected:
  virtual void ProcessPackage(ProcessingDelegatedToMe& current){};
  virtual void ProcessReturnedData(DelegateProcessingToAnother& current){};

  static bool CompareDecreasing(const RankCurrentContactsTuple* const ta, const RankCurrentContactsTuple* const tb)
  {
    return ta->status > tb->status;
  };

  void DetermineInteractionPartitionMapMarchingMatch(localIndex& maxPerNode,
                                                     localIndex& totalCommunications,
                                                     std::map<int, iArray1d>& demandsToSupplies);
  void DetermineInteractionPartitionMapCluster(localIndex& maxPerNode,
                                               localIndex& totalCommunications,
                                               std::map<int, iArray1d>& demandsToSupplies);

  int Status(const int numberOfContacts) const;
  int floorMean, ceilingMax, ceilingMean, status;
  std::vector<RankCurrentContactsTuple> demands, supplies;
  iArray1d allCurrentContacts;
};

#endif /* DISTRIBUTEDBIFROEST_H_ */
