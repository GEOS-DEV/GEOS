/*
 * MasteredBifroest.h
 *
 *  Created on: Aug 15, 2012
 *      Author: johnson346
 */

#ifndef MASTEREDBIFROEST_H_
#define MASTEREDBIFROEST_H_

#include "BifroestBase.h"

class MasteredBifroest: public BifroestBase
{
public:
  MasteredBifroest();
  virtual ~MasteredBifroest();
  virtual void CoordinateShares() {};
protected:
  virtual void ProcessPackage(ProcessingDelegatedToMe& current){};
  virtual void ProcessReturnedData(DelegateProcessingToAnother& current){};
};

#endif /* MASTEREDBIFROEST_H_ */
