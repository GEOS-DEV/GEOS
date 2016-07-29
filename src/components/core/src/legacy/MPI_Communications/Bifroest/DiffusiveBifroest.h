/*
 * DiffusiveBifroest.h
 *
 *  Created on: Aug 15, 2012
 *      Author: johnson346
 */

#ifndef DIFFUSIVEBIFROEST_H_
#define DIFFUSIVEBIFROEST_H_

#include "BifroestBase.h"

class DiffusiveBifroest: public BifroestBase
{
public:
  DiffusiveBifroest();
  virtual ~DiffusiveBifroest();
  virtual void CoordinateShares() {};
protected:
  virtual void ProcessPackage(ProcessingDelegatedToMe& current){};
  virtual void ProcessReturnedData(DelegateProcessingToAnother& current){};
};

#endif /* DIFFUSIVEBIFROEST_H_ */
