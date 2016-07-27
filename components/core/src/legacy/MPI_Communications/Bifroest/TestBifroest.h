/*
 * TestBifroest.h
 *
 *  Created on: Aug 29, 2012
 *      Author: johnson346
 */

#ifndef TESTPROCESSSHARE_H_
#define TESTPROCESSSHARE_H_

#include "BifroestBase.h"

class TestBifroest: public BifroestBase
{
public:
  TestBifroest();
  virtual ~TestBifroest();
protected:
  virtual void ProcessPackage(ProcessingDelegatedToMe& current);
  virtual void ProcessReturnedData(DelegateProcessingToAnother& current);
};


#endif /* TESTPROCESSSHARE_H_ */
