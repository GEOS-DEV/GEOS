/*
 * RestartFile.h
 *
 *  Created on: Jul 18, 2012
 *      Author: settgast1
 */

#ifndef RESTARTFILE_H_
#define RESTARTFILE_H_

#include "PMPIOBase.h"
#include "BinStream.h"

class RestartFile: public PMPIOBase<BinStream>
{
public:
  RestartFile();
  ~RestartFile();

  void Initialize(const PMPIO_iomode_t readwrite);

};

#endif /* RESTARTFILE_H_ */
