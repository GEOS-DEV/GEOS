/*
 * StableTimeStep.h
 *
 *  Created on: Jan 24, 2012
 *      Author: scottjohnson
 */

#ifndef STABLETIMESTEP_H_
#define STABLETIMESTEP_H_

struct StableTimeStep
{
//  std::string m_region;
//  localIndex m_index;
  realT m_maxdt;

  inline bool SetIfSmaller(const realT dt)
  {
    if(dt < m_maxdt) {
      this->m_maxdt = dt;
      return true;
    } else {
      return false;
    }
  }
};

#endif /* STABLETIMESTEP_H_ */
