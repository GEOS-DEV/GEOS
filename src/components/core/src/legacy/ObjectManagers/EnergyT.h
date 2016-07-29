/*
 * EnergyT.h
 *
 *  Created on: Jan 30, 2012
 *      Author: johnson346
 */

#ifndef ENERGYT_H_
#define ENERGYT_H_

#include "Common/intrinsic_typedefs.h"

class EnergyT
{
public:
  EnergyT();
  virtual ~EnergyT();

  /// assignment
  EnergyT& operator=( const EnergyT& rhs )
  {
    m_strainEnergy     = rhs.m_strainEnergy;
    m_strainEnergyMax  = rhs.m_strainEnergyMax;
    m_stressPower      = rhs.m_stressPower;
    m_dissipatedEnergy = rhs.m_dissipatedEnergy;
    m_kineticEnergy    = rhs.m_kineticEnergy;
    m_workDoneOnNodes  = rhs.m_workDoneOnNodes;
    return *this;
  }

  EnergyT& operator=( const realT rhs )
  {
    m_strainEnergy     = rhs;
    m_strainEnergyMax  = rhs;
    m_stressPower      = rhs;
    m_dissipatedEnergy = rhs;
    m_kineticEnergy    = rhs;
    m_workDoneOnNodes  = rhs;
    return *this;
  }


  EnergyT& operator+=( const EnergyT& rhs )
  {
    m_strainEnergy     += rhs.m_strainEnergy;
    m_stressPower      += rhs.m_stressPower;
    m_dissipatedEnergy += rhs.m_dissipatedEnergy;
    m_kineticEnergy    += rhs.m_kineticEnergy;
    if(this->m_strainEnergy > this->m_strainEnergyMax)
      this->m_strainEnergyMax = this->m_strainEnergy;
    m_workDoneOnNodes  += rhs.m_workDoneOnNodes;
    return *this;
  }

  void Zero()
  {
    m_strainEnergy     = 0.0;
    m_strainEnergyMax  = 0.0;
    m_stressPower      = 0.0;
    m_dissipatedEnergy = 0.0;
    m_kineticEnergy    = 0.0;
    m_workDoneOnNodes  = 0.0;
  }

  inline void IncrementStrainEnergy(const realT dE)
  {
    this->m_strainEnergy += dE;
    if(this->m_strainEnergy > this->m_strainEnergyMax)
      this->m_strainEnergyMax = this->m_strainEnergy;
  }

  inline void IncrementDissipatedEnergy(const realT dE)
  {
    this->m_dissipatedEnergy += dE;
  }

  inline void IncrementStressPower(const realT dE)
  {
    this->m_stressPower += dE;
  }

  inline void SetWorkDoneOnNodes(const realT dE)
  {
    this->m_workDoneOnNodes = dE;
  }

  inline void SetKineticEnergy( const realT KE )
  {
    m_kineticEnergy = KE;
  }

  enum
  {
   numVars = 6
  };


  void Increment(const EnergyT& energy);
  void Serialize(double* const buffer) const;
  void Deserialize(const double* const buffer);
  void Increment(const double* const buffer);

  inline realT StrainEnergy() const { return m_strainEnergy;}
  inline realT StressPower() const { return m_stressPower;}
  inline realT DissipatedEnergy() const { return m_dissipatedEnergy;}
  inline realT KineticEnergy() const { return m_kineticEnergy; }
  inline realT WorkDoneOnNodes() const { return m_workDoneOnNodes; }



private:
  realT m_strainEnergy;
  realT m_strainEnergyMax;
  realT m_stressPower;
  realT m_dissipatedEnergy;
  realT m_kineticEnergy;
  realT m_workDoneOnNodes;
};

#endif /* ENERGYT_H_ */
