/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * EnergyT.cpp
 *
 *  Created on: Jan 30, 2012
 *      Author: johnson346
 */

#include "EnergyT.h"
EnergyT::EnergyT():
  m_strainEnergy(0.0),
  m_strainEnergyMax(0.0),
  m_stressPower(0.0),
  m_dissipatedEnergy(0.0),
  m_kineticEnergy(0.0),
  m_workDoneOnNodes(0.0)
{
  // TODO Auto-generated constructor stub
}

EnergyT::~EnergyT()
{
  // TODO Auto-generated destructor stub
}

void EnergyT::Increment(const EnergyT& energy)
{
  this->IncrementDissipatedEnergy(energy.DissipatedEnergy());
  this->IncrementStrainEnergy(energy.StrainEnergy());
  this->IncrementStressPower(energy.StressPower());
}

void EnergyT::Serialize(double* const buffer) const
{
  buffer[0] = this->m_dissipatedEnergy;
  buffer[1] = this->m_strainEnergy;
  buffer[2] = this->m_strainEnergyMax;
  buffer[3] = this->m_stressPower;
  buffer[4] = this->m_kineticEnergy;
  buffer[5] = this->m_workDoneOnNodes;
}

void EnergyT::Deserialize(const double* const buffer)
{
  this->m_dissipatedEnergy = buffer[0];
  this->m_strainEnergy     = buffer[1];
  this->m_strainEnergyMax  = buffer[2];
  this->m_stressPower      = buffer[3];
  this->m_kineticEnergy    = buffer[4];
  this->m_workDoneOnNodes  = buffer[5];
}

void EnergyT::Increment(const double* const buffer)
{
  this->IncrementDissipatedEnergy(buffer[0]);
  this->IncrementStrainEnergy(buffer[1]);
  this->IncrementStressPower(buffer[2]);
}
