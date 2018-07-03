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
 * JointSet.h
 *
 *  Created on: Jun 7, 2012
 *      Author: johnson346
 */

#ifndef JOINTSET_H_
#define JOINTSET_H_

#include "../IO/ticpp/HierarchicalDataNode.h"
#include "Common/Common.h"
#include "FractalVolume.h"
#include "StatisticalDistributionBaseT.h"

class JointSetT
{
public:
  JointSetT();
  virtual ~JointSetT();

  virtual void
  ReadXML(TICPP::HierarchicalDataNode* hdn);

  void NextStrikeDipNormal(R1Tensor& strikeVector,
                           R1Tensor& dipVector,
                           R1Tensor& normalVector);

  realT NextStrikeLengthPowerLaw() {
    return m_strikeDimensionDistribution.PowerLawSample();
  }

  realT NextAspectRatioGaussian() { return m_faultAspectRatioDistribution.NormalSample(); }

  void SamplePositionsFractal(const array<R1Tensor>& ref,
                              const array<R1Tensor>& disp,
                              const R1Tensor& min,
                              const R1Tensor& max,
                              array<R1Tensor>& positions,
                              array<R1Tensor>& normals,
                              array<R1Tensor>& strikes,
                              array<R1Tensor>& dips,
                              const realT weight = 1.0);

  void SampleFrequenciesFractal(const array<R1Tensor>& centroids,
                                const gArray1d& localToGlobal,
                                const R1Tensor& min,
                                const R1Tensor& max,
                                array<real64>& frequencies);

private:
  StatisticalDistributionBaseT m_strikeAngleDistribution;
  StatisticalDistributionBaseT m_dipAngleDistribution;
  StatisticalDistributionBaseT m_strikeDimensionDistribution;
  StatisticalDistributionBaseT m_faultAspectRatioDistribution;

  ///Set FractalVolume
  StatisticalDistributionBaseT m_hypocenterDistribution;
public:
  R1Tensor m_up, m_north;
  unsigned m_seed;
};

#endif /* JOINTSET_H_ */
