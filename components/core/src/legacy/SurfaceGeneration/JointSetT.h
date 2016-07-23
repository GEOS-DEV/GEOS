/*
 * JointSet.h
 *
 *  Created on: Jun 7, 2012
 *      Author: johnson346
 */

#ifndef JOINTSET_H_
#define JOINTSET_H_

#include "Common/Common.h"
#include "IO/ticpp/HierarchicalDataNode.h"
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

  void SamplePositionsFractal(const Array1dT<R1Tensor>& ref,
                              const Array1dT<R1Tensor>& disp,
                              const R1Tensor& min,
                              const R1Tensor& max,
                              Array1dT<R1Tensor>& positions,
                              Array1dT<R1Tensor>& normals,
                              Array1dT<R1Tensor>& strikes,
                              Array1dT<R1Tensor>& dips,
                              const realT weight = 1.0);

  void SampleFrequenciesFractal(const Array1dT<R1Tensor>& centroids,
                                const gArray1d& localToGlobal,
                                const R1Tensor& min,
                                const R1Tensor& max,
                                rArray1d& frequencies);

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
