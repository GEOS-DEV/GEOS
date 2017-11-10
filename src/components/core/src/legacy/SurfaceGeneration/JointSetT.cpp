/*
 * JointSet.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: johnson346
 */

#include "JointSetT.h"
#include "Utilities/GeometryUtilities.h"
#if GPAC_MPI
#include <mpi.h>
#endif

JointSetT::JointSetT() : m_strikeAngleDistribution(), m_dipAngleDistribution(), m_strikeDimensionDistribution(),
                         m_faultAspectRatioDistribution(), m_hypocenterDistribution(), m_up(0.0), m_north(0.0), m_seed(0)

{
  // TODO Auto-generated constructor stub
}

JointSetT::~JointSetT()
{
  // TODO Auto-generated destructor stub
}

void JointSetT::ReadXML( TICPP::HierarchicalDataNode* hdn )
{
  const realT toRad = atan(1) / 45.0;

  //Random seed
  m_seed = hdn->GetAttributeOrDefault<unsigned>("seed", time(NULL));

  //Orientation
  {
	  /**
    ///Up vector
    R1Tensor upDef(0.0);
    upDef(2) = 1.0;
    m_up = hdn->GetAttributeTensorOrDefault("upVector", upDef);
    if (isEqual(Dot(m_up, m_up), 0))
      m_up = upDef;
    else
      m_up.Normalize();

    ///North vector
    R1Tensor northDef(0.0);
    northDef(1) = 1.0;
    m_north = hdn->GetAttributeTensorOrDefault("northVector", northDef);
    if (isEqual(Dot(m_north, m_north), 0))
      m_north = northDef;
    else
      m_north.Normalize();
      **/

    ///Up vector
    if (isEqual(Dot(m_up, m_up), 0)){ // i.e. m_up not defined
      m_up  = R1Tensor(0,0,1); // set default
    }
    m_up = hdn->GetAttributeTensorOrDefault("upVector", m_up); // local override of m_up
    m_up.Normalize();

    ///North vector
    if (isEqual(Dot(m_north, m_north), 0)){ // i.e. m_up not defined
    	m_north  = R1Tensor(0,1,0); // set default
    }
    m_north = hdn->GetAttributeTensorOrDefault("northVector", m_north); // local override of m_up
    m_north.Normalize();

    if( fabs( Dot(m_up,m_north) ) > 1e-8 ) {
        throw GPException("Error JointSetT: Up vector and North vector are not orthogonal.");
    }
  }

  //Strike length distribution
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("exponentStrikeDimension", -3.1);
    this->m_strikeDimensionDistribution.AddParameter(StatisticalDistributionBaseT::POWERLAW_EXPONENT, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minStrikeDimension", 0.0001);
    this->m_strikeDimensionDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxStrikeDimension", 1.0);
    this->m_strikeDimensionDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }

  //Dip length distribution
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("meanAspectRatio", 0.01);
    this->m_faultAspectRatioDistribution.AddParameter(StatisticalDistributionBaseT::MEAN, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("stdevAspectRatio", 0.1);
    this->m_faultAspectRatioDistribution.AddParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minAspectRatio", 0.0001);
    this->m_faultAspectRatioDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxAspectRatio", 1.0);
    this->m_faultAspectRatioDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }

  //Strike angle distribution
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("meanStrikeAngle", 90); tmp *= toRad;
    this->m_strikeAngleDistribution.AddParameter(StatisticalDistributionBaseT::MEAN, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("stdevStrikeAngle", 8); tmp *= toRad;
    this->m_strikeAngleDistribution.AddParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minStrikeAngle", 70); tmp *= toRad;
    this->m_strikeAngleDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxStrikeAngle", 130); tmp *= toRad;
    this->m_strikeAngleDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }

  //Dip angle distribution
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("meanDipAngle", 77); tmp *= toRad;
    this->m_dipAngleDistribution.AddParameter(StatisticalDistributionBaseT::MEAN, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("stdevDipAngle", 7); tmp *= toRad;
    this->m_dipAngleDistribution.AddParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minDipAngle", 45); tmp *= toRad;
    this->m_dipAngleDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxDipAngle", 90); tmp *= toRad;
    this->m_dipAngleDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }

  //Hypocenter distribution
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("hurstExponent", 1.3);
    this->m_hypocenterDistribution.AddParameter(StatisticalDistributionBaseT::HURST_EXPONENT, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("meanFrequency", 5);
    this->m_hypocenterDistribution.AddParameter(StatisticalDistributionBaseT::MEAN, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("stdevFrequency", 4);
    this->m_hypocenterDistribution.AddParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION, tmp);

    this->m_hypocenterDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, 0.0);

    tmp = hdn->GetAttributeOrDefault<realT>("maxFrequency", 20);
    this->m_hypocenterDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }
}

void JointSetT::SamplePositionsFractal(const array<R1Tensor>& ref,
                                       const array<R1Tensor>& disp,
                                       const R1Tensor& min,
                                       const R1Tensor& max,
                                       array<R1Tensor>& positions,
                                       array<R1Tensor>& normals,
                                       array<R1Tensor>& strikes,
                                       array<R1Tensor>& dips,
                                       const realT weight)
{
  FractalVolume spatial;

  //----------GET EXTREMA
  R1Tensor dxv;
  {
    //GET PARAMETERS
    const realT* hurst = this->m_hypocenterDistribution.GetParameter(
        StatisticalDistributionBaseT::HURST_EXPONENT);
    if (!hurst)
      throw GPException("Cannot initialize fractal distribution without Hurst exponent");

    const realT* mean = this->m_hypocenterDistribution.GetParameter(StatisticalDistributionBaseT::MEAN);
    if (!mean)
      throw GPException("Cannot initialize fractal distribution without mean");

    const realT* stdev = this->m_hypocenterDistribution.GetParameter(
        StatisticalDistributionBaseT::STANDARD_DEVIATION);
    if (!stdev)
      throw GPException("Cannot initialize fractal distribution without standard deviation");

  unsigned seed;
  MPI_Allreduce(&m_seed, &seed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2, seed);
  dxv = max;
  dxv -= min;

// #if GPAC_MPI
//     int mpirank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
//     const bool root = mpirank == 0;
//     int mpisize;
//     MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
//     const bool serial = mpisize < 2;
//     if(serial)
//     {
//       spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2, this->m_seed);
//       dxv = max;
//       dxv -= min;
//     }
//     else
//     {
//       if(m_seed == 0)
//       {
//         if(root)
//         {
//           m_seed = spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2);
//           MPI_Bcast((char*)&m_seed, sizeof(unsigned), MPI_CHAR, 0, MPI_COMM_WORLD);
//         }
//         else
//         {
//           MPI_Bcast((char*)&m_seed, sizeof(unsigned), MPI_CHAR, 0, MPI_COMM_WORLD);
//           spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2, m_seed);
//         }
//       }
//       else
//       {
//         spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2, m_seed);
//       }

//       dxv = max;
//       dxv -= min;
//     }
// #else
//     spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2);
//     dxv = max;
//     dxv -= min;
// #endif
  }

  //----------GET DISTRIBUTION
  // note: hard-coding discretization to 20 per smallest dimension
  {
    array<R1Tensor> tpos;
    spatial.Positions(0.05 * dxv.MinVal(), tpos, weight);
    for(localIndex a = 0; a < tpos.size(); ++a)
    {
      R1Tensor nrm, dstrike, ddip;
      NextStrikeDipNormal(dstrike, ddip, nrm);
      ddip *= this->m_strikeDimensionDistribution.NormalSample();
      dstrike *= this->m_faultAspectRatioDistribution.NormalSample();
      const realT rs = dstrike.L2_Norm();
      const realT rd = ddip.L2_Norm();
      const realT radius = rs > rd ? 0.5 * sqrt(rs) : 0.5 * sqrt(rd);
      if( (tpos[a](0)+radius) >= min(0) && (tpos[a](0)-radius) < max(0) &&
          (tpos[a](1)+radius) >= min(1) && (tpos[a](1)-radius) < max(1) &&
          (tpos[a](2)+radius) >= min(2) && (tpos[a](2)-radius) < max(2))
      {
        positions.push_back(tpos[a]);
        normals.push_back(nrm);
        strikes.push_back(dstrike);
        dips.push_back(ddip);
      }
    }
  }
}

void JointSetT::SampleFrequenciesFractal(const array<R1Tensor>& centroids,
                                         const gArray1d& localToGlobal,
                                         const R1Tensor& min,
                                         const R1Tensor& max,
                                         array<real64>& frequencies)
{
  //GET PARAMETERS
  const realT* hurst = this->m_hypocenterDistribution.GetParameter(
      StatisticalDistributionBaseT::HURST_EXPONENT);
  if (!hurst)
    throw GPException("Cannot initialize fractal distribution without Hurst exponent");
  const realT* mean = this->m_hypocenterDistribution.GetParameter(StatisticalDistributionBaseT::MEAN);
  if (!mean)
    throw GPException("Cannot initialize fractal distribution without mean");
  const realT* stdev = this->m_hypocenterDistribution.GetParameter(
      StatisticalDistributionBaseT::STANDARD_DEVIATION);
  if (!stdev)
    throw GPException("Cannot initialize fractal distribution without standard deviation");

  //INITIALIZE FRACTAL
#if 0
  FractalVolume spatial;
  spatial.InitializeFractal(*mean, *stdev, min, max, *hurst, 3, 5, 5, 5, 1.2, this->m_seed);
#endif

  //SAMPLE
  {
    array<real64>::iterator itf = frequencies.begin();
    for (array<R1Tensor>::const_iterator itc = centroids.begin(); itc != centroids.end(); ++itc, ++itf)
    {
#if 0
      *itf = spatial.Value(*itc);
      if(*itf < 0)
      {
        throw GPException("JointSetT::SampleFrequenciesFractal : algorithm generated a negative fracture intensity!");
      }
#else
      *itf = *mean;
#endif
    }
  }
}


//void JointSetT::SampleFrequenciesGaussian(const array<real64>& centroids,
//                                          const gArray1d& localToGlobal,
//                                          const bool serial,
//                                          const globalIndex minGlobalNumber,
//                                          const globalIndex maxGlobalNumber,
//                                          array<real64>& frequencies)
//{
//  //I AM ROOT ... I HAVE TO DO ALL OF THE WORK
//
//  if (serial)
//  {
//    localIndex i = 0;
//    for (gArray1d::const_iterator iter = localToGlobal.begin(); iter != localToGlobal.end();
//        ++iter, ++i)
//    {
//      R1Tensor xneighbor;
//      for (localIndex a = 0; a < nsdof; a++)
//        xneighbor(a) = centroids[nsdof * (*iter) + a];
//      frequencies[*iter] = this->m_hypocenterDistribution.NormalSample();
//    }
//  }
//#if GPAC_MPI
//  else
//  {
//    for (globalIndex global1 = minGlobalNumber; global1 <= maxGlobalNumber; global1++) /* loop over receiver patches */
//    {
//      R1Tensor xneighbor;
//      for (localIndex a = 0; a < nsdof; a++)
//        xneighbor(a) = centroids[nsdof * global1 + a];
//      frequencies[global1] = this->m_hypocenterDistribution.NormalSample();
//    }
//
//    MPI_Bcast(frequencies.data(), frequencies.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//  }
//#else
//  else
//    throw GPException("Cannot be parallel if not compiled for MPI!!");
//#endif
//}

void JointSetT::NextStrikeDipNormal(R1Tensor& strikeVector,
                                    R1Tensor& dipVector,
                                    R1Tensor& normalVector)
{
  const realT dip = this->m_dipAngleDistribution.NormalSample();
  const realT strike = this->m_strikeAngleDistribution.NormalSample();
  R1Tensor east;
  east.Cross(m_north, m_up);
  GeometryUtilities::Strike(east, m_north, strike, strikeVector);
  GeometryUtilities::Dip(m_up, strikeVector, dip, dipVector);
  normalVector.Cross(strikeVector, dipVector);
  normalVector.Normalize();
}
