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

/**
 * @file SpatialPartition.h
 * @author settgast1
 * @date Mar 16, 2011
 */

#ifndef SPATIALPARTITION_H_
#define SPATIALPARTITION_H_

#include <map>
#include "PartitionBase.hpp"


constexpr int nsdof = 3;
namespace geosx
{
//// Planar Sorter
//// Sorts pairs of local and global indexes by the positions of their
// corresponding node points in a plane.
//class PlanarSorter {
//
//public:
//	PlanarSorter(const array<R1Tensor>& refPos, int dim) :
//			dimension(dim), refPositions(refPos) {
//	}
//	;
//
//	// sort operator for pairs containing local indexes (sort based on 1st
// element in pair)
//	bool operator()(const std::pair<localIndex, localIndex>& lhs,
//	const std::pair<localIndex, localIndex>& rhs) {
//		bool rv = false;
//		int a = 0;
//		int b = 2;
//		if (dimension == 0)
//			a = 1;
//		if (dimension == 2)
//			b = 1;
//
//		const R1Tensor& lhsVect = refPositions[lhs.first];
//		const R1Tensor& rhsVect = refPositions[rhs.first];
//
//		if (lhsVect[a] < rhsVect[a]) {
//			rv = true;
//		} else if (isEqual(lhsVect[a], rhsVect[a])
//				&& (lhsVect[b] < rhsVect[b])) {
//			rv = true;
//		};
//
//		return rv;
//	}
//	;
//
//	// sort operator for local indexes
//	bool operator()(const localIndex& lhs, const localIndex& rhs) {
//		bool rv = false;
//		int a = 0;
//		int b = 2;
//		if (dimension == 0)
//			a = 1;
//		if (dimension == 2)
//			b = 1;
//
//		const R1Tensor& lhsVect = refPositions[lhs];
//		const R1Tensor& rhsVect = refPositions[rhs];
//
//		if (lhsVect[a] < rhsVect[a]) {
//			rv = true;
//		} else if (isEqual(lhsVect[a], rhsVect[a])
//				&& (lhsVect[b] < rhsVect[b])) {
//			rv = true;
//		};
//
//		return rv;
//	}
//	;
//
//private:
//	int dimension;
//	const array<R1Tensor>& refPositions;
//};

class SpatialPartition : public PartitionBase
{
public:
  SpatialPartition();
  virtual ~SpatialPartition();

//  void ReadXML( xmlWrapper::xmlNode const & targetNode );

  void InitializePostSubGroups( dataRepository::ManagedGroup * const );
  virtual void InitializeMetis();
  void AddNeighborsMetis(gSet& neighborList);
  virtual bool IsCoordInPartition(const realT& coord, const int dir);
  virtual bool IsCoordInPartition(const R1Tensor& elemCenter);
  virtual bool IsCoordInPartition(const R1Tensor& elemCenter,
                                  const int numDistPartition);
  virtual bool IsCoordInPartitionClosed(const R1Tensor& elemCenter);
  virtual bool IsCoordInPartitionBoundingBox(const R1Tensor& elemCenter);

  virtual bool IsCoordInContactGhostRange(const R1Tensor& elemCenter);

  void setSizes(const R1Tensor& min, const R1Tensor& max);
  void setGlobalDomainSizes(const R1Tensor& min, const R1Tensor& max);
  void getSizes(R1Tensor& min, R1Tensor& max) const;
  void getPartitionSizes(R1Tensor& min, R1Tensor& max) const;
  void getPartitionGeometricalBoundary(R1Tensor& min, R1Tensor& max) const;
//  void UpdatePartitionBoundingBox(NodeManager& nodeManager);
  void GetPartitionBoundingBox(R1Tensor& xmin, R1Tensor& xmax);
  void SetPartitionGeometricalBoundary(R1Tensor& min, R1Tensor& max);

  void setPartitions(unsigned int xPartitions, unsigned int yPartitions,
                     unsigned int zPartitions) {
    m_Partitions.resize(3);
    m_Partitions(0) = xPartitions;
    m_Partitions(1) = yPartitions;
    m_Partitions(2) = zPartitions;
    m_size = 1;
    for (unsigned int i = 0 ; i < nsdof ; i++)
      m_size *= m_Partitions(i);
    SetContactGhostRange(0.0);
  }

  void setPeriodic(unsigned int xPeriodic, unsigned int yPeriodic,
                   unsigned int zPeriodic) {
    m_Periodic(0) = xPeriodic;
    m_Periodic(1) = yPeriodic;
    m_Periodic(2) = zPeriodic;
  }

  void setRadialPeriodic( unsigned int rPeriodic)
  {
    m_Periodic(1) = rPeriodic;
  }

//  virtual void ReadXMLInput(TICPP::HierarchicalDataNode& hdn);

  virtual void SetContactGhostRange(const realT bufferSize);

//  virtual void ResetSinglePartitionGlobalToLocalMap(PhysicalDomainT& domain);

//  virtual void SetPeriodicDomainBoundaryObjects(PhysicalDomainT& domain);
//  virtual void CorrectReferencePositionsForPeriodicBoundaries(
//      PhysicalDomainT& domain);

//  void CreateSinglePartitionGhostObjects(PhysicalDomainT& domain,
//      const bool contactActive, const int elementGhostingDepth);
//  void SetSinglePartitionGhostArrays(PhysicalDomainT& domain);

  int GetColor();

  const array<integer>& GetPartitions() const {
    return m_Partitions;
  }

  const array<integer>& GetCoords() const {
    return m_coords;
  }

  const R1Tensor& xMin() const {
    return m_min;
  }
  const R1Tensor& xMax() const {
    return m_max;
  }

  realT xMin(const int i) const {
    return m_min[i];
  }
  realT xMax(const int i) const {
    return m_max[i];
  }

public:
  array<int> m_Partitions; // number of partitions
  array<int> m_Periodic; // 1 = periodic
  array<int> m_coords; // ijk partition indexes

  R1Tensor m_min; // Minimum extent of partition dimensions (excluding ghost
                  // objects)
  R1Tensor m_max; // Maximum extent of partition dimensions (excluding ghost
                  // objects)

  R1Tensor m_xBoundingBoxMin, m_xBoundingBoxMax;

  array<real64> m_PartitionLocations[3]; // locations of partition boundaries

  R1Tensor m_blockSize; // Length of partition dimensions (excluding ghost
                        // objects)

  R1Tensor m_gridSize; // Total length of problem dimensions (excluding ghost
                       // objects)
  R1Tensor m_gridMin; // Minimum extent of problem dimensions (excluding ghost
                      // objects)
  R1Tensor m_gridMax; // Maximum extent of problem dimensions (excluding ghost
                      // objects)

//	struct PeriodicSet {
//		int m_dimension;
//		array<string> m_setNames;
//		void ReadXML(TICPP::HierarchicalDataNode& hdn) {
//			m_dimension = hdn.GetAttributeValue<int>("dimension");
//			m_setNames = hdn.GetStringVector("setnames");
//		}
//	};
//	std::vector<PeriodicSet> m_periodicSets;

  void AddNeighbors(const unsigned int idim, MPI_Comm& cartcomm,
                    int* ncoords);

//	std::map<array<int>, unsigned int> neighborCommPtrIndx;

//	virtual void WriteSiloDerived(SiloFile& siloFile);
//
//	virtual void ReadSiloDerived(const SiloFile& siloFile);

};
}
#endif /* SPATIALPARTITION_H_ */
