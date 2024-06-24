/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "NewGlobalNumbering.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetSurfaceFilter.h>

#include <vector>

namespace geos::ghosting
{

  struct Exchange
  {
    std::set<NodeGlbIdx> nodes;
    std::set<Edge> edges;
    std::set<Face> faces;
  };

  void to_json(json &j,
               const Exchange &v)
  {
    j = json{{"nodes", v.nodes},
             {"edges", v.edges},
             {"faces", v.faces}};
  }

  void from_json(const json &j,
                 Exchange &v)
  {
    v.nodes = j.at("nodes").get<std::set<NodeGlbIdx>>();
    v.edges = j.at("edges").get<std::set<Edge>>();
    v.faces = j.at("faces").get<std::set<Face>>();
  }

  array1d<std::uint8_t> convertExchange(Exchange const &exchange)
  {
    std::vector<std::uint8_t> const tmp = json::to_cbor(json(exchange));
    array1d<std::uint8_t> result;
    result.reserve(tmp.size());
    for (std::uint8_t const &t : tmp)
    {
      result.emplace_back(t);
    }
    return result;
  }

  Exchange convertExchange(array1d<std::uint8_t> const &input)
  {
    std::vector<std::uint8_t> const tmp(std::cbegin(input), std::cend(input));
    return json::from_cbor(tmp, false).template get<Exchange>();
  }

  /**
   * @brief Extract the cells at the boundary of the mesh.
   * @param mesh The vtk mesh.
   * @return The vtk cell ids.
   */
  std::set<vtkIdType> extractBoundaryCells(vtkSmartPointer<vtkDataSet> mesh)
  {
    // TODO Better handle the boundary information, forgetting about the 3d cells and simply handling the outside shell.
    auto f = vtkDataSetSurfaceFilter::New();
    f->PassThroughCellIdsOn();
    f->PassThroughPointIdsOff();
    f->FastModeOff();

    string const originalCellsKey = "ORIGINAL_CELLS";
    f->SetOriginalCellIdsName(originalCellsKey.c_str());
    auto boundaryMesh = vtkPolyData::New();
    f->UnstructuredGridExecute(mesh, boundaryMesh);
    vtkIdTypeArray const *originalCells = vtkIdTypeArray::FastDownCast(boundaryMesh->GetCellData()->GetArray(originalCellsKey.c_str()));

    std::set<vtkIdType> boundaryCellIdxs;
    for (auto i = 0; i < originalCells->GetNumberOfTuples(); ++i)
    {
      boundaryCellIdxs.insert(originalCells->GetValue(i));
    }

    return boundaryCellIdxs;
  }

  Face reorderFaceNodes(std::vector<NodeGlbIdx> const &nodes, bool &isFlipped, std::uint8_t &start) // TODO unit test
  {
    std::size_t const n = nodes.size();

    std::vector<NodeGlbIdx> f(nodes);

    auto const cit = std::min_element(std::cbegin(nodes), std::cend(nodes));
    auto const minIdx = std::distance(std::cbegin(nodes), cit);

    std::size_t const prevIdx = minIdx == 0 ? n - 1 : minIdx - 1;
    std::size_t const nextIdx = minIdx == intConv<int>(n) - 1 ? 0 : minIdx + 1;

    start = intConv<std::uint8_t>(minIdx);
    isFlipped = nodes[prevIdx] < nodes[nextIdx];

    std::size_t const pivot = isFlipped ? n - 1 - minIdx : minIdx;
    if (isFlipped)
    {
      std::reverse(std::begin(f), std::end(f));
    }
    std::rotate(std::begin(f), std::begin(f) + pivot, std::end(f));

    return f;
  }

  std::vector<NodeLocIdx> resetFaceNodes(std::vector<NodeLocIdx> const &nodes,
                                         bool const &isFlipped,
                                         std::uint8_t const &start)
  {
    std::vector<NodeLocIdx> result(nodes);
    if (isFlipped)
    {
      std::reverse(std::begin(result), std::end(result));
      std::uint8_t const pivot = std::size(result) - 1 - start;
      std::rotate(std::begin(result), std::begin(result) + pivot, std::end(result));
    }
    else
    {
      std::rotate(std::rbegin(result), std::rbegin(result) + start, std::rend(result));
    }
    return result;
  }

  Exchange buildSpecificData(vtkSmartPointer<vtkDataSet> mesh,
                             std::set<vtkIdType> const &cellIds)
  {
    vtkIdTypeArray *gids = vtkIdTypeArray::FastDownCast(mesh->GetPointData()->GetGlobalIds());
    Span<vtkIdType> const globalPtIds((vtkIdType *)gids->GetPointer(0), gids->GetNumberOfTuples());

    std::vector<Edge> tmpEdges;
    std::vector<Face> tmpFaces;

    // Pre-allocation of the temporary vectors.
    {
      std::size_t numEdges = 0;
      std::size_t numFaces = 0;
      for (vtkIdType const &c : cellIds)
      {
        vtkCell *cell = mesh->GetCell(c);
        numEdges += cell->GetNumberOfEdges();
        numFaces += cell->GetNumberOfFaces();
      }
      tmpEdges.reserve(numEdges);
      tmpFaces.reserve(numFaces);
    }

    // Filling the temporary.
    for (auto c = 0; c < mesh->GetNumberOfCells(); ++c)
    {
      vtkCell *cell = mesh->GetCell(c);

      for (auto e = 0; e < cell->GetNumberOfEdges(); ++e)
      {
        vtkCell *edge = cell->GetEdge(e);
        vtkIdType const nli0 = edge->GetPointId(0);
        vtkIdType const nli1 = edge->GetPointId(1);
        vtkIdType const ngi0 = globalPtIds[nli0];
        vtkIdType const ngi1 = globalPtIds[nli1];
        tmpEdges.emplace_back(std::minmax({NodeGlbIdx{ngi0}, NodeGlbIdx{ngi1}}));
      }

      for (auto f = 0; f < cell->GetNumberOfFaces(); ++f)
      {
        vtkCell *face = cell->GetFace(f);
        vtkIdList *pids = face->GetPointIds();
        std::vector<NodeGlbIdx> nodes(pids->GetNumberOfIds());
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
          vtkIdType const nli = face->GetPointId(i);
          vtkIdType const ngi = globalPtIds[nli];
          nodes[i] = NodeGlbIdx{ngi};
        }
        bool isFlipped;
        std::uint8_t start;
        tmpFaces.emplace_back(reorderFaceNodes(nodes, isFlipped, start));
      }
    }

    std::set<NodeGlbIdx> nodes;
    std::transform(std::cbegin(globalPtIds), std::cend(globalPtIds),
                   std::inserter(nodes, std::end(nodes)), [](vtkIdType const &id)
                   { return NodeGlbIdx{id}; });

    // Removing the duplicates by copying into a `std::set`.
    std::set<Edge> edges{tmpEdges.cbegin(), tmpEdges.cend()}; // SortedArray requires the input to be sorted already.
    std::set<Face> faces{tmpFaces.cbegin(), tmpFaces.cend()};

    return {std::move(nodes), std::move(edges), std::move(faces)};
  }

  Exchange buildFullData(vtkSmartPointer<vtkDataSet> mesh)
  {
    std::set<vtkIdType> cellIds;
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i)
    {
      cellIds.insert(cellIds.end(), i);
    }

    return buildSpecificData(mesh, cellIds);
  }

  Exchange buildExchangeData(vtkSmartPointer<vtkDataSet> mesh)
  {
    return buildSpecificData(mesh, extractBoundaryCells(mesh));
  }

  /**
   * @brief
   * @tparam GLB_IDX The global index type of the geometrical quantity considered (typically @c EdgeGlbIdx or @c FaceGlbIdx).
   * @param sizes The sizes of the intersection bucket (from the current MPI rank).
   * @param offsets The offsets for each intersection bucket (from the MPI_Scan process, i.e. the previous ranks).
   * @param curRank Current MPI rank.
   * @return The offset buckets, updated with the sizes of the current MPI rank.
   */
  template <typename GLB_IDX>
  std::map<std::set<MpiRank>, GLB_IDX> updateBucketOffsets(std::map<std::set<MpiRank>, localIndex> const &sizes,
                                                           std::map<std::set<MpiRank>, GLB_IDX> const &offsets,
                                                           MpiRank curRank)
  {
    std::map<std::set<MpiRank>, GLB_IDX> reducedOffsets;

    // We only keep the offsets that are still relevant to the current and higher ranks.
    // Note that `reducedOffsets` will be used by the _current_ rank too.
    // Therefore, we'll send information that may not be useful to higher ranks.
    // This is surely acceptable because the information is tiny.
    // We do this by looping through the buckets, finding the beggest rank, and then only adding the bucket and its offset if biggest rank is >= current rank
    for (auto const &[ranks, offset] : offsets)
    {
      MpiRank const maxConcerned = *std::max_element(ranks.cbegin(), ranks.cend());
      // TODO invert
      //    if( maxConcerned >= curRank )
      //    {
      //      reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, offset );
      //    }
      if (maxConcerned < curRank)
      {
        continue;
      }
      reducedOffsets.emplace_hint(reducedOffsets.end(), ranks, offset);
    }

    // Add the offsets associated to the new size buckets from the current rank.
    GLB_IDX nextOffset{0};
    // Loop through th buckets in either (edges, faces)
    for (auto const &[ranks, size] : sizes)
    {
      // See if this bucket was in the previously indexed data
      // if so, then we just need to update the nextOffset, which will be used for the next new data
      // if the data is new, then we add its offset to the set of indexed data in reducedOffsets
      auto const it = reducedOffsets.find(ranks);
      if (it == reducedOffsets.end())
      {
        reducedOffsets.emplace_hint(reducedOffsets.end(), ranks, nextOffset);
        nextOffset += GLB_IDX{size};
      }
      else
      {
        nextOffset = it->second + GLB_IDX{size}; // Define the new offset from the last
      }
    }

    // Add an extra entry based for the following rank
    // There are edge cases where a rank cannot start counting without this
    reducedOffsets.emplace_hint(reducedOffsets.end(), std::set<MpiRank>{curRank + 1_mpi}, nextOffset);

    return reducedOffsets;
  }

  // TODO Duplicated
  std::map<MpiRank, Exchange> exchange(int commId,
                                       std::vector<NeighborCommunicator> &neighbors,
                                       Exchange const &data)
  {
    MPI_iCommData commData(commId);
    integer const numNeighbors = LvArray::integerConversion<integer>(neighbors.size());
    commData.resize(numNeighbors);

    array1d<std::uint8_t> const cv = convertExchange(data);

    for (integer i = 0; i < numNeighbors; ++i)
    {
      neighbors[i].mpiISendReceiveSizes(cv,
                                        commData.mpiSendBufferSizeRequest(i),
                                        commData.mpiRecvBufferSizeRequest(i),
                                        commId,
                                        MPI_COMM_GEOSX);
    }

    MpiWrapper::waitAll(numNeighbors, commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus());
    MpiWrapper::waitAll(numNeighbors, commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus());

    array1d<array1d<std::uint8_t>> rawExchanged(neighbors.size());

    for (integer i = 0; i < numNeighbors; ++i)
    {
      neighbors[i].mpiISendReceiveData(cv,
                                       commData.mpiSendBufferRequest(i),
                                       rawExchanged[i],
                                       commData.mpiRecvBufferRequest(i),
                                       commId,
                                       MPI_COMM_GEOSX);
    }
    MpiWrapper::waitAll(numNeighbors, commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus());
    MpiWrapper::waitAll(numNeighbors, commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus());

    std::map<MpiRank, Exchange> output;
    for (auto i = 0; i < numNeighbors; ++i)
    {
      output[MpiRank{neighbors[i].neighborRank()}] = convertExchange(rawExchanged[i]);
    }
    return output;
  }

  /**
   * @brief
   * @tparam T Typically NodeGloIdx or a container of NodeGlbIdx (for edges and faces)
   * @param exchanged
   * @param curRank
   * @param neighbors
   * @return
   */
  template <typename T>
  std::map<std::set<MpiRank>, std::set<T>>
  buildIntersectionBuckets(std::map<MpiRank, std::set<T> const &> const &exchanged,
                           MpiRank curRank,
                           std::set<MpiRank> const &neighbors)
  {
    // counts - which ranks use each (node, edge, face)
    std::map<T, std::set<MpiRank>> counts; // TODO Use better intersection algorithms?

    // We "register" all the edges of the current rank: they are the only one we're interested in.
    // i.e. now exchanged is just one of (nodes, edges, faces), we add keys for each of the (nodes, edges, faces)
    // to counts and say it is used by the current rank
    for (T const &node : exchanged.at(curRank))
    {
      counts.emplace_hint(counts.end(), node, std::set<MpiRank>{curRank});
    }

    // We now loop on the neighbor edges.
    // If a neighbor has an edge in common with the current rank, they we store it.
    for (MpiRank const &neighborRank : neighbors) // This does not include the current rank.
    {
      // for each of the (nodes, edges, faces) on the boundary of the other ranks, 
      // add that that entry is also used by the neighbor rank
      for (T const &node : exchanged.at(neighborRank))
      {
        auto it = counts.find(node);
        if (it != counts.cend()) // TODO Extract `counts.cend()` out of the loop.
        {
          it->second.insert(neighborRank);
        }
      }
    }

    // Counts is essentially the inverse relationship to what we want
    // now invert so that we have keys as set of MPI ranks sharing and value as the (node, edge, face)
    std::map<std::set<MpiRank>, std::set<T>> nodeBuckets;
    for (auto const &[node, ranks] : counts)
    {
      if (ranks.find(curRank) != ranks.cend())
      {
        nodeBuckets[ranks].insert(node);
      }
    }
    return nodeBuckets;
  }

  /**
   * @brief Compute the intersection between for the ranks based on the information @p exchanged.
   * @param exchanged Geometrical information provided by the ranks (including the current rank).
   * @param curRank The current MPI rank.
   * @param neighbors excluding current rank
   * @return The intersection buckets.
   */
  Buckets buildIntersectionBuckets(std::map<MpiRank, Exchange> const &exchanged,
                                   MpiRank curRank,
                                   std::set<MpiRank> const &neighbors)
  {
    // create storage for the set of nodes owned by each rank, then fill using exchange info
    // this will be passed to function which computes intersections
    std::map<MpiRank, std::set<NodeGlbIdx> const &> nodeInfo;
    for (auto const &[rank, exchange] : exchanged)
    {
      nodeInfo.emplace(rank, exchange.nodes);
    }
    std::map<std::set<MpiRank>, std::set<NodeGlbIdx>> const nodeBuckets = buildIntersectionBuckets(nodeInfo, curRank, neighbors);

    // Repeat the process for the edges
    std::map<MpiRank, std::set<Edge> const &> edgeInfo;
    for (auto const &[rank, exchange] : exchanged)
    {
      edgeInfo.emplace(rank, exchange.edges);
    }
    std::map<std::set<MpiRank>, std::set<Edge>> const edgeBuckets = buildIntersectionBuckets(edgeInfo, curRank, neighbors);

    // For faces, the algorithm can be a tad simpler because faces can be shared by at most 2 ranks.
    std::map<std::set<MpiRank>, std::set<Face>> faceBuckets;
    std::set<Face> curFaces = exchanged.at(curRank).faces;
    for (MpiRank const &neighborRank : neighbors) // This does not include the current rank.
    {
      for (Face const &face : exchanged.at(neighborRank).faces)
      {
        auto it = curFaces.find(face);
        if (it != curFaces.cend())
        {
          faceBuckets[{curRank, neighborRank}].insert(*it);
          curFaces.erase(it);
        }
      }
    }
    faceBuckets[{curRank}] = curFaces;

    return {nodeBuckets, edgeBuckets, faceBuckets};
  }

  /**
   * @brief Estimate an upper bound of the serialized size of the size @p buckets.
   * @param buckets
   * @return Size in bytes.
   * @details @c MPI_Scan requires a fixed buffer size.
   * An a-priori estimation of the maximum size of the buffer leads to crazy amounts of memory
   * because of the combinatorial pattern of the rank intersections.
   * Therefore, we use an initial MPI communication to get a better estimation.
   * @node We don't care about the @c nodes because their global numbering is already done.
   */
  std::size_t buildMaxBufferSize(Buckets const &buckets) // TODO give the size buckets instead?
  {
    // this is the custom MPI reduction function we use in the sum across ranks
    auto const f = [](auto const &bucket) -> std::size_t
    {
      std::size_t size = std::size(bucket);
      for (auto const &[ranks, _] : bucket)
      {
        size += std::size(ranks) + 1 + 1; // One `+1` for the size of the ranks set, the other one for the offset
      }
      return size;
    };

    return MpiWrapper::sum(f(buckets.edges) + f(buckets.faces)); // Add the ScannedOffsets::{edge, face}Restart
  }

  struct BucketSizes
  {
    using mapping = std::map<std::set<MpiRank>, localIndex>;
    mapping edges;
    mapping faces;
  };

  void to_json(json &j,
               const BucketSizes &v)
  {
    j = json{{"edges", v.edges},
             {"faces", v.faces}};
  }

  void from_json(const json &j,
                 BucketSizes &v)
  {
    v.edges = j.at("edges").get<BucketSizes::mapping>();
    v.faces = j.at("faces").get<BucketSizes::mapping>();
  }

  std::vector<std::uint8_t> serialize(BucketSizes const &sizes)
  {
    return json::to_cbor(json(sizes));
  }

  std::vector<std::uint8_t> serialize(BucketOffsets const &offsets)
  {
    return json::to_cbor(json(offsets));
  }

  /**
   * @brief Compute the sizes of the intersection buckets.
   * @param buckets The intersection buckets.
   * @return The buckets of sizes.
   */
  BucketSizes getBucketSize(Buckets const &buckets)
  {
    BucketSizes output;

    for (auto const &[ranks, edges] : buckets.edges)
    {
      output.edges.emplace_hint(output.edges.end(), ranks, std::size(edges));
    }

    for (auto const &[ranks, faces] : buckets.faces)
    {
      output.faces.emplace_hint(output.faces.end(), ranks, std::size(faces));
    }

    return output;
  }

  /**
   * @brief
   * @tparam V Container of std::uint8_t
   * @param data
   * @return
   */
  template <class V>
  BucketOffsets deserialize(V const &data)
  {
    return json::from_cbor(data, false).template get<BucketOffsets>();
  }

  /**
   * @brief Custom MPI reduction function. Merges the sizes of the bucket from current rank into numbering offsets.
   * @param[in] in Contains the reduced result from the previous ranks. Needs to be unpacked into a @c BucketOffsets instance.
   * @param[inout] inout As an @e input, contains a pointer to the @c BucketSizes instance from the current rank.
   * As an @e output, will contained the (packed) instance of @c BucketOffsets after the @c reduction.
   * The @e output instance will be used by both current and next ranks.
   * @param[in] len Number of data.
   * @param[in] dataType Type of data. Mean to be @c MPI_BYTE for the current case.
   */
  void f(void *in,
         void *inout,
         int *len,
         MPI_Datatype *dataType)
  {
    GEOS_ASSERT_EQ(*dataType, MPI_BYTE);

    MpiRank const curRank{MpiWrapper::commRank()};

    // offsets provided by the previous rank(s)
    BucketOffsets const offsets = deserialize(Span<std::uint8_t>((std::uint8_t *)in, *len));

    // Sizes provided by the current rank, under the form of a pointer to the data.
    // No need to serialize, since we're on the same rank.
    std::uintptr_t addr;
    std::memcpy(&addr, inout, sizeof(std::uintptr_t));
    BucketSizes const *sizes = reinterpret_cast<BucketSizes const *>(addr);

    BucketOffsets updatedOffsets;
    updatedOffsets.edges = updateBucketOffsets<EdgeGlbIdx>(sizes->edges, offsets.edges, curRank);
    updatedOffsets.faces = updateBucketOffsets<FaceGlbIdx>(sizes->faces, offsets.faces, curRank);

    // Serialize the updated offsets, so they get sent to the next rank.
    std::vector<std::uint8_t> const serialized = serialize(updatedOffsets);
    std::memcpy(inout, serialized.data(), serialized.size());
  }

  std::tuple<Buckets, BucketOffsets> doTheNewGlobalNumbering(vtkSmartPointer<vtkDataSet> mesh,
                                                             std::set<MpiRank> const &neighbors)
  {
    // Now we exchange the data with our neighbors.
    MpiRank const curRank{MpiWrapper::commRank()};

    // ncs holds communicators to the neighboring ranks
    std::vector<NeighborCommunicator> ncs;
    ncs.reserve(neighbors.size());
    for (MpiRank const &rank : neighbors)
    {
      ncs.emplace_back(rank.get());
    }

    CommID const commId = CommunicationTools::getInstance().getCommID();

    // On the current rank, we have a map from the neighbor MPI rank to a corresponding set of exchanged data (sets of nodes, edges, faces)
    // buildExchangeData(mesh) just extracts the boundary (faces) of the current rank mesh
    // exchange() actually does the communication, so that the output exchanged maps from nhbr rank to their boundary data
    std::map<MpiRank, Exchange> exchanged = exchange(int(commId), ncs, buildExchangeData(mesh));

    // Then I guess we just put the full mesh data for the current rank into the corresponding slot of exchanged, i think this is just convenient
    exchanged[MpiRank{MpiWrapper::commRank()}] = buildFullData(mesh);

     //std::cout << "exchanged on rank " << curRank << " -> " << json( exchanged[MpiRank{ MpiWrapper::commRank() }] ) << std::endl;

    // Buckets - for each of (nodes, edges, faces), map from set of MPI ranks to the set of (nodes, edges, faces) shared by those ranks
    Buckets const buckets = buildIntersectionBuckets(exchanged, curRank, neighbors);

    // Next step is communication. 
    // We will be passing around info about of the each ranks' buckets, and to do this we need to know the size of the data to be sent
    // Technically this changes, as each rank has a different number of buckets
    // For MPI scan, we need a fixed send/recv size, so we estimate with a large buffer size
    // We dont want to do this too naively, as just taking the total possible cobinations of rank intersections leads to a massive size
    // Instead, we do an MPI operation here to estimate an appropriate size
    // The idea is to just sum the size of buckets across the ranks
    std::size_t const maxBufferSize = 100 * buildMaxBufferSize(buckets);

    std::vector<std::uint8_t> sendBuffer(maxBufferSize);
    std::vector<std::uint8_t> recvBuffer(maxBufferSize);

    // The idea is that each rank needs to know where to start indexing the (nodes, edges, faces) it owns
    // To do this, we need BucketOffsets, which are like buckets, but instead associate a global index to each set of MPI ranks, which marks the index of the first (node, edge, face) in that bucket
    // To get this, we need to know the size of each bucket on the current rank, stored in BucketSizes
    BucketSizes const sizes = getBucketSize(buckets);
    BucketOffsets offsets;
    if (curRank == 0_mpi)
    {
      // The `MPI_Scan` process will not call the reduction operator for rank 0.
      // So we need to reduce for ourselves.
      // updateBucketOffsets just takes the bucket sizes on the current rank,
      // the previously computed offsets from other MPI ranks (none in this case),
      // and then does the additions so that the global indices for the current buckets are known in offsets
      offsets.edges = updateBucketOffsets<EdgeGlbIdx>(sizes.edges, {{{
                                                                         0_mpi,
                                                                     },
                                                                     0_egi}},
                                                      curRank);
      offsets.faces = updateBucketOffsets<FaceGlbIdx>(sizes.faces, {{{
                                                                         0_mpi,
                                                                     },
                                                                     0_fgi}},
                                                      curRank);
      // Still we need to send this reduction to the following rank, by copying to it to the send buffer.
      std::vector<std::uint8_t> const bytes = serialize(offsets);
      std::memcpy(sendBuffer.data(), bytes.data(), bytes.size());
    }
    else
    {
      // For the other ranks, the reduction operator will be called during the `Mpi_Scan` process.
      // So unlike for rank 0, we do not have to do it ourselves.
      // In order to provide the `sizes` to the reduction operator,
      // since `sizes` will only be used on the current rank,
      // we'll provide the information as a pointer to the instance.
      // The reduction operator will then compute the new offsets and send them to the following rank.
      std::uintptr_t const addr = reinterpret_cast<std::uintptr_t>(&sizes);
      std::memcpy(sendBuffer.data(), &addr, sizeof(std::uintptr_t));
    }

    // f is our mpi reduction function, all it essentially does is call the updateBucketOffsets as above,
    // but it also handles the packing/unpacking of transferred data
    // Note that, right now, packing/unpacking is done via json, this should be changed to GEOS packing routines
    MPI_Op op;
    MPI_Op_create(f, false, &op);

    MPI_Scan(sendBuffer.data(), recvBuffer.data(), maxBufferSize, MPI_BYTE, op, MPI_COMM_WORLD);

    if (curRank != 0_mpi)
    {
      offsets = deserialize(recvBuffer);
    }

    // With the offsets for each bucket, we are done with global numbering, 
    // as each (edge, face) is numbered by starting from the offset for that bucket and counting up
    return {std::move(buckets), std::move(offsets)};
  }

} // end of namespace