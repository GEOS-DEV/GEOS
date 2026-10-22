// Compare every returned region/subregion/particle tuple with an independent
// all-pairs radius search. The production routine includes self-neighbors.
void checkAllPairs(SolidMechanicsMPM const &solver, ParticleManager const &pm)
{
  using Neighbor = std::tuple<int, int, int>;
  for(auto const &region: pm.regions) {
    auto const &r = *region;
    for(int pp = 0; pp < r.active.size(); ++pp) {
      int const a = r.active[pp];
      std::vector<Neighbor> expected, actual;
      for(auto const &other: pm.regions) for(int b = 0; b < other->size(); ++b) {
        if(other->deleted[b]) continue;
        double r2 = 0;
        for(int d = 0; d < solver.m_numDims; ++d) {
          double const delta = r.positions[a][d] - other->positions[b][d];
          r2 += delta * delta;
        }
        if(r2 <= solver.m_neighborRadius * solver.m_neighborRadius)
          expected.emplace_back(other->parent.region.id, other->id, b);
      }
      for(int n = 0; n < r.neighbors.m_numParticles[pp]; ++n)
        actual.emplace_back(r.neighbors.m_toParticleRegion[pp][n], r.neighbors.m_toParticleSubRegion[pp][n], r.neighbors.m_toParticleIndex[pp][n]);
      std::sort(actual.begin(), actual.end());
      std::sort(expected.begin(), expected.end());
      assert(actual == expected);
    }
  }
}
int main()
{
  int cases = 0;
  // Minimal reproduction of the reported failure: planar interior particles
  // with no boundary-roundoff ambiguity. The old routine visits z bin 1 even
  // though only bin 0 exists, eventually reading beyond the final subregion.
  {
    SolidMechanicsMPM solver;
    ParticleManager pm;
    for(int sub = 0; sub < 3; ++sub) {
      pm.regions.push_back(std::make_unique<ParticleSubRegion>(2, sub, sub == 2 ? 0 : 1));
      if(sub < 2) {
        auto &r = *pm.regions.back();
        r.positions[0][0] = 0.5; r.positions[0][1] = 0.5;
        r.active.push_back(0);
      }
    }
    solver.computeNeighborList(pm); checkAllPairs(solver, pm); ++cases;
  }
  for(int dimensions: {2, 3}) for(int multiplier: {1, 3}) for(unsigned seed: {1u, 7u}) {
    SolidMechanicsMPM solver;
    solver.m_numDims = dimensions; solver.m_planeStrain = dimensions == 2;
    solver.m_xLocalMaxNoGhost[2] = dimensions == 2 ? 0.0 : 1.0;
    solver.m_binSizeMultiplier = multiplier;
    ParticleManager pm;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(0, 1);
    // Two populated subregions exercise cross-subregion neighbor references;
    // the empty third subregion must not change list offsets or capacities.
    for(int sub = 0; sub < 3; ++sub) {
      int const count = sub == 2 ? 0 : 100;
      pm.regions.push_back(std::make_unique<ParticleSubRegion>(sub == 2 ? 6 : 2, sub, count));
      auto &r = *pm.regions.back();
      for(int p = 0; p < count; ++p) {
        for(int d = 0; d < dimensions; ++d) r.positions[p][d] = distribution(generator);
        r.deleted[p] = p > 20 && p % 11 == 0;
        // Noncontiguous active indices; other nondeleted particles are ghosts.
        if(!r.deleted[p] && p % 5 != 4) r.active.push_back(p);
      }
      if(count) {
        for(int d = 0; d < dimensions; ++d) {
          r.positions[0][d] = 0; r.positions[1][d] = 1;
          r.positions[2][d] = 0.5; r.positions[3][d] = 0.5;
          r.positions[4][d] = 0.5; r.positions[9][d] = 0.5;
          r.positions[14][d] = 0.5; r.positions[19][d] = 0.5;
        }
        // Exact distance cutoff, owned boundary, halo upper/lower boundary,
        // and points one representable value past the halo boundary.
        r.positions[2][0] = 0; r.positions[3][0] = 1;
        r.positions[4][0] = -solver.m_neighborRadius;
        r.positions[9][0] = 1 + solver.m_neighborRadius;
        r.positions[14][0] = std::nextafter(-solver.m_neighborRadius, -1.0);
        r.positions[19][0] = std::nextafter(1 + solver.m_neighborRadius, 2.0);
      }
    }
    solver.computeNeighborList(pm); checkAllPairs(solver, pm); ++cases;
    // Rebuilding discards old capacities. Inactive z storage cannot change XY
    // neighbor membership in the native-2D path.
    for(auto &r: pm.regions) for(int p = 0; p < r->size(); ++p) {
      if(dimensions == 2) r->positions[p][2] = 17.0 * p;
      r->deleted[p] = p % 7 == 0;
    }
    for(auto &r: pm.regions) {
      r->active = Vector<int>();
      for(int p = 0; p < r->size(); ++p) if(!r->deleted[p] && p % 5 != 4) r->active.push_back(p);
    }
    solver.computeNeighborList(pm); checkAllPairs(solver, pm); ++cases;
    // A rank with no owned particles still bins ghost candidates safely.
    for(auto &r: pm.regions) r->active = Vector<int>();
    solver.computeNeighborList(pm); checkAllPairs(solver, pm); ++cases;
  }
  std::cout << "Full production computeNeighborList: " << cases
            << " 2D/3D cases match all-pairs neighbors (mixed/empty subregions, ghosts, deletions, boundary points, rebuilds).\n";
}
