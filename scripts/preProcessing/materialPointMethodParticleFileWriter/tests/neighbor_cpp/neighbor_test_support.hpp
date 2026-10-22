// Serial, bounds-checked framework adapters for the production neighbor routine.
// These adapters exercise control flow and list contents, not MPI/GPU execution.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>
using real64 = double;
using localIndex = int;
#define GEOS_MARK_FUNCTION
#define GEOS_HOST_DEVICE
#define GEOS_HOST
#define GEOS_LOG_RANK(...)
#define GEOS_ERROR_IF(condition, ...) do { if(condition) throw std::runtime_error("neighbor invariant failed"); } while(false)
constexpr int parallelDeviceMemorySpace = 1, hostMemorySpace = 0;
struct serialPolicy {};
template<class...> struct parallelDevicePolicy {};
struct parallelDeviceAtomic {};
struct parallelDeviceReduce {};
struct parallelDeviceMultiReduce {};
inline void parallelDeviceSync() {}
inline double MPI_Wtime() { return 0.0; }
template<class Policy, class F> void forAll(int n, F f) { for(int p = 0; p < n; ++p) f(p); }
template<class T, class U> T dynamicCast(U &v) { return static_cast<T>(v); }

template<class T> struct Vector {
  std::shared_ptr<std::vector<T>> data = std::make_shared<std::vector<T>>();
  Vector() = default;
  explicit Vector(int n): data(std::make_shared<std::vector<T>>(n)) {}
  int size() const { return static_cast<int>(data->size()); }
  T &operator[](int i) const { return data->at(i); }
  void push_back(T value) { data->push_back(value); }
  void move(int, bool = false) const {}
  Vector toView() const { return *this; }
  Vector toViewConst() const { return *this; }
};
template<class T> struct Matrix {
  std::shared_ptr<std::vector<std::vector<T>>> data;
  Matrix(): data(std::make_shared<std::vector<std::vector<T>>>()) {}
  Matrix(int n, int width): data(std::make_shared<std::vector<std::vector<T>>>(n, std::vector<T>(width))) {}
  std::vector<T> &operator[](int i) const { return data->at(i); }
  void move(int, bool = false) const {}
  Matrix toView() const { return *this; }
  Matrix toViewConst() const { return *this; }
};
template<class T> using array1d = Vector<T>;
template<class T> using arrayView1d = Vector<std::remove_const_t<T>>;
template<class T> using array2d = Matrix<T>;
template<class T> using arrayView2d = Matrix<std::remove_const_t<T>>;
template<class T> using SortedArrayView = Vector<std::remove_const_t<T>>;
template<class T> struct Ragged {
  std::shared_ptr<std::vector<std::vector<T>>> rows;
  Ragged(): rows(std::make_shared<std::vector<std::vector<T>>>()) {}
  Ragged(int n, int width): rows(std::make_shared<std::vector<std::vector<T>>>(n, std::vector<T>(width))) {}
  int sizeOfArray(int i) const { return static_cast<int>(rows->at(i).size()); }
  void resizeArray(int i, int n) { rows->at(i).resize(n); }
  struct Row {
    std::vector<T> &data;
    T &operator[](int i) const { return data.at(i); }
  };
  Row operator[](int i) const { return {rows->at(i)}; }
  void move(int, bool = false) const {}
  Ragged toView() const { return *this; }
};
template<class T> using ArrayOfArrays = Ragged<T>;
template<class T> using ArrayOfArraysView = Ragged<T>;
namespace RAJA {
template<class Policy, class T> struct ReduceSum {
  std::shared_ptr<T> value;
  explicit ReduceSum(T v): value(std::make_shared<T>(v)) {}
  void operator+=(T v) const { *value += v; }
  T get() const { return *value; }
};
template<class Policy, class T> struct MultiReduceSum {
  Vector<T> values;
  MultiReduceSum(int n, T v): values(n) { std::fill(values.data->begin(), values.data->end(), v); }
  struct Ref { T &value; void operator+=(T v) const { value += v; } T get() const { return value; } };
  Ref operator[](int i) const { return {values[i]}; }
};
template<class Policy, class T> T atomicAdd(Policy, T *p, T n) { T old = *p; *p += n; return old; }
template<class Policy, class T> void atomicMax(Policy, T *p, T n) { *p = std::max(*p, n); }
}
namespace LvArray {
namespace math { using std::floor; using std::ceil; using std::max; using std::min; }
namespace tensorOps {
template<int N, class A, class B> void copy(A &&a, B const &b) { for(int d = 0; d < N; ++d) a[d] = b[d]; }
template<int N, class A, class B> void subtract(A &&a, B const &b) { for(int d = 0; d < N; ++d) a[d] -= b[d]; }
template<int N, class A> double l2NormSquared(A const &a) { double r = 0; for(int d = 0; d < N; ++d) r += a[d]*a[d]; return r; }
}
}
namespace fields { namespace mpm { struct particleDeleteFlag {}; } }
struct OrderedVariableToManyParticleRelation {
  Vector<int> m_numParticles;
  Ragged<int> m_toParticleRegion, m_toParticleSubRegion, m_toParticleIndex;
  void freeOnDevice() {}
  void resizeFromCapacities(Vector<int> const &n) {
    m_numParticles = Vector<int>(n.size());
    m_toParticleRegion = Ragged<int>(n.size(), 0);
    m_toParticleSubRegion = Ragged<int>(n.size(), 0);
    m_toParticleIndex = Ragged<int>(n.size(), 0);
    for(int i = 0; i < n.size(); ++i) {
      m_numParticles[i] = n[i];
      m_toParticleRegion.resizeArray(i, n[i]);
      m_toParticleSubRegion.resizeArray(i, n[i]);
      m_toParticleIndex.resizeArray(i, n[i]);
    }
  }
};
struct ParticleRegion { int id; int getIndexInParent() const { return id; } };
struct Parent { ParticleRegion region; ParticleRegion &getParent() { return region; } };
struct ParticleSubRegion {
  Parent parent;
  int id;
  Matrix<double> positions;
  Vector<int> deleted, active;
  OrderedVariableToManyParticleRelation neighbors;
  ParticleSubRegion(int region, int subregion, int n): parent{{region}}, id(subregion), positions(n, 3), deleted(n) {}
  int size() const { return deleted.size(); }
  Matrix<double> getParticleCenter() const { return positions; }
  template<class Field> Vector<int> getField() const { return deleted; }
  Parent &getParent() { return parent; }
  int getIndexInParent() const { return id; }
  Vector<int> activeParticleIndices() const { return active; }
  OrderedVariableToManyParticleRelation &neighborList() { return neighbors; }
};
struct ParticleManager {
  std::vector<std::unique_ptr<ParticleSubRegion>> regions;
  int getNumberOfParticles() const { int n = 0; for(auto const &r: regions) n += r->size(); return n; }
  template<class F> void forParticleSubRegions(F f) { for(auto &r: regions) f(*r); }
};
struct SolidMechanicsMPM {
  int m_numberOfSubRegions = 3, m_planeStrain = 1, m_numDims = 2, m_binSizeMultiplier = 1;
  double m_neighborRadius = 0.25;
  std::array<double, 3> m_xLocalMinNoGhost{{0, 0, 0}}, m_xLocalMaxNoGhost{{1, 1, 0}};
  double computeNeighborList(ParticleManager &);
};
