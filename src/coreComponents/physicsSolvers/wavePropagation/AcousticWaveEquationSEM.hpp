/*
 * WaveEquation.hpp
 *
 *  Created on: Jan 12, 2021
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{

class AcousticWaveEquationSEM : public SolverBase
{
public:
  AcousticWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~AcousticWaveEquationSEM() override;

  AcousticWaveEquationSEM() = delete;
  AcousticWaveEquationSEM( AcousticWaveEquationSEM const & ) = delete;
  AcousticWaveEquationSEM( AcousticWaveEquationSEM && ) = default;

  AcousticWaveEquationSEM & operator=( AcousticWaveEquationSEM const & ) = delete;
  AcousticWaveEquationSEM & operator=( AcousticWaveEquationSEM && ) = delete;


  static string catalogName() { return "AcousticSEM"; }

  virtual void initializePreSubGroups( Group * const rootGroup ) override;

  virtual void registerDataOnMesh( Group * const MeshBody ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

  /// Returns the value of a Ricker at time t0 with central Fourier frequency f0
  virtual
  real64 evaluateRicker( real64 const & t0, real64 const & f0 );

  /**@}*/


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr auto varName = "varName";

    static constexpr auto sourceCoordinatesString = "sourceCoordinates";
    static constexpr auto sourceNodeIdsString     = "sourceNodeIds";
    static constexpr auto sourceConstantsString   = "sourceConstants";
    static constexpr auto sourceIsLocalString     = "sourceIsLocal";

    static constexpr auto timeSourceFrequencyString = "timeSourceFrequency";

    static constexpr auto receiverCoordinatesString = "receiverCoordinates";
    static constexpr auto receiverNodeIdsString     = "receiverNodeIds";
    static constexpr auto receiverConstantsString   = "receiverConstants";
    static constexpr auto receiverIsLocalString     = "receiverIsLocal";

    static constexpr auto pressureNp1AtReceiversString   = "pressureNp1AtReceivers";


  } waveEquationViewKeys;


protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups( dataRepository::Group * const problemManager ) override final;

private:

  /**
   * @brief Convert a mesh element point coordinate into a coorinate on the reference element
   * @param coords coordinate of the point
   * @param coordsOnRefElem to contain the coordinate computed in the reference element
   * @param numElem index of the element containing the coords
   * @param numNodesPerElem number of nodes per element
   * @param elemsToNodes map obtain node global index from element index
   */
  //void computeCoordinateOnReferenceElement( array1d< real64 const > const & coords,
  //array1d< real64 > & coordsOnRefElem,
  //					     localIndex const & indexElem,
  //					     localIndex const & numNodesPerElem,
  //					     arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes);

  /// Locates the source term and precomputes the constant part of the source term
  /// And locate receivers and pre_evaluate the basis functions at each receiver coordinate
  void precomputeSourceAndReceiverTerm( MeshLevel & mesh );

  /// Multiply the precomputed term by the ricker and add to the right-hand side
  void addSourceToRightHandSide( real64 const & time, arrayView1d< real64 > const rhs );

  /// Compute the pressure at each receiver coordinate in one time step
  void computeSismoTrace( localIndex const num_timestep, arrayView1d< real64 > const pressure_np1 );

  /// save the sismo trace in file
  void saveSismo( localIndex isismo, real64 val_pressure, char *filename );

  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Indices of the nodes (in the right order) for each source point
  array2d< localIndex > m_sourceNodeIds;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants;

  /// Flag that indicates whether the source is local or not to the MPI rank
  array1d< localIndex > m_sourceIsLocal;

  /// Central frequency for the Ricker time source
  real64 m_timeSourceFrequency;

  /// Coordinates of the receivers in the mesh
  array2d< real64 > m_receiverCoordinates;

  /// Indices of the element nodes (in the right order) for each receiver point
  array2d< localIndex > m_receiverNodeIds;

  /// Basis function evaluated at the receiver for the nodes listed in m_receiverNodeIds
  array2d< real64 > m_receiverConstants;

  /// Flag that indicates whether the receiver is local or not to the MPI rank
  array1d< localIndex > m_receiverIsLocal;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array1d< real64 > m_pressureNp1AtReceivers;



};

namespace extrinsicMeshData
{

EXTRINSIC_MESH_DATA_TRAIT( Pressure_nm1,
                           "pressure_nm1",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Scalar pressure at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_n,
                           "pressure_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Scalar pressure at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_np1,
                           "pressure_np1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS,
                           "rhs",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "RHS" );

EXTRINSIC_MESH_DATA_TRAIT( MassVector,
                           "massVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Mass Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector,
                           "dampingVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Damping Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocity,
                           "mediumVelocity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Medium velocity of the cell" );

EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector,
                           "stiffnessVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Stiffness vector contains R_h*Pressure_n." );

}


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_ */
