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
  
  /// Apply the time source Ricker at the location specified in the xml 
  virtual
  void applyRickerSource( real64 const time, DomainPartition & domain );
  
  /**@}*/


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr auto varName = "varName";

    static constexpr auto sourceCoordinatesString = "sourceCoordinates";
    static constexpr auto sourceNodeIdsString     = "sourceNodeIds";
    static constexpr auto sourceConstantsString   = "sourceConstants";           
    static constexpr auto sourceIsLocalString     = "sourceIsLocal";

    static constexpr auto timeSourceFrequencyString = "timeSourceFrequency";
    
  } waveEquationViewKeys;


  arrayView1d< string const > solidMaterialNames() const { return m_solidMaterialNames; }


protected:

  virtual void postProcessInput() override final;
  
  virtual void initializePostInitialConditionsPreSubGroups( dataRepository::Group * const problemManager ) override final;

  array1d < string const > m_solidMaterialNames;
  
private:

  /// Locates the source term and precomputes the constant part of the source term
  void precomputeSourceTerm( MeshLevel & mesh );

  /// Multiply the precomputed term by the ricker and add to the right-hand side
  void addSourceToRightHandSide( real64 const & time, arrayView1d< real64 > const rhs );
  
  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Indices of the nodes (in the right order) for each source point
  array2d< localIndex > m_sourceNodeIds;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants;
  
  /// Flag that indicates whether the source is local or not
  array1d< localIndex > m_sourceIsLocal;

  /// Central frequency for the Ricker time source
  real64 m_timeSourceFrequency;
  
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
