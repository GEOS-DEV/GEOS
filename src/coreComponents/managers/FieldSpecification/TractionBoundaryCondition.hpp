/*
 * TractionBoundaryCondition.h
 *
 *  Created on: Mar 15, 2021
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_MANAGERS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_H_
#define SRC_CORECOMPONENTS_MANAGERS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_H_

#include "FieldSpecificationBase.hpp"
#include "mesh/FaceManager.hpp"

namespace geosx
{

class TractionBoundaryCondition : public FieldSpecificationBase
{
public:
  TractionBoundaryCondition( string const & name, Group * parent );
  TractionBoundaryCondition() = delete;
  virtual ~TractionBoundaryCondition() = default;
//  TractionBoundaryCondition( const TractionBoundaryCondition &other ) = default;
//  TractionBoundaryCondition( TractionBoundaryCondition &&other ) = default;
//  TractionBoundaryCondition& operator=( const TractionBoundaryCondition &other ) = default;
//  TractionBoundaryCondition& operator=( TractionBoundaryCondition &&other ) = default;


  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Traction"; }

  void launch( real64 const time,
               arrayView1d< globalIndex const > const blockLocalDofNumber,
               globalIndex const dofRankOffset,
               FaceManager const & faceManager,
               SortedArrayView< localIndex const > const & targetSet,
               arrayView1d< real64 > const & localRhs ) const ;

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    constexpr static char const * tractionTypeString() { return "tractionType"; }
    constexpr static char const * inputStressString() { return "inputStress"; }



  };

protected:
  virtual void postProcessInput() override final;

  virtual void initializePreSubGroups() override final;

//  template< int FUNCTION_TYPE >
//  void evaluateMagnitude( FaceManager const & faceManager, real64 const time,  )


  int m_tractionType;

  VoigtTensor m_inputStress;

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_MANAGERS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_H_ */
