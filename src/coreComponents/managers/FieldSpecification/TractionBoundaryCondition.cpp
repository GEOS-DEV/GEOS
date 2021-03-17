/*
 * TractionBoundaryCondition.cpp
 *
 *  Created on: Mar 15, 2021
 *      Author: settgast
 */

#include "TractionBoundaryCondition.hpp"

namespace geosx
{

using namespace dataRepository;

TractionBoundaryCondition::TractionBoundaryCondition( string const & name, Group * parent ):
  FieldSpecificationBase( name, parent ),
  m_tractionType(0),
  m_inputStress{}
{
  registerWrapper( viewKeyStruct::tractionTypeString(), &m_tractionType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of traction boundary condition. Options are:\n"
                    "0 - Traction is applied to the faces as specified from the scale and direction,\n"
                    "1 - Traction is applied to the faces as a pressure specified from the product of scale and the outward face normal, \n"
                    "2 - Traction is applied to the faces as specified by the inner product of input stress and face normal." );

  registerWrapper( viewKeyStruct::inputStressString(), &m_inputStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Input stress for tractionType option 2:\n" );


  getWrapper<string>( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName(catalogName());

  getWrapper<int>( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );


}


void TractionBoundaryCondition::postProcessInput()
{

}

void TractionBoundaryCondition::initializePreSubGroups()
{
  string const & functionName = this->getFunctionName();

}


void TractionBoundaryCondition::launch( real64 const time,
                                        arrayView1d< globalIndex const > const blockLocalDofNumber,
                                        globalIndex const dofRankOffset,
                                        FaceManager const & faceManager,
                                        SortedArrayView< localIndex const > const & targetSet,
                                        arrayView1d< real64 > const & localRhs ) const
{


  arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();
  FunctionManager const & functionManager = getGlobalState().getFunctionManager();

  string const & functionName = this->getFunctionName();

  globalIndex_array nodeDOF;
  real64_array nodeRHS;

  R1Tensor const direction = this->getDirection(0);
  int const tractionType = m_tractionType;
  Tensor<real64, 6> inputStress = m_inputStress;

  real64 tractionMagnitude0;
  array1d< real64 > resultsArray( targetSet.size() );
  bool spatialFunction = false;


  if( functionName.empty() )
  {
    tractionMagnitude0 = this->getScale();
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
    if( function.isFunctionOfTime() == 2 )
    {
      tractionMagnitude0 = this->getScale() * function.evaluate( &time );
    }
    else
    {
      resultsArray.setName( "SolidMechanicsLagrangianFEM::TractionBC function results" );
      function.evaluate( faceManager, time, targetSet, resultsArray );
      spatialFunction = true;
      tractionMagnitude0 = 1e99;
    }
  }

  arrayView1d< real64 const > const results = resultsArray;

  {
    forAll< serialPolicy >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const kf = targetSet[ i ];
      localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );

      // TODO consider dispatch if appropriate
      real64 const tractionMagnitude = spatialFunction ? results[i] : tractionMagnitude0;

      real64 traction[3];
      // TODO consider dispatch if appropriate
      if( tractionType==0)
      {
        traction[0] = tractionMagnitude * direction[0];
        traction[1] = tractionMagnitude * direction[1];
        traction[2] = tractionMagnitude * direction[2];
      }
      else if( tractionType==1 )
      {
        traction[0] = tractionMagnitude * faceNormal(kf,0);
        traction[1] = tractionMagnitude * faceNormal(kf,1);
        traction[2] = tractionMagnitude * faceNormal(kf,2);
      }
      else if( tractionType==2 )
      {
        traction[0] = tractionMagnitude * ( inputStress[0] * faceNormal(kf,0) + inputStress[5] * faceNormal(kf,1) + inputStress[4] * faceNormal(kf,2) );
        traction[1] = tractionMagnitude * ( inputStress[5] * faceNormal(kf,0) + inputStress[1] * faceNormal(kf,1) + inputStress[3] * faceNormal(kf,2) );
        traction[2] = tractionMagnitude * ( inputStress[4] * faceNormal(kf,0) + inputStress[3] * faceNormal(kf,1) + inputStress[2] * faceNormal(kf,2) );
      }

      // TODO replace with proper FEM integration
      traction[0] *= faceArea[kf] / numNodes;
      traction[1] *= faceArea[kf] / numNodes;
      traction[2] *= faceArea[kf] / numNodes;

      for( localIndex a=0; a<numNodes; ++a )
      {
        localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] - dofRankOffset;
        if( dof < 0 || dof >= localRhs.size() )
          continue;
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+0], traction[0] );
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+1], traction[1] );
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+2], traction[2] );
      }
    } );
  }
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, TractionBoundaryCondition, string const &, Group * const )


} /* namespace geosx */
