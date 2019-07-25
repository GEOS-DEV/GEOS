/*
 * kernelMacros.hpp
 *
 *  Created on: Jul 17, 2019
 *      Author: settgast
 */

#ifndef CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELMACROS_HPP_
#define CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELMACROS_HPP_

#define UPDATE_STRESS 0

#define STANDARD_ELEMENT_DNDX_LAYOUT 0
#define STANDARD_ELEMENT_DETJ_LAYOUT 0
#define STANDARD_ELEMENT_MEANSTRESS_LAYOUT 0
#define STANDARD_ELEMENT_DEVIATORSTRESS_LAYOUT 0
#define STANDARD_ELEMENT_TONODESRELATION_LAYOUT 0

#define CALC_SHAPE_FUNCTION_DERIVATIVES 1

#define STORE_NODE_DATA_LOCALLY 0

#if STORE_NODE_DATA_LOCALLY
  #define VELOCITY_ACCESSOR(k, a, b) v_local[ a ][ b ]
  #define POSITION_ACCESSOR(k, a, b) x_local[ b ][ a ]
#else
  #define VELOCITY_ACCESSOR(k, a, b) vel[ TONODESRELATION_ACCESSOR(elemsToNodes, k, a ) ][ b ]
  #define POSITION_ACCESSOR(k, a, b) X[ TONODESRELATION_ACCESSOR(elemsToNodes, k, a ) ][ b ]
#endif

#if CALC_SHAPE_FUNCTION_DERIVATIVES
  #define DNDX_ACCESSOR(dNdX, k, q, n, i) dNdX[i][n]
#elif STANDARD_ELEMENT_DNDX_LAYOUT
  #define DNDX_ACCESSOR(dNdX, k, q, n, i) dNdX(k, q, n, i)
#else
  #define DNDX_ACCESSOR(dNdX, k, q, n, i) dNdX(i, n, q, k)
#endif

#if STANDARD_ELEMENT_DETJ_LAYOUT
  #define DETJ_ACCESSOR(detJ, k, q) detJ(k, q)
#else
  #define DETJ_ACCESSOR(detJ, k, q) detJ(q, k)
#endif

#if STANDARD_ELEMENT_MEANSTRESS_LAYOUT
  #define MEANSTRESS_ACCESSOR(meanStress, k, q) meanStress(k, q)
#else
  #define MEANSTRESS_ACCESSOR(meanStress, k, q) meanStress(q, k)
#endif

#if STANDARD_ELEMENT_DEVIATORSTRESS_LAYOUT
  #define DEVIATORSTRESS_ACCESSOR(deviatorStress, k, q, i) deviatorStress(k, q, i)
#else
  #define DEVIATORSTRESS_ACCESSOR(deviatorStress, k, q, i) deviatorStress(i, q, k)
#endif

#if STANDARD_ELEMENT_TONODESRELATION_LAYOUT
  #define TONODESRELATION_ACCESSOR(toNodes, k, i) toNodes(k, i)
#else
  #define TONODESRELATION_ACCESSOR(toNodes, k, i) toNodes(i, k)
#endif



#endif /* CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELMACROS_HPP_ */
