/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */



#include "EmbeddedSurfaceBlock.hpp"

namespace geos
{ 

localIndex EmbeddedSurfaceBlock::numEmbeddedSurfElem() const{
    return m_numEmbeddedSurfaces;
}

ArrayOfArrays<localIndex> EmbeddedSurfaceBlock::getEmbeddedSurfElemToNodes() const {
    return m_embeddedSurfElemToNodes;
}

ArrayOfArrays<localIndex> EmbeddedSurfaceBlock::getEmbeddedSurfElemTo3dElem() const {
     return m_embeddedSurfElemTo3dElem;
}

ArrayOfArrays<real64> EmbeddedSurfaceBlock::getEmbeddedSurfElemNodes() const {
   return m_embeddedSurfElemNodes;
}

ArrayOfArrays<real64> EmbeddedSurfaceBlock::getEmbeddedSurfElemNormalVectors() const {
    return m_embeddedSurfElemNormals;
}
ArrayOfArrays<real64> EmbeddedSurfaceBlock::getEmbeddedSurfElemTangentialLengthVectors() const{
    return m_embeddedSurfElemLengthVectors;
}
ArrayOfArrays<real64> EmbeddedSurfaceBlock::getEmbeddedSurfElemTangentialWidthVectors() const{
    return m_embeddedSurfElemWidthVectors;
}

void EmbeddedSurfaceBlock::setEmbeddedSurfElemNormalVectors(ArrayOfArrays<real64>&& _normals){
    m_embeddedSurfElemNormals = _normals;
}
void EmbeddedSurfaceBlock::setEmbeddedSurfElemTangentialLengthVectors(ArrayOfArrays<real64> && _lengthVectors){
    m_embeddedSurfElemLengthVectors = _lengthVectors;
}
void EmbeddedSurfaceBlock::setEmbeddedSurfElemTangentialWidthVectors(ArrayOfArrays<real64> && _widthVectors){ 
    m_embeddedSurfElemWidthVectors= _widthVectors;
}

void EmbeddedSurfaceBlock::setNumEmbeddedSurfElem(localIndex _numEmbeddedSurfaces){
    
    m_numEmbeddedSurfaces = _numEmbeddedSurfaces;
}

void EmbeddedSurfaceBlock::setEmbeddedSurfElemToNodes(ArrayOfArrays<localIndex> && _embeddedSurfElemToNodes){
    m_embeddedSurfElemToNodes = _embeddedSurfElemToNodes;
}

void EmbeddedSurfaceBlock::setEmbeddedSurfElemTo3dElem(ArrayOfArrays<localIndex> && _embeddedSurfElemTo3dElem){
    m_embeddedSurfElemTo3dElem = _embeddedSurfElemTo3dElem;
}

void EmbeddedSurfaceBlock::setEmbeddedSurfElemNodes(ArrayOfArrays<real64> && _embeddedSurfElemNodes){
    m_embeddedSurfElemNodes = _embeddedSurfElemNodes;
}

}