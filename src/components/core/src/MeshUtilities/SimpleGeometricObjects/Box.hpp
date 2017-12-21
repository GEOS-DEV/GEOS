/*
 * Box.hpp
 *
 *  Created on: Aug 4, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class Box : public SimpleGeometricObjectBase
{
public:
  Box( const std::string& name,
       ManagedGroup * const parent );

  virtual ~Box();

  static string CatalogName() { return "Box"; }

  virtual void FillDocumentationNode() override;

  virtual void ReadXML_PostProcess() override final;

  bool IsCoordInObject( const R1Tensor& coord ) const override final;
private:
  R1Tensor m_min;
  R1Tensor m_max;
  realT m_strikeAngle=0.0;
  R1Tensor m_boxCenter={0.0,0.0,0.0};
  realT m_cosStrike=0.0, m_sinStrike=0.0;

  struct viewKeyStruct
  {
    dataRepository::ViewKey xmin = { "xmin" };
    dataRepository::ViewKey xmax = { "xmax" };
  } viewKeys;


};
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
        */
