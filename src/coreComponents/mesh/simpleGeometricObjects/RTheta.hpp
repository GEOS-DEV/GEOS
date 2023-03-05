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

/**
 * @file RTheta.hpp
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_RTHETA_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_RTHETA_HPP_

#include "SimpleGeometricObjectBase.hpp"
#include "BoundedPlanarObject.hpp"

namespace geosx
{

/**
 * @class RTheta
 * @brief Class to represent a geometric disc in GEOSX.
 */
class RTheta : public BoundedPlanarObject
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  RTheta( const string & name,
                Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~RTheta() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "RTheta"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the center of the RTheta.
   * @return the center of the RTheta
   */
  virtual R1Tensor & getCenter() override final {return m_center;}

  /**
   * @copydoc getCenter()
   */
  virtual R1Tensor const & getCenter() const override final {return m_center;}

  /**
   * @brief Update the geometric function describing the boundary of the object.
   * @param coefficients define all the coefficients of the radius function.
   * @return void
   */
  void setRThetaFunction(array1d<real64> & coefficients)
  {
    m_radius.m_coefficients = coefficients;
  }

  /**
   * @brief Update a single coeff. of the geometric function describing the boundary of the object.
   * @param entry define the entry to be modified
   * @param value the new value for that entry
   * @return void
   */
  void setRThetaFunction(integer const entry, real64 const value)
  {
    if(m_radius.m_coefficients.size() < entry+1){
      std::cout<<"ERROR - COEFF VECTOR NOT LARGE ENOUGH\n";
      return;
      //GEOSX_ERROR("Vector of coefficients in the RTheta object is not large enough. Requested change to entry "<<entry<<" array size: "<<m_radius.m_coefficients.size()<<"\n");
    }
    m_radius.m_coefficients[entry] = value;
  }

  real64 getRadius(real64 angle) const
  {
    return m_radius.getRadius(angle);
  }



  class VariableRadius
  {  
    public: 

    VariableRadius()
    {
      m_coefficients.resize(6);
      m_coefficients[0] = 0.0;
      m_coefficients[1] = 0.0;
      m_coefficients[2] = 1.0;
      m_coefficients[3] = 0.0;
      m_coefficients[4] = 0.0;
      m_coefficients[5] = 0.0;
    }
  
    real64 getRadius(real64 angle) const
    {
      real64 radius = 0;
      integer count = 0;
      for(auto coeff:m_coefficients)
      {
        radius = radius + coeff*cos(count*angle);
        count++;
      }
      return radius;
    }

    void updateCoefficients(integer index, real64 value)
    {
      m_coefficients[index]=value;
    }

    array1d<real64> m_coefficients;
  
  };  


protected:

  /**
   * @brief This function provides the capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postProcessInput() override final;

private:

  /// center of the RTheta in (x,y,z) coordinates
  R1Tensor m_center;
  /// Dimensions of the RTheta's radius (as a function of theta)
  VariableRadius m_radius;
  /// tolerance to determine if a point sits on the RTheta or not
  real64 m_tolerance;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * centerString() { return "center"; }
    // static constexpr char const * radiusString() { return "radius"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
  };

  /// @endcond

};
} /* namespace geosx */

#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_RTheta_HPP_*/
