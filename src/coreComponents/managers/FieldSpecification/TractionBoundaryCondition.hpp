/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TractionBoundaryCondition.hpp
 */


#ifndef GEOSX_MANAGERS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_HPP
#define GEOSX_MANAGERS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_HPP

#include "FieldSpecificationBase.hpp"
#include "mesh/FaceManager.hpp"

namespace geosx
{

class TableFunction;

/**
 * @class TractionBoundaryCondition
 * Holds data and methods to apply a traction boundary condition
 */
class TractionBoundaryCondition : public FieldSpecificationBase
{
public:
  /// @copydoc FieldSpecificationBase(string const &, Group *)
  TractionBoundaryCondition( string const & name, Group * parent );

  /// deleted default constructor
  TractionBoundaryCondition() = delete;

  /// default destructor
  virtual ~TractionBoundaryCondition() = default;

  /// deleted copy constructor
  TractionBoundaryCondition( TractionBoundaryCondition const & other ) = delete;

  /// defaulted move constructor
  TractionBoundaryCondition( TractionBoundaryCondition && other ) = default;

  /// deleted copy assignment operator
  TractionBoundaryCondition & operator=( TractionBoundaryCondition const & other ) = delete;

  /// deleted move assignment operator
  TractionBoundaryCondition & operator=( TractionBoundaryCondition && other ) = delete;


  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Traction"; }


  /**
   * @brief Setup and Launche of the traction BC kernel.
   * @param time The time that should be used to evaluate any time tables.
   * @param blockLocalDofNumber Array of block local DOF numbers for the displacement.
   * @param dofRankOffset The rank offset for the DOF.
   * @param faceManager Reference to the face manager (Tractions are applied on faces)
   * @param targetSet The set of faces to apply the BC to.
   * @param localRhs The RHS of the system to add contributions to.
   */
  void launch( real64 const time,
               arrayView1d< globalIndex const > const blockLocalDofNumber,
               globalIndex const dofRankOffset,
               FaceManager const & faceManager,
               SortedArrayView< localIndex const > const & targetSet,
               arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for tractionType
    constexpr static char const * tractionTypeString() { return "tractionType"; }

    /// @return The key for inputStress.
    constexpr static char const * inputStressString() { return "inputStress"; }

//    /// @return The key for the function describing the components of stress.
//    constexpr static char const * stressFunctionString() { return "stressFunctions"; }

  };

protected:
  virtual void postProcessInput() override final;

  virtual void initializePreSubGroups() override final;

  /// The type of traction to be applied, i.e. how to generate the traction.
  int m_tractionType;

  /// single specified value for stress used to generate the traction if m_tractionType==2.
  VoigtTensor m_inputStress;

//  /// names of the functions used to specify stress for the generation of tractions.
//  array1d<string> m_stressFunctionNames;
//
//  bool m_useStressFunctions;
//
//  TableFunction const * m_stressFunctions[6];

};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FIELDSPECIFICATION_TRACTIONBOUNDARYCONDITION_HPP */
