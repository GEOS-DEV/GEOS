/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file NumericalMethodsManager.hpp
 */

#ifndef GEOSX_MANAGERS_NUMERICALMETHODSMANAGER_HPP_
#define GEOSX_MANAGERS_NUMERICALMETHODSMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "finiteElement/basis/BasisFunctionManager.hpp"
#include "finiteElement/quadrature/QuadratureRuleManager.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{}
}

/**
 * @class NumericalMethodsManager
 *
 * Contains the discretization methods and their components for application to the mesh.
 */
class NumericalMethodsManager : public dataRepository::Group
{
public:
  /// Deleted default constructor.
  NumericalMethodsManager() = delete;

  /**
   * @brief Constructor
   *
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  NumericalMethodsManager( string const & name, Group * const parent );

  virtual ~NumericalMethodsManager() override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// Contains the keys for the object names in the data repository.
  struct groupKeysStruct
  {
    /// Name for the basis function manager.
    static constexpr auto basisFunctions = "BasisFunctions";

    /// Name for the quadrature rule manager.
    static constexpr auto quadratureRules = "QuadratureRules";

    /// Name for the finite element discretization manager.
    static constexpr auto finiteElementDiscretizations = "FiniteElements";

    /// Name for the finite volume manager.
    static constexpr auto finiteVolumeManager = "FiniteVolume";

  };

  /**
   * @brief @return Returns reference to const BasisFunctionManager m_basisFunctions.
   */
  BasisFunctionManager const & getBasisFunctions() const { return m_basisFunctions; }

  /**
   * @brief @return Returns reference to BasisFunctionManager m_basisFunctions.
   */
  BasisFunctionManager & getBasisFunctions()       { return m_basisFunctions; }

  /**
   * @brief @return Returns reference to const QuadratureRuleManager m_quadratureRules.
   */
  QuadratureRuleManager const & getQuadratureRules() const { return m_quadratureRules; }

  /**
   * @brief @return Returns reference to QuadratureRuleManager m_quadratureRules.
   */
  QuadratureRuleManager & getQuadratureRules()       { return m_quadratureRules; }

  /**
   * @brief @return Returns reference to const FiniteElementDiscretizationManager m_finiteElementDiscretizationManager.
   */
  FiniteElementDiscretizationManager const & getFiniteElementDiscretizationManager() const { return m_finiteElementDiscretizationManager; }

  /**
   * @brief @return Returns reference to FiniteElementDiscretizationManager m_finiteElementDiscretizationManager.
   */
  FiniteElementDiscretizationManager & getFiniteElementDiscretizationManager()       { return m_finiteElementDiscretizationManager; }

  /**
   * @brief @return Returns reference to FiniteVolumeManager m_finiteVolumeManager.
   */
  FiniteVolumeManager & getFiniteVolumeManager()       { return m_finiteVolumeManager; }

  /**
   * @brief @return Returns reference to const FiniteVolumeManager m_finiteVolumeManager.
   */
  FiniteVolumeManager const & getFiniteVolumeManager() const { return m_finiteVolumeManager; }

private:
  /// Contains basis functions available for use in the various discretization packages.
  BasisFunctionManager m_basisFunctions;

  /// Contains quadrature rules available for use in the various discretization packages.
  QuadratureRuleManager m_quadratureRules;

  /// Contains the finite element discretizations
  FiniteElementDiscretizationManager m_finiteElementDiscretizationManager;

  /// Contains the finite volume discretizations.
  FiniteVolumeManager m_finiteVolumeManager;

};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_NUMERICALMETHODSMANAGER_HPP_ */
