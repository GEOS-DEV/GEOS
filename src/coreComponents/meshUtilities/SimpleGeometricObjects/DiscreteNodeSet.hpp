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

/*
 * @file DiscreteNodeSet.hpp
 *
 *  Created on: Feb 24, 2020
 *      Author: ron
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_DISCRETENODESET_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_DISCRETENODESET_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class DiscreteNodeSet : public SimpleGeometricObjectBase
{
public:
	DiscreteNodeSet( const std::string& name,
		       Group * const parent );

	virtual ~DiscreteNodeSet() override;

	static string CatalogName() { return "DiscreteNodeSet"; }

	bool IsCoordInObject( const R1Tensor& coord ) const override final;

	/*
	 * @brief Parse a table file.
	 * @tparam T The type for table or axis values.
	 * @param[in] target The place to store values.
	 * @param[in] filename The name of the file to read.
     * @param[in] delimiter The delimiter used for file entries.
     */
	  template< typename T >
	  void parse_file( array1d<T> & target, string const & filename, char delimiter );

protected:
	virtual void PostProcessInput() override final;

private:

	string m_nodeSetFile;

	struct viewKeyStruct
	{
		static constexpr auto nodeSetFileString = "nodeSetFile";
	};


};

} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_DISCRETENODESET_HPP_ */
