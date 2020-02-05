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
 * @file Tetrahedron.h
 */

#include "finiteElement/ElementLibrary/FiniteElementBase.h"

#ifndef TETRAHEDRON_H_
#define TETRAHEDRON_H_

namespace geosx
{
class TetrahedralElement : public FiniteElementBase
{
public:
  TetrahedralElement( BasisBase const & basis,
                      QuadratureBase const & quadrature,
                      const int num_zero_energy_modes );

  ~TetrahedralElement() override;

  static string CatalogName() { return "C3D4"; }

  void reinit( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X, arraySlice1d< localIndex const, -1 > const & mapped_support_points ) override
  { return reinitPrivate( X, mapped_support_points ); }

  void reinit( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X, arraySlice1d< localIndex const, 0 > const & mapped_support_points ) override
  { return reinitPrivate( X, mapped_support_points ); }

private:
  template< int USD >
  void reinitPrivate( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X, arraySlice1d< localIndex const, USD > const & mapped_support_points )
  {
    const unsigned int q = 0;

    localIndex const n0 = mapped_support_points[ 0 ];
    localIndex const n1 = mapped_support_points[ 1 ];
    localIndex const n2 = mapped_support_points[ 2 ];
    localIndex const n3 = mapped_support_points[ 3 ];

    real64 V, a[4], b[4], c[4];
    constexpr real64 sixth = 1.0 / 6.0;

    //a[1]=(y4-y2)*(z3-z2)-(y3-y2)*(z4-z2)
    a[0] = (X( n3, 1 ) - X( n1, 1 )) * (X( n2, 2 ) - X( n1, 2 )) - (X( n2, 1 ) - X( n1, 1 )) * (X( n3, 2 ) - X( n1, 2 ));
    //      b[1]=x32 z42 - x42 z32
    b[0] = (X( n2, 0 ) - X( n1, 0 )) * (X( n3, 2 ) - X( n1, 2 )) - (X( n3, 0 ) - X( n1, 0 )) * (X( n2, 2 ) - X( n1, 2 ));
    //      c[1]=x42 y32 - x32 y42
    c[0] = (X( n3, 0 ) - X( n1, 0 )) * (X( n2, 1 ) - X( n1, 1 )) - (X( n2, 0 ) - X( n1, 0 )) * (X( n3, 1 ) - X( n1, 1 ));
    //      a[2]=y31 z43 - y34 z13
    a[1] = (X( n2, 1 ) - X( n0, 1 )) * (X( n3, 2 ) - X( n2, 2 )) - (X( n2, 1 ) - X( n3, 1 )) * (X( n0, 2 ) - X( n2, 2 ));
    //      b[2]=x43 z31 - x13 z34
    b[1] = (X( n3, 0 ) - X( n2, 0 )) * (X( n2, 2 ) - X( n0, 2 )) - (X( n0, 0 ) - X( n2, 0 )) * (X( n2, 2 ) - X( n3, 2 ));
    //      c[2]=x31 y43 - x34 y13
    c[1] = (X( n2, 0 ) - X( n0, 0 )) * (X( n3, 1 ) - X( n2, 1 )) - (X( n2, 0 ) - X( n3, 0 )) * (X( n0, 1 ) - X( n2, 1 ));
    //      a[3]=y24 z14 - y14 z24
    a[2] = (X( n1, 1 ) - X( n3, 1 )) * (X( n0, 2 ) - X( n3, 2 )) - (X( n0, 1 ) - X( n3, 1 )) * (X( n1, 2 ) - X( n3, 2 ));
    //      b[3]=x14 z24 - x24 z14
    b[2] = (X( n0, 0 ) - X( n3, 0 )) * (X( n1, 2 ) - X( n3, 2 )) - (X( n1, 0 ) - X( n3, 0 )) * (X( n0, 2 ) - X( n3, 2 ));
    //      c[3]=x24 y14 - x14 y24
    c[2] = (X( n1, 0 ) - X( n3, 0 )) * (X( n0, 1 ) - X( n3, 1 )) - (X( n0, 0 ) - X( n3, 0 )) * (X( n1, 1 ) - X( n3, 1 ));
    //      a[4]=y13 z21 - y12 z31
    a[3] = (X( n0, 1 ) - X( n2, 1 )) * (X( n1, 2 ) - X( n0, 2 )) - (X( n0, 1 ) - X( n1, 1 )) * (X( n2, 2 ) - X( n0, 2 ));
    //      b[4]=x21 z13 - x31 z12
    b[3] = (X( n1, 0 ) - X( n0, 0 )) * (X( n0, 2 ) - X( n2, 2 )) - (X( n2, 0 ) - X( n0, 0 )) * (X( n0, 2 ) - X( n1, 2 ));
    //      c[4]=x13 y21 - x12 y31
    c[3] = (X( n0, 0 ) - X( n2, 0 )) * (X( n1, 1 ) - X( n0, 1 )) - (X( n0, 0 ) - X( n1, 0 )) * (X( n2, 1 ) - X( n0, 1 ));
    //
    //6V=x21 (y23 z34 - y34 z23) + x32 (y34 z12 - y12 z34) + x43 (y12 z23 - y23
    // z12),
    V = (X( n1, 0 ) - X( n0, 0 )) * ((X( n1, 1 ) - X( n2, 1 )) * (X( n2, 2 ) - X( n3, 2 )) - (X( n2, 1 ) - X( n3, 1 )) * (X( n1, 2 ) - X( n2, 2 ))) + (X( n2, 0 ) - X( n1, 0 )) *
        ((X( n2, 1 ) - X( n3, 1 )) * (X( n0, 2 ) - X( n1, 2 )) - (X( n0, 1 ) - X( n1, 1 )) * (X( n2, 2 ) - X( n3, 2 ))) + (X( n3, 0 ) - X( n2, 0 )) *
        ((X( n0, 1 ) - X( n1, 1 )) * (X( n1, 2 ) - X( n2, 2 )) - (X( n1, 1 ) - X( n2, 1 )) * (X( n0, 2 ) - X( n1, 2 )));
    V *= sixth;

    data[q].jacobian_determinant = V;
    for( int iNd=0 ; iNd<4 ; ++iNd )
    {
      data[ q ].mapped_gradients( iNd, 0 ) = a[ iNd ] * sixth / data[ q ].jacobian_determinant;
      data[ q ].mapped_gradients( iNd, 1 ) = b[ iNd ] * sixth / data[ q ].jacobian_determinant;
      data[ q ].mapped_gradients( iNd, 2 ) = c[ iNd ] * sixth / data[ q ].jacobian_determinant;
    }
  }

};
}

#endif /* TETRAHEDRON_H_ */
