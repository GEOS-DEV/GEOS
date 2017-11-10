//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MaterialBlock.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------
MaterialBlock::MaterialBlock( double rho, double vs, double vp, double xmin, 
                              double xmax, double ymin, double ymax, double zmin, double zmax,
			      double qs, double qp, double freq )
{
   m_rho = rho;
   m_vp  = vp;
   m_vs  = vs;
   m_xmin = xmin;
   m_xmax = xmax;
   m_ymin = ymin;
   m_ymax = ymax;
   m_zmin = zmin;
   m_zmax = zmax;
   m_tol = 1e-5;
   m_vpgrad  = 0;
   m_vsgrad  = 0;
   m_rhograd = 0;
   m_qs = qs;
   m_qp = qp;
   m_freq = freq;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_absoluteDepth( bool absDepth )
{
   m_absoluteDepth = absDepth;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_gradients( double rhograd, double vsgrad, double vpgrad )
{
   m_rhograd = rhograd;
   m_vsgrad  = vsgrad;
   m_vpgrad  = vpgrad;
}

//-----------------------------------------------------------------------
bool MaterialBlock::inside_block( double x, double y, double z )
{
   return m_xmin-m_tol <= x && x <= m_xmax+m_tol && m_ymin-m_tol <= y && 
    y <= m_ymax+m_tol &&  m_zmin-m_tol <= z && z <= m_zmax+m_tol;
}

//-----------------------------------------------------------------------
void MaterialBlock::set_material_properties( array<realT> & rho,
					     array<realT>& cs, array<realT> & cp,
					     array<R1Tensor>& coord, realT zsurf )
{
   for( localIndex ind = 0 ; ind < rho.size() ; ind++ )
      if(inside_block(coord[ind][0],coord[ind][1],coord[ind][2]))
      {
	 if( m_rho != -1 )
	    rho[ind] = m_rho + m_rhograd*(coord[ind][2]-zsurf);
	 if( m_vs != -1 )
	    cs[ind]  = m_vs + m_vsgrad*(coord[ind][2]-zsurf);
	 if( m_vp != -1 )
	    cp[ind]  = m_vp + m_vpgrad*(coord[ind][2]-zsurf);
      }
}



