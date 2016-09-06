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
/**
 * @file ContactManagerT.h
 * @author Scott Johnson
 * @date Jun 15, 2011
 */
#ifndef CONTACTMANAGERT_H_
#define CONTACTMANAGERT_H_

#include "ContactManagerBaseT.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "dataRepository/SynchronizedGroup.hpp"

/**
 * @author Scott Johnson
 * @brief Class to manage collections of face-face contacts
 */
class ContactManagerT : public ContactManagerBaseT {
public:
	/**
	 * @brief Contact manager constructor
	 * @author Scott Johnson
	 */
  ContactManagerT();//ExternalFaceManagerT* fm);

	/**
	 * @brief Contact manager destructor
	 * @author Scott Johnson
	 */
  virtual ~ContactManagerT();

  localIndex NumberOfPolygons() const;

public:
  void SetDefaultPolygonDimensions(const R1Tensor& xmin, const R1Tensor& xmax);

  inline realT DefaultPolygonDimension() const { return dd;}

  inline R1Tensor DefaultPolygonCenter() const { return ctr;}

  ///maps the discrete element to the indices in the _external_ face manager
  OrderedVariableOneToManyRelation& m_contactToIntersectionPolygonPointsMap;

  ///stateless cache of the points comprising the contact patches
  Array1dT<R1Tensor> m_intersectionPolygonPoints;

private:
  // FOR DEFAULT VISUALIZATION
  R1Tensor ctr;
  realT dd;
};

#endif /* CONTACTMANAGERT_H_ */
