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

/// Very simple sparse matrix data structure

#include <iostream>
#include <fstream>
#include "RCVSparse.h"

// sort functions
inline bool rowColumnOrder(const rcv& A, const rcv& B){  return (A.r < B.r)|| ( (A.r==B.r) && (A.c < B.c) );}
inline bool columnRowOrder(const rcv& A, const rcv& B){  return (A.c < B.c)|| ( (A.c==B.c) && (A.r < B.r) );}

/// Sum multiple entries with same row and column in K
/// and reduce the vector size
///
/// K = the sparse array stored in row, column value format
void ConsolidateSparseArray(array<rcv>& K){
  	sort (K.begin(), K.end(), rowColumnOrder); 
  	localIndex i = 0;
  	localIndex ii = 1; 
  	
  	while(ii < K.size()){
  	  	if(  K[i].c == K[ii].c && K[i].r == K[ii].r ){
  	  	  K[i].v +=K[ii].v;	
  	  	} else {
  	  	  ++i;
  	  	  if(i != ii){
  	  	    K[i] = K[ii];	
  	  	  } 	
  	  	}
  	  	++ii;
  	}
  	K.resize(i+1); 
}

void WriteSparseArrayToFile(std::string filename, const array<rcv>& K){
    std::ofstream outputfile(filename.c_str());
    for(unsigned int i =0; i < K.size(); ++i){
     	  outputfile <<  K[i].r << " " <<  K[i].c << " "  <<  K[i].v << "\n";  
    }
    outputfile.close();
}
