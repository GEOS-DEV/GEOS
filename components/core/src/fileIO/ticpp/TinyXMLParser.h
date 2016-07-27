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
#ifndef __TINYXMLPARSER__
#define __TINYXMLPARSER__
#include "HierarchicalDataNode.hpp"
#include "tinyxml.h"
#include <sstream>
/*
 The MIT License
 
 Copyright (c) 2011 Grainflow Dynamics, Inc. (www.grainflow.com)
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */
namespace geosx
{
namespace TICPP
{
  ///This encapsulates the TinyXML implementation of the HDN
  /**
   This assumes that all data associated with the element is in the XML element "attributes"
   All sub-elements, therefore, are necessarily mapped to child HDN's
   */
  class TinyXMLParser
  {
  public:
    TinyXMLParser()
    {
    }

    ~TinyXMLParser()
    {
    }

  private:
    static TiXmlElement*
    Save(HierarchicalDataNode* hdn)
    {
      //set the name of the element
      TiXmlElement* me = new TiXmlElement(hdn->Heading());

      //create attributes associated with my node
      {
        std::string key, value;
        bool ok = true;
        for (ok = hdn->Next(key, value, true); ok; hdn->Next(key, value))
          me->SetAttribute(key.c_str(), value.c_str());
      }

      //create sub-elements associated with me
      {
        HierarchicalDataNode* chdn;
        for (chdn = hdn->Next(true); chdn; chdn = hdn->Next())
        {
          me->LinkEndChild(Save(chdn));
        }
      }
      return me;
    }

  public:
    static void Save(const char* pFilename, HierarchicalDataNode* hdn)

    {

      SLIC_ERROR_IF( !pFilename || hdn->Level() != 0, "Cannot save a file to an invalid file name or a hierarchical data node that is not root" );
      //create the document
      TiXmlDocument doc(pFilename);
      TiXmlDeclaration* decl = new TiXmlDeclaration("1.0", "", "");
      doc.LinkEndChild(decl);

      //create my node (and all my children)
      doc.LinkEndChild(Save(hdn));

      //save me!
      doc.SaveFile(pFilename);
    } // end Save

  private:
    static bool Load(TiXmlElement* me, HierarchicalDataNode* hdn)
    {
      //get my header
      hdn->SetHeading(me->Value());

      //get my attributes
      {
        TiXmlAttribute* pAttrib = 0;
        for (pAttrib = me->FirstAttribute(); pAttrib; pAttrib = pAttrib->Next())
        {
          hdn->AddAttributePair(pAttrib->Name(), pAttrib->Value());
        }
      }

      //get sub-elements
      {
        TiXmlElement* child = 0;
        for (child = me->FirstChildElement(); child; child
            = child->NextSiblingElement())
        {
          HierarchicalDataNode* chdn = hdn->NewChild(child->Value());
          Load(&(*child), chdn);
        }
      }
      return true;
    }
  public:

    static bool Load(const char* pFilename, HierarchicalDataNode* hdn)


    {
      //serialize world into the data structure
      TiXmlDocument doc(pFilename);
      if (!doc.LoadFile())
      {
        std::ostringstream oss;
        oss<<"TinyXml failed to load file "+ std::string(pFilename) +"\n"
           <<doc.ErrorDesc()
           <<"\nrow: "<<doc.ErrorRow()<<" column: "<<doc.ErrorCol();

        SLIC_ERROR( oss.str().c_str() );
        printf("TinyXml failed to LoadFile.\n%s\nrow: %d column: %d\n",doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
        return false;
      }
      TiXmlHandle hDoc(&doc);
      // block: name
      {
        TiXmlElement* pElem = hDoc.FirstChildElement().ToElement();
        // should always have a valid root but handle gracefully if it does
        SLIC_ERROR_IF( !pElem, "No valid root element during load of XML file");
        return Load(pElem, hdn);
      }
    }
    
    /// load from an xml file with multiple (or no) root elements
    static bool Load(const char* pFilename, std::vector<HierarchicalDataNode>& hdnv, int hdnLevel)


    {
        //serialize world into the data structure
        TiXmlDocument doc(pFilename);
        if (!doc.LoadFile())
        {
          std::ostringstream oss;
          oss<<"TinyXml failed to load file "+ std::string(pFilename) +"\n"
             <<doc.ErrorDesc()
             <<"\nrow: "<<doc.ErrorRow()<<" column: "<<doc.ErrorCol();

          SLIC_ERROR( oss.str().c_str() );
          printf("TinyXml failed to LoadFile.\n%s\nrow: %d column: %d\n",doc.ErrorDesc(), doc.ErrorRow(), doc.ErrorCol());
          return false;
        }
        TiXmlHandle hDoc(&doc);
        // block: name
        {
          TiXmlElement* pElem = hDoc.FirstChildElement().ToElement();
          bool rv = true;
          int count = 0;
          while (pElem && rv)
          {
            TICPP::HierarchicalDataNode hdn( hdnLevel );
            hdnv.push_back(hdn);
            rv = Load(pElem, &hdnv.back());
            ++count;
            pElem = hDoc.ChildElement(count).ToElement();
          }
          return rv;
        }
    }
    
  };
}
}
#endif
