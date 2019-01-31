/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * @file Segment.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SEGMENT_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SEGMENT_HPP

#include "dataRepository/ManagedGroup.hpp"

namespace geosx
{

class Segment : public dataRepository::ManagedGroup
{
public:

  explicit Segment( string const & name, dataRepository::ManagedGroup * const parent );
  ~Segment() override;

  Segment() = delete;
  Segment( Segment const &) = delete;
  Segment( Segment && ) = delete;

  string const & getSegmentName() const    { return m_segmentName; }
  void setSegmentName(string const & name) { m_segmentName = name; }

  struct viewKeyStruct
  {

    static constexpr auto segmentNameString      = "segmentName";

    using ViewKey = dataRepository::ViewKey;
    
    ViewKey segmentName      = { segmentNameString };
    
  } viewKeysSegment;

private:

  // name 
  string   m_segmentName;

  /*
  
    Below: the static properties attached to segments in AD-GPRS

    // input data 
    real64 length;	     // length of the segment
    real64 depth;	     // depth of the segment (defined at the center)
    SegmentType type;	     // type (TUBING or ANNULUS) of the segment, needed for DF model, default: TUBING
    real64 diameter;	     // diameter of the segment (or tubing if annulus exist)
    real64 roughness;	     // roughness of the segment
    real64 diameter_annulus; // diameter of the annulus
    real64 Tformation;	     // temperature of surroundings
    real64 Uto;		     // heat transfer coefficient
    
    // computed properties
    real64 incl;                // computed inclination of the segment
    real64 area;                // computed area of the segment
    real64 frictionCoefficient; // computed friction coefficient of the segment
    real64 volume;	        // computed volume of the segment
  
  */

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
