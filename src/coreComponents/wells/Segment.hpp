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

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP

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
    static constexpr auto nextSegmentIndexString = "nextSegmentIndex";
    static constexpr auto prevSegmentIndexString = "prevSegmentIndex";

    using ViewKey = dataRepository::ViewKey;
    
    ViewKey segmentName      = { segmentNameString };
    ViewKey nextSegmentIndex = { nextSegmentIndexString };
    ViewKey prevSegmentIndex = { prevSegmentIndexString };
    
  } viewKeysSegment;

private:

  string   m_segmentName;
  localIndex nextSegmentIndex;
  localIndex prevSegmentIndex;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_PERFORATION_HPP
