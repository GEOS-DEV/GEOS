//
//  DerivedClassFinal.hpp
//  testRTTypes
//
//  Created by Omar Duran on 12/16/24.
//

#ifndef DerivedClassFinal_hpp
#define DerivedClassFinal_hpp

#include <stdio.h>
#include "BaseClass.hpp"

class Derived final : public BaseClass
{
public:
  explicit Derived();

  virtual ~Derived() noexcept override
  {};
};

#endif /* DerivedClassFinal_hpp */
