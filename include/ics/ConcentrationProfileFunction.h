//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CONCENTRATIONPROFILEFUNCTION_H
#define CONCENTRATIONPROFILEFUNCTION_H

//MOOSE includes
#include "EquilibriumProfileFunction.h"


class ConcentrationProfileFunction;

template<>
InputParameters validParams<ConcentrationProfileFunction>();

class ConcentrationProfileFunction : public EquilibriumProfileFunction
{
  public:
  
    ConcentrationProfileFunction(const InputParameters & parameters);
    virtual Real value(const Point &p) override;

  protected:

    Real h(const Point &p);  // Returns the interpolation function
    
    const Real _xB_alpha;    //Alpha phase composition IC;
    const Real _xB_beta;     // Beta phase compoistion IC;

};
#endif //CONCENTRATIONPROFILEFUNCTION_H
