//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef EQUILIBRIUMPROFILEFUNCTION_H
#define EQUILIBRIUMPROFILEFUNCTION_H

//MOOSE includes
#include "InitialCondition.h"
#include "FEProblem.h"
#include "MooseMesh.h"

class EquilibriumProfileFunction;

template<>
InputParameters validParams<EquilibriumProfileFunction>();

class EquilibriumProfileFunction : public InitialCondition
{
  public:
    EquilibriumProfileFunction(const InputParameters & parameters);
    
    //This member function is re-dedfined  from the base class
    virtual Real value(const Point & p) override;

  protected:
    const Real _W_height;  //Barrier height;
    const Real _kappa;    // interface energy coefficent
};
#endif //EQUILIBRIUMPROFILEFUNCTION_H
