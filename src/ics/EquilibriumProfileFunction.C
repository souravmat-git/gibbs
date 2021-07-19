//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "EquilibriumProfileFunction.h"
#include "FEProblem.h"
#include "MooseMesh.h"

registerMooseObject("gibbsApp", EquilibriumProfileFunction);

template <>
InputParameters validParams<EquilibriumProfileFunction>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addParam<Real>("barrier_height", 1.0, "Double-well barrier height");
  params.addParam<Real>("kappa", 1.0, "Interface energy coeffecient");
  return params;
}
EquilibriumProfileFunction::EquilibriumProfileFunction(const InputParameters & parameters)
  : InitialCondition(parameters),
   _W_height(getParam<Real>("barrier_height")),
   _kappa(getParam<Real>("kappa"))
{
}

Real
EquilibriumProfileFunction::value(const Point &p)
{
  const Real _factor = (std::sqrt(_W_height)/(std::sqrt(2.0*_kappa)));
  //p(0) is the x-coordinate
 return (0.5*(1.0 - std::tanh(p(0)*_factor)));
}
