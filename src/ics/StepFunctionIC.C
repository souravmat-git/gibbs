//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StepFunctionIC.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include <cmath>
using namespace std;

registerMooseObject("gibbsApp", StepFunctionIC);

template <>
InputParameters
validParams<StepFunctionIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addClassDescription(
      "Creates a step function for composition");
  params.addRequiredParam<Real>("phase_comp_alpha", "The value on left (xmin) boundary.");
  params.addRequiredParam<Real>("phase_comp_beta", "The value on right (xmax) boundary.");
  params.addRequiredParam<Real>("phase_comp_gamma", "The value on right (xmax) boundary.");
  params.addRequiredParam<Real>("phase_length", "Length of the intermetallic phase");
  return params;
}

StepFunctionIC::StepFunctionIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _xlength(_fe_problem.mesh().dimensionWidth(0)),
    _ylength(_fe_problem.mesh().dimensionWidth(1)),
    _phase_comp_alpha(getParam<Real>("phase_comp_alpha")),
    _phase_comp_beta(getParam<Real>("phase_comp_beta")),
    _phase_comp_gamma(getParam<Real>("phase_comp_gamma")),
    _phase_len(getParam<Real>("phase_length"))
{
}

Real
StepFunctionIC::value(const Point & p)
{
  // p(0) to the x-coordinate of the problem
  
  Real _x1 = (_xlength -_phase_len)*0.5;
  Real _x2 = (_xlength + _phase_len)*0.5;
  
  if (p(0) < _x1 && p(1) < _ylength)
  {
    return _phase_comp_alpha;
  } 
  else if ( p(0)>= _x1 &&  p(0) < _x2 && p(1) < _ylength)
  {
    return _phase_comp_gamma;
  }
  else 
    return _phase_comp_beta;
  
}


