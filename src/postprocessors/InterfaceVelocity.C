/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "InterfaceVelocity.h"

registerMooseObject("gibbsApp", InterfaceVelocity);

template <>
InputParameters
validParams<InterfaceVelocity>()
{
  InputParameters params = validParams<ChangeOverTimePostprocessor>();
  params.addParam<Real>("initial_velocity", 0.0, "Initial velocity");
  return params;
}

InterfaceVelocity::InterfaceVelocity(const InputParameters & parameters)
  : ChangeOverTimePostprocessor(parameters),
  _v0(getParam<Real>("initial_velocity"))
{
}

Real
InterfaceVelocity::getValue()
{
  _change_in_position = ChangeOverTimePostprocessor::getValue(); 
  return (_v0 + (_change_in_position/_dt));
}
