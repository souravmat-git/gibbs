//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Moles.h"
registerMooseObject("gibbsApp", Moles);

template <>
InputParameters
validParams<Moles>()
{
  InputParameters params = validParams<ElementIntegralVariablePostprocessor>();
  params.addParam<Real>("volume", 1.0 , "Volume of the system (3D)");
  return params;
}

Moles::Moles(const InputParameters & parameters)
  :ElementIntegralVariablePostprocessor(parameters),
  _volume(getParam<Real>("volume")) 
{

}

Real
Moles::computeQpIntegral()
{
  return (_u[_qp]/_volume);
}
