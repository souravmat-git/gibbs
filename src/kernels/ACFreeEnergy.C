//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*In this code we explicitly code the free energy
//*The form of the free energy is f(eta) phi^(4)/4 - phi^(2)/2
//The non-linear variable that this kernel operates on is eta

#include "ACFreeEnergy.h"

registerMooseObject("gibbsApp", ACFreeEnergy);

template <>
InputParameters
validParams<ACFreeEnergy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Double well barrier with constant Mobility and Interfacial parameter");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  return params;
}

ACFreeEnergy::ACFreeEnergy(const InputParameters & parameters)
  : Kernel(parameters),
    _L(getMaterialProperty<Real>("mob_name"))
{
}

Real
ACFreeEnergy::computeQpResidual()
{
  return (_L[_qp] * _test[_i][_qp] * ( _u[_qp] * _u[_qp] * _u[_qp] - _u[_qp]));
}

Real
ACFreeEnergy::computeQpJacobian()
{
  return  (3.0* _L[_qp] * _test[_i][_qp] *_phi[_j][_qp] * ( _u[_qp] * _u[_qp] - 1));
}
