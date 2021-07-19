//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "DrivingForceLarcheCahn.h"

registerMooseObject("gibbsApp", DrivingForceLarcheCahn);

template<>
InputParameters
validParams<DrivingForceLarcheCahn>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("stress_x", "Elastic stress in the x-direction");  
  return params;
}

DrivingForceLarcheCahn::DrivingForceLarcheCahn(const InputParameters & parameters)
  :AuxKernel(parameters),
  _stress_x(coupledValue("stress_x")),
  _dh(getMaterialProperty<Real>("dh")),
  _comp_energy_alpha(getMaterialProperty<Real>("comp_energy_alpha")),
  _comp_energy_beta(getMaterialProperty<Real>("comp_energy_beta"))
{
}
 
Real
DrivingForceLarcheCahn::computeValue()
{
  
  return _dh[_qp]*(_comp_energy_beta[_qp] - _comp_energy_alpha[_qp]);
}
