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

#include "DrivingForceKHS.h"

registerMooseObject("gibbsApp", DrivingForceKHS);

template<>
InputParameters
validParams<DrivingForceKHS>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("disp_x", "Displacement in the x-direction");  
  return params;
}

DrivingForceKHS::DrivingForceKHS(const InputParameters & parameters)
  :AuxKernel(parameters),
  _grad_ux(coupledGradient("disp_x")),
  _dh(getMaterialProperty<Real>("dh")),
  _elastic_energy_alpha(getMaterialProperty<Real>("elastic_energy_alpha")),
  _elastic_energy_beta(getMaterialProperty<Real>("elastic_energy_beta"))
{
}
 
Real
DrivingForceKHS::computeValue()
{
  
  return _dh[_qp]*(_elastic_energy_beta[_qp] - _elastic_energy_alpha[_qp]);
}
