//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalDrivingForce.h"

registerMooseObject("gibbsApp",ChemicalDrivingForce);

template<>
InputParameters
validParams<ChemicalDrivingForce>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("xB_alpha", "Component B in alpha phase");
  params.addRequiredCoupledVar("xB_beta", "Component B in beta phase");
  return params;
}

ChemicalDrivingForce::ChemicalDrivingForce(const InputParameters & parameters)
  :AuxKernel(parameters),
  _xB_alpha(coupledValue("xB_alpha")),
  _xB_beta(coupledValue("xB_beta")),
  _dh(getMaterialProperty<Real>("dh")),
  _f_alpha(getMaterialProperty<Real>("f_alpha")),
  _f_beta(getMaterialProperty<Real>("f_beta")), 
  _B_diff_pot_alpha(getMaterialProperty<Real>("B_diff_pot_alpha")),
  _B_diff_pot_beta(getMaterialProperty<Real>("B_diff_pot_beta")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}
 
Real
ChemicalDrivingForce::computeValue()
{
 
  const Real omega_beta = (_f_beta[_qp] -_B_diff_pot_beta[_qp] * _xB_beta[_qp]);                                       
  const Real omega_alpha = (_f_alpha[_qp]- _B_diff_pot_alpha[_qp] *_xB_alpha[_qp]);
  
  return _nd_factor[_qp] * (omega_beta - omega_alpha)*_dh[_qp];
}
