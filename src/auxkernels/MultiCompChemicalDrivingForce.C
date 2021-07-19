//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiCompChemicalDrivingForce.h"

registerMooseObject("gibbsApp",MultiCompChemicalDrivingForce);

template<>
InputParameters
validParams<MultiCompChemicalDrivingForce>()
{
  InputParameters params = validParams<ChemicalDrivingForce>();
  params.addRequiredCoupledVar("xC_alpha","Component C in alpha phase");
  params.addRequiredCoupledVar("xC_beta", "Component C in beta phase");
  params.addRequiredParam<MaterialPropertyName>("C_diff_pot_alpha", "mu(C) - mu(A)");
  params.addRequiredParam<MaterialPropertyName>("C_diff_pot_beta",  "mu(C) -mu(A)");
  return params;
}

MultiCompChemicalDrivingForce::MultiCompChemicalDrivingForce(const InputParameters & parameters)
  :ChemicalDrivingForce(parameters),
  _xC_alpha(coupledValue("xC_alpha")),
  _xC_beta(coupledValue("xC_beta")),
  _C_diff_pot_alpha(getMaterialProperty<Real>("C_diff_pot_alpha")),
  _C_diff_pot_beta(getMaterialProperty<Real>("C_diff_pot_beta"))
{
}
 
Real
MultiCompChemicalDrivingForce::computeValue()
{
 
  const Real omega_beta = (_f_beta[_qp] -_B_diff_pot_beta[_qp] * _xB_beta[_qp]
                                        -_C_diff_pot_beta[_qp] * _xC_beta[_qp]);
                                       
  const Real omega_alpha= (_f_alpha[_qp]- _B_diff_pot_alpha[_qp] *_xB_alpha[_qp]
                                       - _C_diff_pot_alpha[_qp]* _xC_alpha[_qp]);                             
  
  return (omega_beta - omega_alpha)*_dh[_qp];
}
