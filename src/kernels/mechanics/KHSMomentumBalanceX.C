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

#include "KHSMomentumBalanceX.h"
registerMooseObject("gibbsApp", KHSMomentumBalanceX);

template <>
InputParameters
validParams<KHSMomentumBalanceX>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  params.addRequiredCoupledVar("eta", "Phase field variable");
  return params;
}

KHSMomentumBalanceX::KHSMomentumBalanceX(const InputParameters & parameters)
  :Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _mat_const_beta(getMaterialProperty<Real>("mat_const_beta")),
  _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh"))
{
}

Real
KHSMomentumBalanceX::computeQpResidual(){
    Real _interpolated_stress = (_h[_qp]* _sx_beta[_qp] + (1.0-_h[_qp]) * _sx_alpha[_qp]);        
    return (-_grad_test[_i][_qp](0) * _interpolated_stress);          
}

Real
KHSMomentumBalanceX::computeQpJacobian(){
  Real _interpolated_modulus = (_h[_qp]* _mat_const_beta[_qp] + (1.0-_h[_qp]) * _mat_const_alpha[_qp]);
  return (-_grad_test[_i][_qp](0) * _interpolated_modulus * _grad_phi[_j][_qp](0));
}

Real
KHSMomentumBalanceX::computeQpOffDiagJacobian(unsigned int jvar){
 if (jvar == _eta_var){
 
   Real  _diff_stress = (_sx_beta[_qp] - _sx_alpha[_qp]);
   return (-_grad_test[_i][_qp](0) * _diff_stress * _dh[_qp] * _phi[_j][_qp]);
 }  
 else
    return 0.0;
}
