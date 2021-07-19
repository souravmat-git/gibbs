//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseStrainConstraint.h"
registerMooseObject("gibbsApp", TwoPhaseStrainConstraint);

template <>
InputParameters
validParams<TwoPhaseStrainConstraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (1-h(eta))*ex_alpha + h(eta)*ex_beta - ex = 0."
                             "non-linear variable of this kernel is s_x");
  params.addRequiredCoupledVar("eta", "to distinguish phases");
  params.addRequiredCoupledVar("ux", "Displacement in the x-direction");
  return params;
}

TwoPhaseStrainConstraint::TwoPhaseStrainConstraint(const InputParameters & parameters)
  :  Kernel(parameters),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _grad_ux(coupledGradient("ux")),
    _ux_var(coupled("ux")),
    _ex_alpha(getMaterialProperty<Real>("ex_alpha")),
    _ex_beta(getMaterialProperty<Real>("ex_beta")),
    _compliance_alpha(getMaterialProperty<Real>("compliance_alpha")),
    _compliance_beta(getMaterialProperty<Real>("compliance_beta")),
    _h(getMaterialProperty<Real>("h")),
    _dh(getMaterialProperty<Real>("dh"))
{
}

Real
TwoPhaseStrainConstraint::computeQpResidual(){
  // Variable on which the kernel operates: s_x
  // Note: ex_beta(sx), ex_alpha(sx)
  return _test[_i][_qp]*(_h[_qp]*_ex_beta[_qp] + (1.0-_h[_qp])* _ex_alpha[_qp]
                         -_grad_ux[_qp](0));
}

Real
TwoPhaseStrainConstraint::computeQpJacobian(){
  //Derivative with respect to s_x
  return (_test[_i][_qp] * (_h[_qp] * _compliance_beta[_qp] 
                    + (1.0-_h[_qp]) * _compliance_alpha[_qp])* _phi[_j][_qp]);
}

Real
TwoPhaseStrainConstraint::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _eta_var){
    return _test[_i][_qp]*(_ex_beta[_qp] - _ex_alpha[_qp])* _dh[_qp] * _phi[_j][_qp];
  }
  else if (jvar == _ux_var){
    return  -(_test[_i][_qp]* _grad_phi[_j][_qp](0));
  }
  else
    return 0.0;
}
