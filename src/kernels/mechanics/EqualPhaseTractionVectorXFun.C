//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EqualPhaseTractionVectorXFun.h"
#include "Function.h"
registerMooseObject("gibbsApp", EqualPhaseTractionVectorXFun);

template <>
InputParameters
validParams<EqualPhaseTractionVectorXFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Enforces equality of traction along x-direction");
  params.addRequiredCoupledVar(
      "ex_beta", "Strain in the beta phase"); // xB_{beta} is the coupled variable
  //params.addRequiredParam<FunctionName>("eta", "To distinguish phases");
  return params;
}

EqualPhaseTractionVectorXFun::EqualPhaseTractionVectorXFun(const InputParameters & parameters)
  : Kernel(parameters),
   //_eta(getParam<Function>("eta")),
   _ex_beta(coupledValue("ex_beta")),
   _ex_beta_var(coupled("ex_beta")),
   //Material properties required for Equality of traction vectors
   _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
   _sx_beta(getMaterialProperty<Real>("sx_beta")),
   _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
   _mat_const_beta(getMaterialProperty<Real>("mat_const_beta"))
{
}

Real
EqualPhaseTractionVectorXFun::computeQpResidual(){
  return _test[_i][_qp] * (_sx_beta[_qp] - _sx_alpha[_qp]);
}

Real
EqualPhaseTractionVectorXFun::computeQpJacobian(){
   // derivative with repect to _ux_alpha phase
  return -_test[_i][_qp] * _mat_const_alpha[_qp] * _phi[_j][_qp];
}

Real
EqualPhaseTractionVectorXFun::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ex_beta_var){
    // derivative with respect to the coupled-variable ex_beta
    return _test[_i][_qp] * _mat_const_beta[_qp] * _phi[_j][_qp];
  }   
  else
    return 0.0;     
}
