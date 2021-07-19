//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoDimEqualPhaseTractionVectorYFun.h"
registerMooseObject("gibbsApp", TwoDimEqualPhaseTractionVectorYFun);

template <>
InputParameters
validParams<TwoDimEqualPhaseTractionVectorYFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Enforces equality of traction along x-direction");
  params.addRequiredCoupledVar(
      "eyy_beta", "Strain component yy in the beta phase"); 
  return params;
}

TwoDimEqualPhaseTractionVectorYFun::TwoDimEqualPhaseTractionVectorYFun(const InputParameters & parameters)
  : Kernel(parameters),
   //_eyy_beta(coupledValue("eyy_beta")),
   _eyy_beta_var(coupled("eyy_beta")),
   //Material properties required for Equality of traction vectors
   _syy_alpha(getMaterialProperty<Real>("syy_alpha")),
   _syy_beta(getMaterialProperty<Real>("syy_beta")),
   _C22_alpha(getMaterialProperty<Real>("C22_alpha")),
   _C22_beta(getMaterialProperty<Real>("C22_beta"))
{
}

Real
TwoDimEqualPhaseTractionVectorYFun::computeQpResidual(){
  return _test[_i][_qp] * (_syy_beta[_qp] - _syy_alpha[_qp]);
}

Real
TwoDimEqualPhaseTractionVectorYFun::computeQpJacobian(){
   // derivative with repect to _eyy_alpha phase
  return -_test[_i][_qp] * _C22_alpha[_qp] * _phi[_j][_qp];
}

Real
TwoDimEqualPhaseTractionVectorYFun::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _eyy_beta_var){
    // derivative with respect to the coupled-variable ux_beta
    return _test[_i][_qp] * _C22_beta[_qp] * _phi[_j][_qp];
  } 
  else
    return 0.0;     
}
