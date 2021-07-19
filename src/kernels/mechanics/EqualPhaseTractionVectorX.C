//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EqualPhaseTractionVectorX.h"
registerMooseObject("gibbsApp", EqualPhaseTractionVectorX);

template <>
InputParameters
validParams<EqualPhaseTractionVectorX>()
{
  InputParameters params = validParams<ObtainUnitNormalBase>();
  params.addClassDescription("Enforces equality of traction along x-direction");
  params.addRequiredCoupledVar(
      "ex_beta", "Strain in the beta phase"); 
  return params;
}

EqualPhaseTractionVectorX::EqualPhaseTractionVectorX(const InputParameters & parameters)
  : ObtainUnitNormalBase(parameters),
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
EqualPhaseTractionVectorX::computeQpResidual(){
  return _test[_i][_qp] * (_sx_beta[_qp] - _sx_alpha[_qp])*ObtainUnitNormalBase::nx();
}

Real
EqualPhaseTractionVectorX::computeQpJacobian(){
   // derivative with repect to _ex_alpha phase
  return -_test[_i][_qp] * _mat_const_alpha[_qp] * _phi[_j][_qp]*ObtainUnitNormalBase::nx();
}

Real
EqualPhaseTractionVectorX::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ex_beta_var){
    // derivative with respect to the coupled-variable ux_beta
    return _test[_i][_qp] * _mat_const_beta[_qp] * _phi[_j][_qp]*ObtainUnitNormalBase::nx();
  } 
  else
    return 0.0;     
}
