//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoDimEqualPhaseTractionVectorX.h"
registerMooseObject("gibbsApp", TwoDimEqualPhaseTractionVectorX);

template <>
InputParameters
validParams<TwoDimEqualPhaseTractionVectorX>()
{
  InputParameters params = validParams<ObtainUnitNormalBase>();
  params.addClassDescription("Enforces equality of traction along x-direction");
  params.addRequiredCoupledVar(
      "exy_beta", "Strain component xy in the beta phase"); 
  return params;
}

TwoDimEqualPhaseTractionVectorX::TwoDimEqualPhaseTractionVectorX(const InputParameters & parameters)
  : ObtainUnitNormalBase(parameters),
   //_exy_beta(coupledValue("exy_beta")),
   _exy_beta_var(coupled("exy_beta")),
   //Material properties required for Equality of traction vectors
   _sxy_alpha(getMaterialProperty<Real>("sxy_alpha")),
   _sxy_beta(getMaterialProperty<Real>("sxy_beta")),
   _C66_alpha(getMaterialProperty<Real>("C66_alpha")),
   _C66_beta(getMaterialProperty<Real>("C66_beta"))
{
}

Real
TwoDimEqualPhaseTractionVectorX::computeQpResidual(){
  return _test[_i][_qp] * (_sxy_beta[_qp] - _sxy_alpha[_qp])*ObtainUnitNormalBase::ny();
}

Real
TwoDimEqualPhaseTractionVectorX::computeQpJacobian(){
   // derivative with repect to _exy_alpha phase
  return -_test[_i][_qp] * _C66_alpha[_qp] * _phi[_j][_qp]*ObtainUnitNormalBase::ny();
}

Real
TwoDimEqualPhaseTractionVectorX::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _exy_beta_var){
    // derivative with respect to the coupled-variable ux_beta
    return _test[_i][_qp] * _C66_beta[_qp] * _phi[_j][_qp]*ObtainUnitNormalBase::ny();
  } 
  else
    return 0.0;     
}
