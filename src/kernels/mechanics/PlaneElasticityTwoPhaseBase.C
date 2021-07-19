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

#include "PlaneElasticityTwoPhaseBase.h"
registerMooseObject("gibbsApp", PlaneElasticityTwoPhaseBase);

template <>
InputParameters
validParams<PlaneElasticityTwoPhaseBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  params.addRequiredParam<MaterialPropertyName>("h", "Interpolation function");
  params.addClassDescription("This kernel is the base kernel for plane elasticity two-phase problems");
  return params;
}

PlaneElasticityTwoPhaseBase::PlaneElasticityTwoPhaseBase(const InputParameters & parameters)
  :Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  //Interpolation function and its derivate
  _h_name(getParam<MaterialPropertyName>("h")),
  _h(getMaterialProperty<Real>(_h_name)),
  _dh(getMaterialProperty<Real>("dh")),
  
   //Traction force normal to the x-plane
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
   //Traction force normal to the x-plane
  _sxy_alpha(getMaterialProperty<Real>("sxy_alpha")),
  _sxy_beta(getMaterialProperty<Real>("sxy_beta")),
   //Traction force normal to the y-plane
  _sy_alpha(getMaterialProperty<Real>("sy_alpha")),
  _sy_beta(getMaterialProperty<Real>("sy_beta")),
  
  //C11
  _C11_alpha(getMaterialProperty<Real>("C11_alpha")),
  _C11_beta(getMaterialProperty<Real>("C11_beta")),
  //C12
  _C12_alpha(getMaterialProperty<Real>("C12_alpha")),
  _C12_beta(getMaterialProperty<Real>("C12_beta")),
  //_C22
  _C22_alpha(getMaterialProperty<Real>("C22_alpha")),
  _C22_beta(getMaterialProperty<Real>("C22_beta")),
  //_C66
  _C66_alpha(getMaterialProperty<Real>("C66_alpha")),
  _C66_beta(getMaterialProperty<Real>("C66_beta"))
{
}


Real
PlaneElasticityTwoPhaseBase::C11_interp() const{
  return (_h[_qp] * _C11_beta[_qp] + (1.0 - _h[_qp])* _C11_alpha[_qp]);
}

Real
PlaneElasticityTwoPhaseBase::C12_interp() const{
  return (_h[_qp] * _C12_beta[_qp] + (1.0 - _h[_qp])* _C12_alpha[_qp]);
}

Real
PlaneElasticityTwoPhaseBase::C22_interp() const{
  return (_h[_qp] * _C22_beta[_qp] + (1.0 - _h[_qp])* _C22_alpha[_qp]);
}

Real
PlaneElasticityTwoPhaseBase::C66_interp() const{
  return (_h[_qp] * _C66_beta[_qp] + (1.0 - _h[_qp])* _C66_alpha[_qp]);
}


Real
PlaneElasticityTwoPhaseBase::computeQpResidual(){
  return 0;
}

Real
PlaneElasticityTwoPhaseBase::computeQpJacobian(){
  //Derivative with respect to the non-linear variable
  return 0;
}

Real
PlaneElasticityTwoPhaseBase::computeQpOffDiagJacobian(unsigned int /*jvar*/){ 
    return 0;
}
