//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ObtainUnitNormalBase.h"
registerMooseObject("gibbsApp", ObtainUnitNormalBase);

template<>
InputParameters
validParams<ObtainUnitNormalBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Returns the unit normal.");
  params.addParam<Real>("tol", 1e-45, "Tolerance for the gradient");
  params.addParam<Real>("val", 1.0, "value of unit normal if gradient is below the tolerance");
  params.addRequiredCoupledVar("eta", "Phase-field variable");
  return params;
}

ObtainUnitNormalBase::ObtainUnitNormalBase(const InputParameters & parameters)
 : Kernel(parameters),
  _tol(getParam<Real>("tol")),
  _val(getParam<Real>("val")),
  _grad_eta(coupledGradient("eta"))
{
}


Real 
ObtainUnitNormalBase::nx() const{
  if (_grad_eta[_qp].norm() > _tol) 
   return -_grad_eta[_qp](0)/_grad_eta[_qp].norm();
  else 
   return _val;
}

Real 
ObtainUnitNormalBase::ny() const{
  if (_grad_eta[_qp].norm() > _tol) 
   return -_grad_eta[_qp](1)/_grad_eta[_qp].norm();
  else
   return 0;
}

Real 
ObtainUnitNormalBase::nz() const{
  if (_grad_eta[_qp].norm() > _tol) 
   return -_grad_eta[_qp](2)/_grad_eta[_qp].norm();
  else
  return _val;
}


Real
ObtainUnitNormalBase::computeQpResidual(){  
  return 0;
}

Real
ObtainUnitNormalBase::computeQpJacobian(){
  return 0;
}

Real
ObtainUnitNormalBase::computeQpOffDiagJacobian(unsigned int /*jvar*/){
  return 0;
}

