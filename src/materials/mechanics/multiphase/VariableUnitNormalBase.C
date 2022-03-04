//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableUnitNormalBase.h"
registerMooseObject("gibbsApp", VariableUnitNormalBase);

//template<>
InputParameters
VariableUnitNormalBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Base classes for defining unit normals");
  params.addRequiredCoupledVar("phase_alpha", "Phase-field variable for alpha");
  params.addRequiredCoupledVar("phase_beta", "Phase-field variable for beta");
  params.addRequiredCoupledVar("phase_gamma", "Phase-field variable for gamma");
  return params;
}

VariableUnitNormalBase::VariableUnitNormalBase(const InputParameters & parameters)
 : Material(parameters),
  _grad_phialpha(coupledGradient("phase_alpha")),
  _grad_phibeta(coupledGradient("phase_beta")),
  _grad_phigamma(coupledGradient("phase_gamma"))
{
}

Real
VariableUnitNormalBase::unalpha_norm_sq() const {
   return _grad_phialpha[_qp].norm_sq();
}

Real
VariableUnitNormalBase::unbeta_norm_sq() const {
  return _grad_phibeta[_qp].norm_sq();
}

Real
VariableUnitNormalBase::ungamma_norm_sq() const {
  return _grad_phigamma[_qp].norm_sq();
}

void
VariableUnitNormalBase::computeQpProperties(){
}
