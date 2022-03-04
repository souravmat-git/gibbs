//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableUnitNormal_MP.h"
registerMooseObject("gibbsApp", VariableUnitNormal_MP);

//template<>
InputParameters
VariableUnitNormal_MP::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Base classes for defining unit normals");
  params.addRequiredCoupledVar("phase_alpha", "Phase-field variable for alpha");
  params.addRequiredCoupledVar("phase_beta", "Phase-field variable for beta");
  params.addRequiredCoupledVar("phase_gamma", "Phase-field variable for gamma");
  return params;
}

VariableUnitNormal_MP::VariableUnitNormal_MP(const InputParameters & parameters)
 : Material(parameters),
  _grad_phialpha(coupledGradient("phase_alpha")),
  _grad_phibeta(coupledGradient("phase_beta")),
  _grad_phigamma(coupledGradient("phase_gamma")),
  _n1(declareProperty<RealVectorValue>("n_ab")),
  _n2(declareProperty<RealVectorValue>("n_bc"))
{
}

void
VariableUnitNormal_MP::computeQpProperties(){

  const Real _un1_norm_sq = _grad_phibeta[_qp].norm_sq();
  const Real _un2_norm_sq = _grad_phigamma[_qp].norm_sq();

  if (_un1_norm_sq < std::pow(libMesh::TOLERANCE,3.0))
      _n1[_qp].zero();
  else
    _n1[_qp] = - _grad_phibeta[_qp]/std::sqrt(_un1_norm_sq);

  if (_un2_norm_sq < std::pow(libMesh::TOLERANCE,3.0))
    _n2[_qp].zero();
  else
    _n2[_qp] = - _grad_phigamma[_qp]/std::sqrt(_un2_norm_sq);

}
