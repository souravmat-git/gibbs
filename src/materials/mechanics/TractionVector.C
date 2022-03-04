//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TractionVector.h"
registerMooseObject("gibbsApp", TractionVector);

//template<>
InputParameters
TractionVector::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the traction vector.");
  params.addRequiredCoupledVar("eta", "Phase-field variable");
  params.addRequiredParam<MaterialPropertyName>("tx_name", "Traction vector in x -direction");
  params.addRequiredParam<MaterialPropertyName>("ty_name", "Traction vectorin x -direction");
  params.addRequiredParam<std::string>("phase_name", "traction vector for a given phase");
  return params;
}

TractionVector::TractionVector(const InputParameters & parameters)
 : Material(parameters),
  _phase_name(getParam<std::string>("phase_name")),
  _grad_eta(coupledGradient("eta")),
  _sxx(getMaterialProperty<Real>("sxx_" + _phase_name)),
  _syy(getMaterialProperty<Real>("syy_" + _phase_name)),
  _sxy(getMaterialProperty<Real>("sxy_" + _phase_name)),
  _tx_name(getParam<MaterialPropertyName>("tx_name")),
  _ty_name(getParam<MaterialPropertyName>("ty_name")),
  _tx(declareProperty<Real>(_tx_name)),
  _ty(declareProperty<Real>(_ty_name))
{
}

void
TractionVector::computeQpProperties()
{
  //Gradient components obtained from phase-field
  const Real px = _grad_eta[_qp](0);
  const Real py = _grad_eta[_qp](1);
  const Real pz = _grad_eta[_qp](2);

  //Norm of the gradient vector
  const Real mag = std::sqrt((px*px + py*py + pz*pz));

  //Calculate the unit normal
  const Real nx = px/mag;
  const Real ny = py/mag;
  //const Real nz = pz/mag;

  _tx[_qp] = _sxx[_qp] * nx + _sxy[_qp] * ny;
  _ty[_qp] = _sxy[_qp] * nx + _syy[_qp] * ny;
}
