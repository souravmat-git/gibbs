//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DeformationVector.h"
registerMooseObject("gibbsApp", DeformationVector);

template<>
InputParameters
validParams<DeformationVector>()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the traction vector.");
  params.addRequiredCoupledVar("eta", "Phase-field variable");
  params.addRequiredParam<MaterialPropertyName>("ax_name", "Deformation vector in x -direction");
  params.addRequiredParam<MaterialPropertyName>("ay_name", "Deformation vectorin x -direction");
  params.addRequiredParam<std::string>("phase_name", "traction vector for a given phase");
  return params;
}

DeformationVector::DeformationVector(const InputParameters & parameters)
 : Material(parameters),
  _phase_name(getParam<std::string>("phase_name")),
  _grad_eta(coupledGradient("eta")),
  _exx(getMaterialProperty<Real>("exx_" + _phase_name)),
  _eyy(getMaterialProperty<Real>("eyy_" + _phase_name)),
  _exy(getMaterialProperty<Real>("exy_" + _phase_name)),
  _ax_name(getParam<MaterialPropertyName>("ax_name")),
  _ay_name(getParam<MaterialPropertyName>("ay_name")),
  _ax(declareProperty<Real>(_ax_name)),
  _ay(declareProperty<Real>(_ay_name))       
{
}

void
DeformationVector::computeQpProperties()
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
  
  _ax[_qp] = _exx[_qp] * nx + _exy[_qp] * ny;
  _ay[_qp] = _exy[_qp] * nx + _eyy[_qp] * ny;      
}
