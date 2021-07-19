//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RelativeDisplacementPerUnitLength.h"
registerMooseObject("gibbsApp", RelativeDisplacementPerUnitLength);

template<>
InputParameters
validParams<RelativeDisplacementPerUnitLength>()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the traction vector.");
  params.addRequiredCoupledVar("eta", "Phase-field variable");
  params.addRequiredCoupledVar("disp_x", "Displacement in the x-direction");
  params.addCoupledVar("disp_y", 0.0,"Displacement in the y-direction");
  params.addRequiredParam<MaterialPropertyName>("dlx_name", "Deformation vector in x -direction");
  params.addRequiredParam<MaterialPropertyName>("dly_name", "Deformation vectorin x -direction");
  params.addRequiredParam<std::string>("phase_name", "traction vector for a given phase");
  return params;
}

RelativeDisplacementPerUnitLength::RelativeDisplacementPerUnitLength(const InputParameters & parameters)
 : Material(parameters),
  _phase_name(getParam<std::string>("phase_name")),
  _grad_eta(coupledGradient("eta")),
  _grad_ux(coupledGradient("disp_x")),
  _grad_uy(coupledGradient("disp_y")),
  _dlx_name(getParam<MaterialPropertyName>("dlx_name")),
  _dly_name(getParam<MaterialPropertyName>("dly_name")),
  _dlx(declareProperty<Real>(_dlx_name)),
  _dly(declareProperty<Real>(_dly_name))       
{
}

void
RelativeDisplacementPerUnitLength::computeQpProperties(){  
  //Gradient components obtained from phase-field
   Real px = _grad_eta[_qp](0);
   Real py = _grad_eta[_qp](1);
   Real pz = _grad_eta[_qp](2);
  
  //Norm of the gradient vector
   Real mag = std::sqrt((px*px + py*py + pz*pz)); 
  
  //Calculate the unit normal
   Real nx = px/mag;
   Real ny = py/mag;
   //Real nz = pz/mag;
  
  _dlx[_qp] = _grad_ux[_qp](0) * nx + _grad_ux[_qp](1) * ny;
  _dly[_qp] = _grad_uy[_qp](0) * nx + _grad_uy[_qp](1) * ny;      
}
