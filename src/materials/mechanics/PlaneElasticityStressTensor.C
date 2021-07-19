//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "PlaneElasticityStressTensor.h"
registerMooseObject("gibbsApp", PlaneElasticityStressTensor);

template <>
InputParameters
validParams<PlaneElasticityStressTensor>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("disp_x", "Displacement of a node in the x-direction");
  params.addRequiredCoupledVar("disp_y", "Displacement of a node in the y-direction");
  params.addRequiredParam<std::string>("phase_name", "Phase name of the material");
  return params;
}

PlaneElasticityStressTensor::PlaneElasticityStressTensor(const InputParameters & parameters)
 : Material(parameters),
  _grad_ux(coupledGradient("disp_x")),
  _grad_uy(coupledGradient("disp_y")),
  _phase_name(getParam<std::string>("phase_name")),
  _C11(getMaterialProperty<Real>("C11_"+ _phase_name)),
  _C12(getMaterialProperty<Real>("C12_"+ _phase_name)),
  _C22(getMaterialProperty<Real>("C22_"+ _phase_name)),
  _C66(getMaterialProperty<Real>("C66_"+ _phase_name)),
  //Declare the values of the material properties to compute
  _sxx_val(declareProperty<Real>("sx_" + _phase_name)),
  _syy_val(declareProperty<Real>("sy_" + _phase_name)),
  _sxy_val(declareProperty<Real>("sxy_" + _phase_name))
{
} 

void
PlaneElasticityStressTensor::computeQpProperties()
{
  _sxx_val[_qp] = _C11[_qp] * _grad_ux[_qp](0) + _C12[_qp] * _grad_uy[_qp](1);
  _syy_val[_qp] = _C12[_qp] * _grad_ux[_qp](0) + _C22[_qp] * _grad_uy[_qp](1);
  _sxy_val[_qp] = _C66[_qp] * (_grad_ux[_qp](1) + _grad_uy[_qp](0));
}
