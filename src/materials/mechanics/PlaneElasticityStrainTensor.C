//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "PlaneElasticityStrainTensor.h"
registerMooseObject("gibbsApp", PlaneElasticityStrainTensor);

template <>
InputParameters
validParams<PlaneElasticityStrainTensor>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("disp_x", "Displacement of a node in the x-direction");
  params.addRequiredCoupledVar("disp_y", "Displacement of a node in the y-direction");
  params.addRequiredParam<std::string>("phase_name", "Phase name of the material");
  return params;
}

PlaneElasticityStrainTensor::PlaneElasticityStrainTensor(const InputParameters & parameters)
  : Material(parameters),
  _grad_ux(coupledGradient("disp_x")),
  _grad_uy(coupledGradient("disp_y")),
  _phase_name(getParam<std::string>("phase_name")),
  //Declare the values of the material properties to compute
  _exx_val(declareProperty<Real>("ex_" + _phase_name)),
  _eyy_val(declareProperty<Real>("ey_" + _phase_name)),
  _exy_val(declareProperty<Real>("exy_" + _phase_name))
{
} 

void
PlaneElasticityStrainTensor::computeQpProperties()
{
  _exx_val[_qp] = _grad_ux[_qp](0);
  _eyy_val[_qp] = _grad_uy[_qp](1);
  _exy_val[_qp] = 0.5*(_grad_ux[_qp](1) + _grad_uy[_qp](0));
}
