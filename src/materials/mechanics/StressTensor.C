//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "StressTensor.h"
registerMooseObject("gibbsApp", StressTensor);

//template <>
InputParameters
StressTensor::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("disp_x", "Displacement of a node in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("youngs_modulus", "Youngs modulus");
  params.addRequiredParam<MaterialPropertyName>("sx", "Stress in x-direction");
  return params;
}

StressTensor::StressTensor(const InputParameters & parameters)
  : Material(parameters),
   _grad_ux(coupledGradient("disp_x")),
   _E(getMaterialProperty<Real>(getParam<MaterialPropertyName>("youngs_modulus"))),
   _sx_val(declareProperty<Real>(getParam<MaterialPropertyName>("sx")))
{
}

void
StressTensor::computeQpProperties()
{
    //x = 0, y = 1, z=2
   _sx_val[_qp] = _E[_qp]*_grad_ux[_qp](0);
}
