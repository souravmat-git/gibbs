//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StrainTensor.h"
registerMooseObject("gibbsApp", StrainTensor);

template <>
InputParameters
validParams<StrainTensor>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("disp_x", "Displacement of a node in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("ex", "The compatible strain in the x-direction");
  return params;
}

StrainTensor::StrainTensor(const InputParameters & parameters)
  : Material(parameters),
   _grad_ux(coupledGradient("disp_x")),
   _ex_val(declareProperty<Real>(getParam<MaterialPropertyName>("ex")))
{
} 

void
StrainTensor::computeQpProperties()
{
    //x = 0, y = 1, z=2
   _ex_val[_qp] = _grad_ux[_qp](0);  
}
