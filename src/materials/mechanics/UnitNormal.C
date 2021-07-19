//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "UnitNormal.h"
registerMooseObject("gibbsApp", UnitNormal);

template<>
InputParameters
validParams<UnitNormal>()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the projection tensor from the PFV");
  params.addRequiredParam<RealVectorValue>("n", "Constant unit normal");
  params.addParam<MaterialPropertyName>("unit_normal_name","n", "Unit normal name");
  
  return params;
}

UnitNormal::UnitNormal(const InputParameters & parameters)
 : Material(parameters),
  _n(getParam<RealVectorValue>("n")),
  _unit_normal_name(getParam<MaterialPropertyName>("unit_normal_name")),
  _n_val(declareProperty<RealVectorValue>(_unit_normal_name))                
{
}

void
UnitNormal::computeQpProperties()
{ 
  _n_val[_qp] = _n;
}
