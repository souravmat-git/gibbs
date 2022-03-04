//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryConstantKineticMaterial.h"
registerMooseObject("gibbsApp", BinaryConstantKineticMaterial);

//template <>
InputParameters
BinaryConstantKineticMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("dL_BB_muB","L_BB w.r.t B");
  params.addRequiredParam<Real>("L_BB_val", "Input value for L_BB");
  params.addClassDescription("Constant mobility material for A-B alloy");
  return params;
}

BinaryConstantKineticMaterial::BinaryConstantKineticMaterial(const InputParameters & parameters)
  :Material(parameters),
  _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
  _dL_BB_muB_name(getParam<MaterialPropertyName>("dL_BB_muB")),
  _L_BB(declareProperty<Real>(_L_BB_name)),
  _dL_BB_muB(declareProperty<Real>(_dL_BB_muB_name)),
  _L_BB_val(getParam<Real>("L_BB_val"))
{
}

void
BinaryConstantKineticMaterial::computeQpProperties()
{
  //The Onsagermobility value, here asssumed to be a constant
  //and non-dimensional
  _L_BB[_qp] = _L_BB_val;

  //Its first derivative is zero
  _dL_BB_muB[_qp] = 0;
}
