//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConstantMaterialModulus.h"
registerMooseObject("gibbsApp", ConstantMaterialModulus);

template <>
InputParameters
validParams<ConstantMaterialModulus>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<std::string>("base_name", "name of the phase");
  //params.addRequiredParam<MaterialPropertyName>("youngs_modulus", "Youngs modulus");
  //params.addRequiredParam<MaterialPropertyName>("poissons_ratio", "Poissons ratio");
  //params.addRequiredParam<MaterialPropertyName>("mat_const","material constant for a phase");
  params.addClassDescription("Constant material modulus for 1D problems");
  return params;
}

ConstantMaterialModulus::ConstantMaterialModulus(const InputParameters & parameters)
  : Material(parameters),
   _phase_name(getParam<std::string>("base_name")),
   _stiffness(getMaterialProperty<RankFourTensor>(_phase_name + "_elasticity_tensor")),
   _mat_const(declareProperty<Real>(_phase_name + "_mat_const"))
{
}

void
ConstantMaterialModulus::computeQpProperties(){

  //Fetch lames parameters;
  //Initialize the Lames constant //C_1122 = C12
 Real _lambda = _stiffness[_qp](0,0,1,1);
 //mu = C_2323 = C_44 (Voigt Notation)
 Real  _mu    = _stiffness[_qp](1,2,1,2);

 //Then declare mat constant as
 _mat_const[_qp] = (_lambda + 2 * _mu);

}
