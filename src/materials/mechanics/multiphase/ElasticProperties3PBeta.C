//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElasticProperties3PBeta.h"
registerMooseObject("gibbsApp", ElasticProperties3PBeta);

//template <>
InputParameters
ElasticProperties3PBeta::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<std::string>("base_name","The phase name");
  //params.addRequiredParam<MaterialPropertyName>("strain_jump", "Jump in strain");
  //params.addRequiredParam<MaterialPropertyName>("h", "Interpolation Function");
  return params;
}

ElasticProperties3PBeta::ElasticProperties3PBeta(const InputParameters & parameters)
  : Material(parameters),
   _phase_name(getParam<std::string>("base_name")),
   _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
   _a_alpha_beta(getMaterialProperty<Real>("a_alpha_beta")),
   _a_beta_gamma(getMaterialProperty<Real>("a_alpha_beta")),
    //interpolation function
   _h_alpha(getMaterialProperty<Real>("h_alpha")),
   _h_gamma(getMaterialProperty<Real>("h_gamma")),
   //material constant
   _mat_const(getMaterialProperty<Real>(_phase_name + "_mat_const")),
   //Calculate strain, stress, and strain energy
   _elastic_strain_val(declareProperty<Real>(_phase_name + "_elastic_strain")),
   _elastic_stress_val(declareProperty<Real>(_phase_name + "_elastic_stress")),
   _strain_energy_val(declareProperty<Real>(_phase_name + "_strain_energy"))
{
}

void
ElasticProperties3PBeta::computeQpProperties()
{
  //Calculat elastic strain
   _elastic_strain_val[_qp] = _total_strain[_qp](0,0)
                            - _h_alpha[_qp] * _a_alpha_beta[_qp]
                            + _h_gamma[_qp] * _a_beta_gamma[_qp];

   //Calculate elastic stress
  _elastic_stress_val[_qp] = _mat_const[_qp] * _elastic_strain_val[_qp];

  //Calculate strain energy
  _strain_energy_val[_qp]  = 0.5 * _elastic_strain_val[_qp] * _elastic_stress_val[_qp];
}
