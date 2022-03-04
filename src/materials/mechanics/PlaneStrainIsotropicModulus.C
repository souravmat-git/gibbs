//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* Calculates the elastic constants for a given youngs modulus and poissons

#include "PlaneStrainIsotropicModulus.h"
registerMooseObject("gibbsApp", PlaneStrainIsotropicModulus);

//template <>
InputParameters
PlaneStrainIsotropicModulus::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("youngs_modulus", "Youngs modulus");
  params.addRequiredParam<MaterialPropertyName>("poisson_ratio", "Poisson ration");
  params.addRequiredParam<std::string>("phase_name", "Phase name of the material");
  return params;
}

PlaneStrainIsotropicModulus::PlaneStrainIsotropicModulus(const InputParameters & parameters)
  : Material(parameters),
    _phase_name(getParam<std::string>("phase_name")),
    //Input properties
    _E(getMaterialProperty<Real>("youngs_modulus")),
    _nu(getMaterialProperty<Real>("poisson_ratio")),
    //Declare the values of the material properties to compute
    _C11_val(declareProperty<Real>("C11_" + _phase_name)),
    _C12_val(declareProperty<Real>("C12_" + _phase_name)),
    _C22_val(declareProperty<Real>("C22_" + _phase_name)),
    _C66_val(declareProperty<Real>("C66_" + _phase_name))
{
}

void
PlaneStrainIsotropicModulus::computeQpProperties(){

  _C11_val[_qp] = _E[_qp]* (1.0 - _nu[_qp])/((1.0 + _nu[_qp])* (1.0 - 2.0*_nu[_qp]));
  _C22_val[_qp] = _E[_qp]* (1.0 - _nu[_qp])/((1.0 + _nu[_qp])* (1.0 - 2.0*_nu[_qp])) ;
  _C12_val[_qp] = (_nu[_qp]*_E[_qp])/((1.0 + _nu[_qp]) * (1.0 - 2.0*_nu[_qp]));
  _C66_val[_qp] = _E[_qp]/(2.0*(1.0 + _nu[_qp]));

}
