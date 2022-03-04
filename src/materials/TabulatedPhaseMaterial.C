//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TabulatedPhaseMaterial.h"
//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", TabulatedPhaseMaterial);

//template <>
InputParameters
TabulatedPhaseMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<Real>("molar_volume", 1.0, "Molar volume of phase (m^{3})");
  params.addParam<Real>("char_energy", 1e9, "Characteristic energy J/mol");
  params.addClassDescription("This material class is the base class..."
                              "that returns the tabulated properties");
  return params;
}

TabulatedPhaseMaterial::TabulatedPhaseMaterial(const InputParameters & parameters)
  : Material(parameters),
   _Vm(getParam<Real>("molar_volume")),
   _Ec(getParam<Real>("char_energy"))
{
}
