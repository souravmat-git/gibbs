//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TabulatedKineticMaterial.h"
//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", TabulatedKineticMaterial);

//template <>
InputParameters
TabulatedKineticMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<Real>("char_time", 1e-9, "Characteristic time of the sim s");
  params.addParam<Real>("char_length",1e-6, "Characteristic length of the sim m");
  params.addParam<Real>("char_energy", 1e4, "Characteristic energy J/mol");
  params.addParam<Real>("molar_volume", 1.0, "Molar volume of phase (m^{3})");
  params.addClassDescription("This material class is the base class..."
                              "that returns the tabulated properties");
  return params;
}

TabulatedKineticMaterial::TabulatedKineticMaterial(const InputParameters & parameters)
  : Material(parameters),
   _tc(getParam<Real>("char_time")),
   _xc(getParam<Real>("char_length")),
   _Ec(getParam<Real>("char_energy")),
   _Vm(getParam<Real>("molar_volume"))
{
}
