//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThermoChemicalProperties.h"

registerMooseObject("gibbsApp",ThermoChemicalProperties);

template <>
InputParameters
validParams<ThermoChemicalProperties>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription("Base class from which binary, ternary..."
                              "and quartenary data can be extracted"); 
  return params;
}

ThermoChemicalProperties::ThermoChemicalProperties(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}

//This is a destructor
ThermoChemicalProperties::~ThermoChemicalProperties() 
{}

