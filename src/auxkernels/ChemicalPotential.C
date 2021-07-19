//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalPotential.h"

registerMooseObject("gibbsApp",ChemicalPotential);

template<>
InputParameters
validParams<ChemicalPotential>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<MaterialPropertyName>("dh", "interpolation");
  return params;
}

ChemicalPotential::ChemicalPotential(const InputParameters & parameters)
  :AuxKernel(parameters),
  _A_chem_pot_alpha(getMaterialProperty<Real>("A_chem_pot_alpha")),
  _A_chem_pot_beta(getMaterialProperty<Real>("A_chem_pot_beta")),
  _dh(getMaterialProperty<Real>("dh")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}
 
Real
ChemicalPotential::computeValue()
{                             
  return _nd_factor[_qp] * (_A_chem_pot_beta[_qp] - _A_chem_pot_alpha[_qp])*_dh[_qp];
}
