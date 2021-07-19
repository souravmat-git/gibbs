//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "CoherentDrivingForce.h"

registerMooseObject("gibbsApp", CoherentDrivingForce);

template<>
InputParameters
validParams<CoherentDrivingForce>()
{
  InputParameters params = validParams<AuxKernel>();

  return params;
}

CoherentDrivingForce::CoherentDrivingForce(const InputParameters & parameters)
  :AuxKernel(parameters),
  _dh(getMaterialProperty<Real>("dh")),
  _comp_energy_alpha(getMaterialProperty<Real>("comp_energy_alpha")),
  _comp_energy_beta(getMaterialProperty<Real>("comp_energy_beta")),
  _A_chem_pot_alpha(getMaterialProperty<Real>("A_chem_pot_alpha")),
  _A_chem_pot_beta(getMaterialProperty<Real>("A_chem_pot_beta"))
{
}
 
Real
CoherentDrivingForce::computeValue()
{
  
  //This is the driving force due to chemical and mechanical energy
  return _dh[_qp]*((_A_chem_pot_beta[_qp]  + _comp_energy_beta[_qp]) 
                 - (_A_chem_pot_alpha[_qp] + _comp_energy_alpha[_qp]));
}
