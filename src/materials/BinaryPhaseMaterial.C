//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryPhaseMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", BinaryPhaseMaterial);

template <>
InputParameters
validParams<BinaryPhaseMaterial>()
{
  InputParameters params = validParams<TabulatedPhaseMaterial>();
  params.addRequiredParam<MaterialPropertyName>("free_energy", "Free energy of a given phase");
  params.addRequiredParam<MaterialPropertyName>("B_diff_pot", "Diffusion potential of comp B");
  params.addRequiredParam<MaterialPropertyName>("B_therm_factor","Thermodyanmic factor");
  params.addRequiredCoupledVar("mol_fraction_B", "Mole fraction of the solute component");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

BinaryPhaseMaterial::BinaryPhaseMaterial(const InputParameters & parameters)
  :TabulatedPhaseMaterial(parameters),
  _f_phase_name(getParam<MaterialPropertyName>("free_energy")),
  _B_diff_pot_phase_name(getParam<MaterialPropertyName>("B_diff_pot")),
  _B_therm_factor_phase_name(getParam<MaterialPropertyName>("B_therm_factor")),
  _f_phase_val(declareProperty<Real>(_f_phase_name)),
  _B_diff_pot_phase_val(declareProperty<Real>(_B_diff_pot_phase_name)),
  _B_therm_factor_phase_val(declareProperty<Real>(_B_therm_factor_phase_name)),
  _mol_fraction_B(coupledValue("mol_fraction_B")),
  _table_object(getUserObject<BinaryPhaseData>("table_object"))
{
}

void 
BinaryPhaseMaterial::computeQpProperties()
{

  //Note that we expect the data to be in J/mol
  //For phase-field simulations, we need to convert the data to J/m^3;
  //And then, non-dimensionalize the energy with a characteristic energy _Ec;

  //use a user object that returns the free energy in non-dimensional form
  _f_phase_val[_qp] = ((_table_object.free_energy(_mol_fraction_B[_qp]))/_Vm)/_Ec;
  
  //return the diffusion potential of comp B in non-dimensional form
  _B_diff_pot_phase_val[_qp] = ((_table_object.B_diff_pot(_mol_fraction_B[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor in non-dimensial form
  _B_therm_factor_phase_val[_qp] =((_table_object.thermodynamic_factor(_mol_fraction_B[_qp]))/_Vm)/_Ec;
}
