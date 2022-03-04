//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryPhaseMaterial.h"
//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", TernaryPhaseMaterial);

//template <>
InputParameters
TernaryPhaseMaterial::validParams()
{
  InputParameters params = TabulatedPhaseMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("free_energy", "Free energy of a given phase");
  params.addRequiredParam<MaterialPropertyName>("B_diff_pot", "Diffusion potential of comp B");
  params.addRequiredParam<MaterialPropertyName>("C_diff_pot", "Diffusion potential of comp C");
  params.addRequiredParam<MaterialPropertyName>("B_therm_factor", "Thermodynamic factor w.r.t B");
  params.addRequiredParam<MaterialPropertyName>("BC_therm_factor", "Thermodynamic factor w.r.t BC");
  params.addRequiredParam<MaterialPropertyName>("C_therm_factor", "Thermodynamic factor w.r,t C");
  params.addCoupledVar("mol_fraction_B", "Mole fraction of the component B");
  params.addCoupledVar("mol_fraction_C", "Mole fraction of component C ");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values");
  return params;
}

TernaryPhaseMaterial::TernaryPhaseMaterial(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
    _f_energy_name(getParam<MaterialPropertyName>("free_energy")),
    _B_diff_pot_name(getParam<MaterialPropertyName>("B_diff_pot")),
    _C_diff_pot_name(getParam<MaterialPropertyName>("C_diff_pot")),
    _B_therm_factor_name(getParam<MaterialPropertyName>("B_therm_factor")),
    _BC_therm_factor_name(getParam<MaterialPropertyName>("BC_therm_factor")),
    _C_therm_factor_name(getParam<MaterialPropertyName>("C_therm_factor")),
    _f_energy_val(declareProperty<Real>(_f_energy_name)),
    _B_diff_pot_val(declareProperty<Real>(_B_diff_pot_name)),
    _C_diff_pot_val(declareProperty<Real>(_C_diff_pot_name)),
    _B_therm_factor_val(declareProperty<Real>(_B_therm_factor_name)),
    _BC_therm_factor_val(declareProperty<Real>(_BC_therm_factor_name)),
    _C_therm_factor_val(declareProperty<Real>(_C_therm_factor_name)),
    _mol_fraction_B(coupledValue("mol_fraction_B")),
    _mol_fraction_C(coupledValue("mol_fraction_C")),
   _table_object(getUserObject<TernaryPhaseData>("table_object"))
{
}

void
TernaryPhaseMaterial::computeQpProperties()
{

  //Note that we expect the data to be in J/mol
  //For phase-field simulations, we need to convert the data to J/m^3;
  //And then, non-dimensionalize the energy density with a characteristic energy density _Ec;

  //use a user object that returns the free energy in non-dimensional form
  _f_energy_val[_qp] = ((_table_object.free_energy(_mol_fraction_B[_qp],_mol_fraction_C[_qp]))/_Vm)/_Ec;

  //return the diffusion potential of comp B in non-dimensional form
  _B_diff_pot_val[_qp] = ((_table_object.B_diff_pot(_mol_fraction_B[_qp],_mol_fraction_C[_qp]))/_Vm)/_Ec;

   //return the diffusion potential of comp C in non-dimensional form
  _C_diff_pot_val[_qp] = ((_table_object.C_diff_pot(_mol_fraction_B[_qp], _mol_fraction_C[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor w.r.t B in non-dimensial form
  _B_therm_factor_val[_qp] =((_table_object.thermodynamic_factor_B(_mol_fraction_B[_qp], _mol_fraction_C[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor w.r.t BC in non-dimensial form
  _BC_therm_factor_val[_qp] =((_table_object.thermodynamic_factor_BC(_mol_fraction_B[_qp], _mol_fraction_C[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor w.r.t B in non-dimensial form
  _C_therm_factor_val[_qp] =((_table_object.thermodynamic_factor_C(_mol_fraction_B[_qp], _mol_fraction_C[_qp]))/_Vm)/_Ec;
}
