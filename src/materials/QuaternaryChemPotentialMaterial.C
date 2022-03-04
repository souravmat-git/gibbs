//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryChemPotentialMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", QuaternaryChemPotentialMaterial);

//template <>
InputParameters
QuaternaryChemPotentialMaterial::validParams()
{
  InputParameters params = TabulatedPhaseMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot", "Chemical potential of comp A");
  params.addRequiredParam<MaterialPropertyName>("AB_therm_factor", "Thermodynamic factor w.r.t AB");
  params.addRequiredParam<MaterialPropertyName>("AC_therm_factor", "Thermodynamic factor w.r.t AC");
  params.addRequiredParam<MaterialPropertyName>("AD_therm_factor", "Thermodynamic factor w.r.t AD");
  params.addCoupledVar("mol_fraction_B", "Mole fraction of the component B");
  params.addCoupledVar("mol_fraction_C", "Mole fraction of component C ");
  params.addCoupledVar("mol_fraction_D", "Mole fraction of component D ");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values");
  return params;
}

QuaternaryChemPotentialMaterial::QuaternaryChemPotentialMaterial(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
   _A_chem_pot_name(getParam<MaterialPropertyName>("A_chem_pot")),
   _AB_therm_factor_name(getParam<MaterialPropertyName>("AB_therm_factor")),
   _AC_therm_factor_name(getParam<MaterialPropertyName>("AC_therm_factor")),
   _AD_therm_factor_name(getParam<MaterialPropertyName>("AD_therm_factor")),
   _A_chem_pot_val(declareProperty<Real>(_A_chem_pot_name)),
   _AB_therm_factor_val(declareProperty<Real>(_AB_therm_factor_name)),
   _AC_therm_factor_val(declareProperty<Real>(_AC_therm_factor_name)),
   _AD_therm_factor_val(declareProperty<Real>(_AD_therm_factor_name)),
   _mol_fraction_B(coupledValue("mol_fraction_B")),
   _mol_fraction_C(coupledValue("mol_fraction_C")),
   _mol_fraction_D(coupledValue("mol_fraction_D")),
   _table_object(getUserObject<QuaternaryChemPotentialData>("table_object"))
{
}

void
QuaternaryChemPotentialMaterial::computeQpProperties()
{

  //Note that we expect the data to be in J/mol
  //For phase-field simulations, we need to convert the data to J/m^3;
  //And then, non-dimensionalize the energy density with a characteristic energy density _Ec;

  //return the chemical potential of comp A in non-dimensional form
  _A_chem_pot_val[_qp] =
  ((_table_object.A_chem_pot(_mol_fraction_B[_qp],_mol_fraction_C[_qp],_mol_fraction_D[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor w.r.t AB in non-dimensial form
  _AB_therm_factor_val[_qp] =
  ((_table_object.thermodynamic_factor_AB(_mol_fraction_B[_qp], _mol_fraction_C[_qp], _mol_fraction_D[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor w.r.t AC in non-dimensial form
  _AC_therm_factor_val[_qp] =
  ((_table_object.thermodynamic_factor_AC(_mol_fraction_B[_qp], _mol_fraction_C[_qp], _mol_fraction_D[_qp]))/_Vm)/_Ec;

  //return the thermodynamic factor w.r.t AD in non-dimensial form
  _AD_therm_factor_val[_qp] =
  ((_table_object.thermodynamic_factor_AD(_mol_fraction_B[_qp], _mol_fraction_C[_qp], _mol_fraction_D[_qp]))/_Vm)/_Ec;

}
