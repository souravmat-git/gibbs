//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryConjugatePhaseMaterial.h"
//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", TernaryConjugatePhaseMaterial);

//template <>
InputParameters
TernaryConjugatePhaseMaterial::validParams()
{
  InputParameters params = TabulatedPhaseMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot", "Chemical potential of the depe");
  params.addRequiredParam<MaterialPropertyName>("B_mole_fraction", "Mole fraction of comp B");
  params.addRequiredParam<MaterialPropertyName>("C_mole_fraction", "Mole fraction of comp C");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf", "Inverse of thermodynamic factor w.r.t B");
  params.addRequiredParam<MaterialPropertyName>("inv_BC_tf", "Inverse of thermodynamic factor w.r.t BC");
  params.addRequiredParam<MaterialPropertyName>("inv_C_tf", "Inverse of thermodynamic factor w.r.t C");
  params.addCoupledVar("B_diff_pot", "Diffusion potential  of component B");
  params.addCoupledVar("C_diff_pot", "Diffusion potential  of component C ");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values");
  return params;
}

TernaryConjugatePhaseMaterial::TernaryConjugatePhaseMaterial(const InputParameters & parameters)
  : TabulatedPhaseMaterial(parameters),
    _A_chem_pot_name(getParam<MaterialPropertyName>("A_chem_pot")),
    _xB_name(getParam<MaterialPropertyName>("B_mole_fraction")),
    _xC_name(getParam<MaterialPropertyName>("C_mole_fraction")),
    _inv_B_tf_name(getParam<MaterialPropertyName>("inv_B_tf")),
    _inv_BC_tf_name(getParam<MaterialPropertyName>("inv_BC_tf")),
    _inv_C_tf_name(getParam<MaterialPropertyName>("inv_C_tf")),
    //Assign the material property name to material property value
    _A_chem_pot_val(declareProperty<Real>(_A_chem_pot_name)),
    _xB_val(declareProperty<Real>(_xB_name)),
    _xC_val(declareProperty<Real>(_xC_name)),
    _inv_B_tf_val(declareProperty<Real>(_inv_B_tf_name)),
    _inv_BC_tf_val(declareProperty<Real>(_inv_BC_tf_name)),
    _inv_C_tf_val(declareProperty<Real>(_inv_C_tf_name)),
    //Coupled variables
    _B_diff_pot(coupledValue("B_diff_pot")),
    _C_diff_pot(coupledValue("C_diff_pot")),
   _table_object(getUserObject<TernaryConjugatePhaseData>("table_object"))
{
}

void
TernaryConjugatePhaseMaterial::computeQpProperties()
{

  //Note that we expect the data to be in non-dimensional form

  _A_chem_pot_val[_qp] = _table_object.A_chem_pot(_B_diff_pot[_qp],_C_diff_pot[_qp]);

  //return the mole fraction of comp B in non-dimensional form
  _xB_val[_qp] =     _table_object.xB(_B_diff_pot[_qp],_C_diff_pot[_qp]);

   //turn the mole fraction of comp C in non-dimensional form
  _xC_val[_qp] =    _table_object.xC(_B_diff_pot[_qp],_C_diff_pot[_qp]);

  //return the thermodynamic factor w.r.t B in non-dimensial form
  _inv_B_tf_val[_qp] = _table_object.inv_therm_factor_B(_B_diff_pot[_qp], _C_diff_pot[_qp]);

  //return the thermodynamic factor w.r.t BC in non-dimensial form
  _inv_BC_tf_val[_qp] = _table_object.inv_therm_factor_BC(_B_diff_pot[_qp], _C_diff_pot[_qp]);

  //return the thermodynamic factor w.r.t BC in non-dimensial form
  _inv_C_tf_val[_qp] = _table_object.inv_therm_factor_C(_B_diff_pot[_qp], _C_diff_pot[_qp]);
}
