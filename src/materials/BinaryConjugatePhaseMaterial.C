//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryConjugatePhaseMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", BinaryConjugatePhaseMaterial);

template <>
InputParameters
validParams<BinaryConjugatePhaseMaterial>()
{
  InputParameters params = validParams<TabulatedPhaseMaterial>();
  params.addRequiredParam<MaterialPropertyName>("A_chem_pot", "Chemical potential of the dependent comp");
  params.addRequiredParam<MaterialPropertyName>("xB", "Mole fraction of comp B");
  params.addRequiredParam<MaterialPropertyName>("inv_B_tf", "Thermodyanmic factor");
  //params.addRequiredParam<MaterialPropertyName>("inv_B_td", "Third derivative");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

BinaryConjugatePhaseMaterial::BinaryConjugatePhaseMaterial(const InputParameters & parameters)
  :TabulatedPhaseMaterial(parameters),
  _A_chem_pot_name(getParam<MaterialPropertyName>("A_chem_pot")),
  _xB_name(getParam<MaterialPropertyName>("xB")),
  _inv_B_tf_name(getParam<MaterialPropertyName>("inv_B_tf")),
  //_inv_B_td_name(getParam<MaterialPropertyName>("inv_B_td")),
  _A_chem_pot_val(declareProperty<Real>(_A_chem_pot_name)),
  _xB_val(declareProperty<Real>(_xB_name)),
  _inv_B_tf_val(declareProperty<Real>(_inv_B_tf_name)),
  //_inv_B_td_val(declareProperty<Real>(_inv_B_td_name)),
  _B_diff_pot(coupledValue("B_diff_pot")),
  _table_object(getUserObject<BinaryConjugatePhaseData>("table_object"))
{
}

void 
BinaryConjugatePhaseMaterial::computeQpProperties()
{
  //use a user object that returns the Chemical potential of the dependent component
  _A_chem_pot_val[_qp] = _table_object.A_chem_pot(_B_diff_pot[_qp])/_Vm;
  
  //return the mole fraction of B in non-dimensinla form
  _xB_val[_qp] = _table_object.xB(_B_diff_pot[_qp]);

  //return the inverse of the thermodynamic factor of B
  _inv_B_tf_val[_qp] = _table_object.inv_therm_factor_B(_B_diff_pot[_qp]);
  
  //return the inverse of the third derivative of B 
  //_inv_B_td_val[_qp] = _table_object.inv_third_deriv_B(_B_diff_pot[_qp]);
}
