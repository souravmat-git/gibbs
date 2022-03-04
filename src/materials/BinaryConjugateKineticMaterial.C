//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryConjugateKineticMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", BinaryConjugateKineticMaterial);

//template <>
InputParameters
BinaryConjugateKineticMaterial::validParams()
{
  InputParameters params = TabulatedKineticMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addParam<MaterialPropertyName>("dL_BB_muB", 0.0, "L_BB w.r.t B");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values");
  return params;
}

BinaryConjugateKineticMaterial::BinaryConjugateKineticMaterial(const InputParameters & parameters)
  :TabulatedKineticMaterial(parameters),
  _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
  _dL_BB_muB_name(getParam<MaterialPropertyName>("dL_BB_muB")),
  _L_BB_val(declareProperty<Real>(_L_BB_name)),
  _dL_BB_muB_val(declareProperty<Real>(_dL_BB_muB_name)),
  _B_diff_pot(coupledValue("B_diff_pot")),
  _table_object(getUserObject<BinaryConjugateMobilityData>("table_object"))
{
}

void
BinaryConjugateKineticMaterial::computeQpProperties()
{
  //use a user object that returns the Chemical potential of the dep in non-dimensional form
  _L_BB_val[_qp] = _table_object.L_BB(_B_diff_pot[_qp]);

  //return the mole fraction of B in non-dimensinla form
  _dL_BB_muB_val[_qp] = _table_object.dL_BB_muB(_B_diff_pot[_qp]);
}
