//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryKineticMaterial.h"
//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", BinaryKineticMaterial);

//template <>
InputParameters
BinaryKineticMaterial::validParams()
{
  InputParameters params = TabulatedKineticMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addParam<MaterialPropertyName>("dL_BB_xB", 0.0, "L_BB w.r.t mole fraction of B");
  params.addRequiredCoupledVar("mol_fraction_B", "Mole fraction of comp B");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values");
  return params;
}

BinaryKineticMaterial::BinaryKineticMaterial(const InputParameters & parameters)
  :TabulatedKineticMaterial(parameters),
  _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
  _dL_BB_xB_name(getParam<MaterialPropertyName>("dL_BB_xB")),
  _L_BB_val(declareProperty<Real>(_L_BB_name)),
  _dL_BB_xB_val(declareProperty<Real>(_dL_BB_xB_name)),
  _xB(coupledValue("mol_fraction_B")),
  _table_object(getUserObject<BinaryMobilityData>("table_object"))
{
}

void
BinaryKineticMaterial::computeQpProperties()
{
  //use a user object that returns the Onsager mobility as a function of composition
  _L_BB_val[_qp] = _table_object.L_BB(_xB[_qp]);

  //return the interpolated first derivative wrt mole fraction of B
  _dL_BB_xB_val[_qp] = _table_object.dL_BB_xB(_xB[_qp]);
}
