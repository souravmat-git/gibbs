//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryKineticMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", TernaryKineticMaterial);

template <>
InputParameters
validParams<TernaryKineticMaterial>()
{
  InputParameters params = validParams<TabulatedKineticMaterial>();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("L_BC", "Onsager mobility L_BC");
  params.addRequiredParam<MaterialPropertyName>("L_CC", "Onsager mobility L_CC");
  params.addParam<MaterialPropertyName>("dL_BB_xB", 0.0, "L_BB w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BC_xB",0.0, "L_BC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CC_xB",0.0, "L_CC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BB_xC", 0.0, "L_BB w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BC_xC", 0.0, "L_BC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CC_xC", 0.0, "L_CC w.r.t C");
  params.addCoupledVar("mol_fraction_B", "Mole fraction of the component B");
  params.addCoupledVar("mol_fraction_C", "Mole fraction of component C ");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryKineticMaterial::TernaryKineticMaterial(const InputParameters & parameters)
  : TabulatedKineticMaterial(parameters),
    _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
    _L_BC_name(getParam<MaterialPropertyName>("L_BC")),
    _L_CC_name(getParam<MaterialPropertyName>("L_CC")),
    _dL_BB_xB_name(getParam<MaterialPropertyName>("dL_BB_xB")),
    _dL_BC_xB_name(getParam<MaterialPropertyName>("dL_BC_xB")),
    _dL_CC_xB_name(getParam<MaterialPropertyName>("dL_CC_xB")),
    _dL_BB_xC_name(getParam<MaterialPropertyName>("dL_BB_xC")),
    _dL_BC_xC_name(getParam<MaterialPropertyName>("dL_BC_xC")),
    _dL_CC_xC_name(getParam<MaterialPropertyName>("dL_CC_xC")),
    //Declare properties 
    _L_BB_val(declareProperty<Real>(_L_BB_name)),
    _L_BC_val(declareProperty<Real>(_L_BC_name)),
    _L_CC_val(declareProperty<Real>(_L_CC_name)),
    //Derivatives wr.t xB, xC
    _dL_BB_xB_val(declareProperty<Real>(_dL_BB_xB_name)),
    _dL_BB_xC_val(declareProperty<Real>(_dL_BC_xB_name)),
    //Derivatives wr,rt xB,xC
    _dL_BC_xB_val(declareProperty<Real>(_dL_CC_xB_name)),
    _dL_BC_xC_val(declareProperty<Real>(_dL_BB_xC_name)),
    //L_cc wr.r xB,xC
    _dL_CC_xB_val(declareProperty<Real>(_dL_BC_xC_name)),
    _dL_CC_xC_val(declareProperty<Real>(_dL_CC_xC_name)),
    //Each property is coupled to
    _mol_fraction_B(coupledValue("mol_fraction_B")),
    _mol_fraction_C(coupledValue("mol_fraction_C")),
   _table_object(getUserObject<TernaryMobilityData>("table_object"))
{
}

void 
TernaryKineticMaterial::computeQpProperties()
{
    //Diagonal terms of the matrix
    _L_BB_val[_qp] = _table_object.L_BB(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    _L_BC_val[_qp] = _table_object.L_BC(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    _L_CC_val[_qp] = _table_object.L_CC(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    
    
    _dL_BB_xB_val[_qp] = _table_object.dL_BB_xB(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    _dL_BC_xB_val[_qp] = _table_object.dL_BC_xB(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    _dL_CC_xB_val[_qp] = _table_object.dL_CC_xB(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    
      
    _dL_BB_xC_val[_qp] = _table_object.dL_BB_xC(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    _dL_BC_xC_val[_qp] = _table_object.dL_BC_xC(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
    _dL_CC_xC_val[_qp] = _table_object.dL_CC_xC(_mol_fraction_B[_qp],_mol_fraction_C[_qp]);
}
