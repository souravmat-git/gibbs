//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryConjugateKineticMaterial.h"

//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", TernaryConjugateKineticMaterial);

template <>
InputParameters
validParams<TernaryConjugateKineticMaterial>()
{
  InputParameters params = validParams<TabulatedKineticMaterial>();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("L_BC", "Onsager mobility L_BC");
  params.addRequiredParam<MaterialPropertyName>("L_CC", "Onsager mobility L_CC");
  params.addParam<MaterialPropertyName>("dL_BB_muB", 0.0, "L_BB w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BC_muB",0.0, "L_BC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CC_muB",0.0, "L_CC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BB_muC", 0.0, "L_BB w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BC_muC", 0.0, "L_BC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CC_muC", 0.0, "L_CC w.r.t C");
  params.addCoupledVar("B_diff_pot", "Diffusion potential of component B");
  params.addCoupledVar("C_diff_pot", "Diffusion potential of component C");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryConjugateKineticMaterial::TernaryConjugateKineticMaterial(const InputParameters & parameters)
  : TabulatedKineticMaterial(parameters),
    _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
    _L_BC_name(getParam<MaterialPropertyName>("L_BC")),
    _L_CC_name(getParam<MaterialPropertyName>("L_CC")),
    _dL_BB_muB_name(getParam<MaterialPropertyName>("dL_BB_muB")),
    _dL_BC_muB_name(getParam<MaterialPropertyName>("dL_BC_muB")),
    _dL_CC_muB_name(getParam<MaterialPropertyName>("dL_CC_muB")),
    _dL_BB_muC_name(getParam<MaterialPropertyName>("dL_BB_muC")),
    _dL_BC_muC_name(getParam<MaterialPropertyName>("dL_BC_muC")),
    _dL_CC_muC_name(getParam<MaterialPropertyName>("dL_CC_muC")),
    //Declare properties 
    _L_BB_val(declareProperty<Real>(_L_BB_name)),
    _L_BC_val(declareProperty<Real>(_L_BC_name)),
    _L_CC_val(declareProperty<Real>(_L_CC_name)),
    //Derivatives wr.t xB, xC
    _dL_BB_muB_val(declareProperty<Real>(_dL_BB_muB_name)),
    _dL_BB_muC_val(declareProperty<Real>(_dL_BC_muB_name)),
    //Derivatives wr,rt xB,xC
    _dL_BC_muB_val(declareProperty<Real>(_dL_CC_muB_name)),
    _dL_BC_muC_val(declareProperty<Real>(_dL_BB_muC_name)),
    //L_cc wr.r xB,xC
    _dL_CC_muB_val(declareProperty<Real>(_dL_BC_muC_name)),
    _dL_CC_muC_val(declareProperty<Real>(_dL_CC_muC_name)),
    //Each property is coupled to
    _B_diff_pot(coupledValue("B_diff_pot")),
    _C_diff_pot(coupledValue("C_diff_pot")),
   _table_object(getUserObject<TernaryConjugateMobilityData>("table_object"))
{
}

void 
TernaryConjugateKineticMaterial::computeQpProperties()
{
    //Diagonal terms of the matrix
    _L_BB_val[_qp] = _table_object.L_BB(_B_diff_pot[_qp], _C_diff_pot[_qp]);
    _L_BC_val[_qp] = _table_object.L_BC(_B_diff_pot[_qp], _C_diff_pot[_qp]);
    _L_CC_val[_qp] = _table_object.L_CC(_B_diff_pot[_qp], _C_diff_pot[_qp]);
    
    
    //Derivative with respect to muB
    _dL_BB_muB_val[_qp] = _table_object.dL_BB_muB(_B_diff_pot[_qp], _C_diff_pot[_qp]);
    _dL_BC_muB_val[_qp] = _table_object.dL_BC_muB(_B_diff_pot[_qp], _C_diff_pot[_qp]);   
    _dL_CC_muB_val[_qp] = _table_object.dL_CC_muB(_B_diff_pot[_qp], _C_diff_pot[_qp]);
    
      
    _dL_BB_muC_val[_qp] = _table_object.dL_BB_muC(_B_diff_pot[_qp], _C_diff_pot[_qp]); 
    _dL_BC_muC_val[_qp] = _table_object.dL_BC_muC(_B_diff_pot[_qp], _C_diff_pot[_qp]);
    _dL_CC_muC_val[_qp] = _table_object.dL_CC_muC(_B_diff_pot[_qp], _C_diff_pot[_qp]);
}
