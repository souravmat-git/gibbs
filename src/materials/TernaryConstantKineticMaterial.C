//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryConstantKineticMaterial.h"
registerMooseObject("gibbsApp", TernaryConstantKineticMaterial);

template <>
InputParameters
validParams<TernaryConstantKineticMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("L_BC", "Onsager mobility L_BC");
  params.addRequiredParam<MaterialPropertyName>("L_CC", "Onsager mobility L_CC");
  params.addRequiredParam<MaterialPropertyName>("dL_BB_muB","L_BB w.r.t B");
  params.addRequiredParam<MaterialPropertyName>("dL_BC_muB","L_BC w.r.t B");
  params.addRequiredParam<MaterialPropertyName>("dL_CC_muB","L_CC w.r.t C");
  params.addRequiredParam<MaterialPropertyName>("dL_BB_muC","L_BB w.r.t C");
  params.addRequiredParam<MaterialPropertyName>("dL_BC_muC","L_BC w.r.t C");
  params.addRequiredParam<MaterialPropertyName>("dL_CC_muC", "L_CC w.r.t C");
  params.addRequiredParam<Real>("L_BB_val", "L_BB value");
  params.addRequiredParam<Real>("L_BC_val", "L_BC_value");
  params.addRequiredParam<Real>("L_CC_val", "L_CC_value" );
  params.addClassDescription("Constant Onsager matrix for A-B-C alloy"); 
  return params;
}

TernaryConstantKineticMaterial::TernaryConstantKineticMaterial(const InputParameters & parameters)
  : Material(parameters),
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
    _L_BB(declareProperty<Real>(_L_BB_name)),
    _L_BC(declareProperty<Real>(_L_BC_name)),
    _L_CC(declareProperty<Real>(_L_CC_name)),
    //Derivatives wr.t xB, xC
    _dL_BB_muB(declareProperty<Real>(_dL_BB_muB_name)),
    _dL_BB_muC(declareProperty<Real>(_dL_BC_muB_name)),
    //Derivatives wr,rt xB,xC
    _dL_BC_muB(declareProperty<Real>(_dL_CC_muB_name)),
    _dL_BC_muC(declareProperty<Real>(_dL_BB_muC_name)),
    //L_cc wr.r xB,xC
    _dL_CC_muB(declareProperty<Real>(_dL_BC_muC_name)),
    _dL_CC_muC(declareProperty<Real>(_dL_CC_muC_name)),
    _L_BB_val(getParam<Real>("L_BB_val")),
    _L_BC_val(getParam<Real>("L_BC_val")),
    _L_CC_val(getParam<Real>("L_CC_val"))
{
}

void 
TernaryConstantKineticMaterial::computeQpProperties()
{
    //Onsager matrix for a A-B-C alloy
    _L_BB[_qp] = _L_BB_val;
    _L_BC[_qp] = _L_BC_val;
    _L_CC[_qp] = _L_CC_val;
    
    
    //Derivative with respect to muB
    _dL_BB_muB[_qp] = 0;
    _dL_BC_muB[_qp] = 0;   
    _dL_CC_muB[_qp] = 0;
    //Derivative with respect to muC  
    _dL_BB_muC[_qp] = 0; 
    _dL_BC_muC[_qp] = 0;
    _dL_CC_muC[_qp] = 0;
}
