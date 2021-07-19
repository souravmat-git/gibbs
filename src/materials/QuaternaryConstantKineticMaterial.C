//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryConstantKineticMaterial.h"
registerMooseObject("gibbsApp", QuaternaryConstantKineticMaterial);

template <>
InputParameters
validParams<QuaternaryConstantKineticMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("L_CC", "Onsager mobility L_CC");
  params.addRequiredParam<MaterialPropertyName>("L_DD", "Onsager mobility L_DD");
  params.addRequiredParam<MaterialPropertyName>("L_BC", "Onsager mobility L_BC");
  params.addRequiredParam<MaterialPropertyName>("L_BD", "Onsager mobility L_BD");
  params.addRequiredParam<MaterialPropertyName>("L_CD", "Onsager mobility L_CD");
  //w.r.t B
  params.addParam<MaterialPropertyName>("dL_BB_muB",0.0,"L_BB w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CC_muB",0.0,"L_CC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_DD_muB",0.0,"L_DD w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BC_muB",0.0,"L_BC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BD_muB",0.0,"L_BD w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CD_muB",0.0,"L_CD w.r.t B");
  //w.r.t C
  params.addParam<MaterialPropertyName>("dL_BB_muC",0.0,"L_BB w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CC_muC",0.0,"L_CC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_DD_muC",0.0,"L_DD w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BC_muC",0.0,"L_BC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BD_muC",0.0,"L_BD w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CD_muC",0.0,"L_CD w.r.t C");
  //w.r.t D
  params.addParam<MaterialPropertyName>("dL_BB_muD",0.0,"L_BB w.r.t D");
  params.addParam<MaterialPropertyName>("dL_CC_muD",0.0,"L_CC w.r.t D");
  params.addParam<MaterialPropertyName>("dL_DD_muD",0.0,"L_DD w.r.t D");
  params.addParam<MaterialPropertyName>("dL_BC_muD",0.0,"L_BC w.r.t D");
  params.addParam<MaterialPropertyName>("dL_BD_muD",0.0,"L_BD w.r.t D");
  params.addParam<MaterialPropertyName>("dL_CD_muD",0.0,"L_CD w.r.t D");
  //Input values
  params.addRequiredParam<Real>("L_BB_val", "L_BB value");
  params.addRequiredParam<Real>("L_CC_val", "L_CCue");
  params.addRequiredParam<Real>("L_DD_val", "L_DDue" );
  params.addRequiredParam<Real>("L_BC_val", "L_BC value");
  params.addRequiredParam<Real>("L_BD_val", "L_BDue");
  params.addRequiredParam<Real>("L_CD_val", "L_CDue" );
  params.addClassDescription("Constant Onsager matrix for A-B-C-D alloy"); 
  return params;
}

QuaternaryConstantKineticMaterial::QuaternaryConstantKineticMaterial(const InputParameters & parameters)
  : Material(parameters),
   _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
   _L_CC_name(getParam<MaterialPropertyName>("L_CC")),
   _L_DD_name(getParam<MaterialPropertyName>("L_DD")),
   _L_BC_name(getParam<MaterialPropertyName>("L_BC")),
   _L_BD_name(getParam<MaterialPropertyName>("L_BD")),
   _L_CD_name(getParam<MaterialPropertyName>("L_CD")),
   //Derivative with respect to muB   
   _dL_BB_muB_name(getParam<MaterialPropertyName>("dL_BB_muB")),
   _dL_CC_muB_name(getParam<MaterialPropertyName>("dL_CC_muB")),
   _dL_DD_muB_name(getParam<MaterialPropertyName>("dL_DD_muB")),
   _dL_BC_muB_name(getParam<MaterialPropertyName>("dL_BC_muB")),
   _dL_BD_muB_name(getParam<MaterialPropertyName>("dL_BD_muB")),
   _dL_CD_muB_name(getParam<MaterialPropertyName>("dL_CD_muB")),
    //Derivative with respect to muC  
   _dL_BB_muC_name(getParam<MaterialPropertyName>("dL_BB_muC")),
   _dL_CC_muC_name(getParam<MaterialPropertyName>("dL_CC_muC")),
   _dL_DD_muC_name(getParam<MaterialPropertyName>("dL_DD_muC")),
   _dL_BC_muC_name(getParam<MaterialPropertyName>("dL_BC_muC")),
   _dL_BD_muC_name(getParam<MaterialPropertyName>("dL_BD_muC")),
   _dL_CD_muC_name(getParam<MaterialPropertyName>("dL_CD_muC")),
   //Derivative with respect to muD  
   _dL_BB_muD_name(getParam<MaterialPropertyName>("dL_BB_muD")),
   _dL_CC_muD_name(getParam<MaterialPropertyName>("dL_CC_muD")),
   _dL_DD_muD_name(getParam<MaterialPropertyName>("dL_DD_muD")),
   _dL_BC_muD_name(getParam<MaterialPropertyName>("dL_BC_muD")),
   _dL_BD_muD_name(getParam<MaterialPropertyName>("dL_BD_muD")),
   _dL_CD_muD_name(getParam<MaterialPropertyName>("dL_CD_muD")),
    //Declare properties 
   _L_BB(declareProperty<Real>(_L_BB_name)),
   _L_CC(declareProperty<Real>(_L_CC_name)),
   _L_DD(declareProperty<Real>(_L_DD_name)),
   _L_BC(declareProperty<Real>(_L_BC_name)),
   _L_BD(declareProperty<Real>(_L_BD_name)),
   _L_CD(declareProperty<Real>(_L_CD_name)),
    //Derivatives wr.t muB
    _dL_BB_muB(declareProperty<Real>(_dL_BB_muB_name)),
    _dL_CC_muB(declareProperty<Real>(_dL_CC_muB_name)),
    _dL_DD_muB(declareProperty<Real>(_dL_DD_muB_name)),
    _dL_BC_muB(declareProperty<Real>(_dL_BC_muB_name)),
    _dL_BD_muB(declareProperty<Real>(_dL_BD_muB_name)),
    _dL_CD_muB(declareProperty<Real>(_dL_CD_muB_name)),
    //Derivatives wr,rt xB,xC
    _dL_BB_muC(declareProperty<Real>(_dL_BB_muC_name)),
    _dL_CC_muC(declareProperty<Real>(_dL_CC_muC_name)),
    _dL_DD_muC(declareProperty<Real>(_dL_DD_muC_name)),
    _dL_BC_muC(declareProperty<Real>(_dL_BC_muC_name)),
    _dL_BD_muC(declareProperty<Real>(_dL_BD_muC_name)),
    _dL_CD_muC(declareProperty<Real>(_dL_CD_muC_name)),
    //L_cc wr.r xB,xC
    _dL_BB_muD(declareProperty<Real>(_dL_BB_muD_name)),
    _dL_CC_muD(declareProperty<Real>(_dL_CC_muD_name)),
    _dL_DD_muD(declareProperty<Real>(_dL_DD_muD_name)),
    _dL_BC_muD(declareProperty<Real>(_dL_BC_muD_name)),
    _dL_BD_muD(declareProperty<Real>(_dL_BD_muD_name)),
    _dL_CD_muD(declareProperty<Real>(_dL_CD_muD_name)),
    //Input values
    _L_BB_val(getParam<Real>("L_BB_val")),
    _L_CC_val(getParam<Real>("L_CC_val")),
    _L_DD_val(getParam<Real>("L_DD_val")),
    _L_BC_val(getParam<Real>("L_BC_val")),
    _L_BD_val(getParam<Real>("L_BD_val")),
    _L_CD_val(getParam<Real>("L_CD_val"))
{
}

void 
QuaternaryConstantKineticMaterial::computeQpProperties()
{
    //Onsager matrix for a A-B-C-D alloy
    _L_BB[_qp] = _L_BB_val;
    _L_CC[_qp] = _L_CC_val;
    _L_DD[_qp] = _L_DD_val;
    _L_BC[_qp] = _L_BC_val;
    _L_BD[_qp] = _L_BD_val;
    _L_CD[_qp] = _L_CD_val;
    
    //Derivative with respect to muB
    _dL_BB_muB[_qp] = 0;
    _dL_CC_muB[_qp] = 0;   
    _dL_DD_muB[_qp] = 0;
    
    _dL_BC_muB[_qp] = 0;
    _dL_BD_muB[_qp] = 0;   
    _dL_CD_muB[_qp] = 0;
        
    //Derivative with respect to muC
    _dL_BB_muC[_qp] = 0;
    _dL_CC_muC[_qp] = 0;   
    _dL_DD_muC[_qp] = 0;
    
    _dL_BC_muC[_qp] = 0;
    _dL_BD_muC[_qp] = 0;   
    _dL_CD_muC[_qp] = 0;
    
    //Derivative with respect to muD
    _dL_BB_muD[_qp] = 0;
    _dL_CC_muD[_qp] = 0;   
    _dL_DD_muD[_qp] = 0;
    
    _dL_BC_muD[_qp] = 0;
    _dL_BD_muD[_qp] = 0;   
    _dL_CD_muD[_qp] = 0;
}
