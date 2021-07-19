//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "MultiCompMultiPhaseBase.h"

registerMooseObject("gibbsApp", MultiCompMultiPhaseBase);

template <>
InputParameters
validParams<MultiCompMultiPhaseBase>()
{
  InputParameters params = validParams<MultiPhaseBase>();
  params.addClassDescription("Base class for continuity eqn for a A-B-C-D alloy");
  //B thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_B_tf_gamma",  0.0, "Thermodynamic factor B in gamma phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_delta",  0.0, "Thermodynamic factor B in delta phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_epsilon",0.0, "Thermodynamic factor B in epsilon phase");
  //BC thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_BC_tf_gamma", 0.0, "Thermodynamic factor BC in gamma phase");
  params.addParam<MaterialPropertyName>("inv_BC_tf_delta", 0.0, "Thermodynamic factor BC in delta phase");
  params.addParam<MaterialPropertyName>("inv_BC_tf_epsilon", 0.0, "Thermodynamic factor BC in epsilon phase");
  //BD thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_BD_tf_gamma", 0.0, "Thermodynamic factor BD in gamma phase");
  params.addParam<MaterialPropertyName>("inv_BD_tf_delta", 0.0, "Thermodynamic factor BD in delta phase");
  params.addParam<MaterialPropertyName>("inv_BD_tf_epsilon", 0.0, "Thermodynamic factor BD in epsilon phase");
  //C thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_C_tf_gamma", 0.0, "Thermodynamic factor C in gamma phase");
  params.addParam<MaterialPropertyName>("inv_C_tf_delta", 0.0, "Thermodynamic factor C in delta phase");
  params.addParam<MaterialPropertyName>("inv_C_tf_epsilon", 0.0, "Thermodynamic factor C in epsilon phase");
  //CD thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_CD_tf_gamma", 0.0, "Thermodynamic factor CD in gamma phase");
  params.addParam<MaterialPropertyName>("inv_CD_tf_delta", 0.0, "Thermodynamic factor CD in delta phase");
  params.addParam<MaterialPropertyName>("inv_CD_tf_epsilon",0.0, "Thermodynamic factor CD in epsilon phase");
  //BD thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_D_tf_gamma", 0.0, "Thermodynamic factor D in gamma phase");
  params.addParam<MaterialPropertyName>("inv_D_tf_delta", 0.0, "Thermodynamic factor D in delta phase");
  params.addParam<MaterialPropertyName>("inv_D_tf_epsilon", 0.0, "Thermodynamic factor D in epsilon phase");  
  //***************************************Onsager mobility assumed only for three phases***************//
  //L_BB
  params.addParam<MaterialPropertyName>("L_BB_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("L_BB_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("L_BB_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");  
  //L_BC
  params.addParam<MaterialPropertyName>("L_BC_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("L_BC_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("L_BC_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase");  
  //L_BD
  params.addParam<MaterialPropertyName>("L_BD_gamma", 0.0, "Onsager mobility of comp BD in gamma phase");
  params.addParam<MaterialPropertyName>("L_BD_delta", 0.0, "Onsager mobility of comp BD in delta phase");
  params.addParam<MaterialPropertyName>("L_BD_epsilon", 0.0, "TOnsager mobility of comp BD in epsilon phase");  
  //L_CC
  params.addParam<MaterialPropertyName>("L_CC_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("L_CC_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("L_CC_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase"); 
  //L_CD
  params.addParam<MaterialPropertyName>("L_CD_gamma", 0.0, "Onsager mobility of comp CD in gamma phase");
  params.addParam<MaterialPropertyName>("L_CD_delta", 0.0, "Onsager mobility of comp CD in delta phase");
  params.addParam<MaterialPropertyName>("L_CD_epsilon", 0.0, "TOnsager mobility of comp CD in epsilon phase");  
  //L_DD
  params.addParam<MaterialPropertyName>("L_DD_gamma", 0.0, "Onsager mobility of comp D in gamma phase");
  params.addParam<MaterialPropertyName>("L_DD_delta", 0.0, "Onsager mobility of comp D in delta phase");
  params.addParam<MaterialPropertyName>("L_DD_epsilon", 0.0, "Onsager mobility of comp D in epsilon phase"); 
  //*************************Derivatives of L wr.t B********************************************************//  
   //dL_BB_muB
  params.addParam<MaterialPropertyName>("dL_BB_muB_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BB_muB_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("dL_BB_muB_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");  
  //dL_BC_muB
  params.addParam<MaterialPropertyName>("dL_BC_muB_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BC_muB_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("dL_BC_muB_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase");  
  //dL_BD_muB
  params.addParam<MaterialPropertyName>("dL_BD_muB_gamma", 0.0, "Onsager mobility of comp BD in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BD_muB_delta", 0.0, "Onsager mobility of comp BD in delta phase");
  params.addParam<MaterialPropertyName>("dL_BD_muB_epsilon", 0.0, "TOnsager mobility of comp BD in epsilon phase");  
  //dL_CC_muB
  params.addParam<MaterialPropertyName>("dL_CC_muB_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CC_muB_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("dL_CC_muB_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase");  
  //dL_CD_muB
  params.addParam<MaterialPropertyName>("dL_CD_muB_gamma", 0.0, "Onsager mobility of comp CD in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CD_muB_delta", 0.0, "Onsager mobility of comp CD in delta phase");
  params.addParam<MaterialPropertyName>("dL_CD_muB_epsilon", 0.0, "TOnsager mobility of comp CD in epsilon phase");  
  //dL_DD_muB
  params.addParam<MaterialPropertyName>("dL_DD_muB_gamma", 0.0, "Onsager mobility of comp D in gamma phase");
  params.addParam<MaterialPropertyName>("dL_DD_muB_delta", 0.0, "Onsager mobility of comp D in delta phase");
  params.addParam<MaterialPropertyName>("dL_DD_muB_epsilon", 0.0, "Onsager mobility of comp D in epsilon phase");  
  //*************************Derivatives of L w.r.t muC********************************************************//  
   //dL_BB_muC
  params.addParam<MaterialPropertyName>("dL_BB_muC_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BB_muC_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("dL_BB_muC_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");  
  //dL_BC_muC
  params.addParam<MaterialPropertyName>("dL_BC_muC_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BC_muC_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("dL_BC_muC_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase"); 
  //dL_BD_muC
  params.addParam<MaterialPropertyName>("dL_BD_muC_gamma", 0.0, "Onsager mobility of comp BD in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BD_muC_delta", 0.0, "Onsager mobility of comp BD in delta phase");
  params.addParam<MaterialPropertyName>("dL_BD_muC_epsilon", 0.0, "TOnsager mobility of comp BD in epsilon phase");  
  //dL_CC_muC
  params.addParam<MaterialPropertyName>("dL_CC_muC_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CC_muC_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("dL_CC_muC_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase");  
  //dL_CD_muC
  params.addParam<MaterialPropertyName>("dL_CD_muC_gamma", 0.0, "Onsager mobility of comp CD in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CD_muC_delta", 0.0, "Onsager mobility of comp CD in delta phase");
  params.addParam<MaterialPropertyName>("dL_CD_muC_epsilon", 0.0, "TOnsager mobility of comp CD in epsilon phase");  
  //dL_DD_muC
  params.addParam<MaterialPropertyName>("dL_DD_muC_gamma", 0.0, "Onsager mobility of comp D in gamma phase");
  params.addParam<MaterialPropertyName>("dL_DD_muC_delta", 0.0, "Onsager mobility of comp D in delta phase");
  params.addParam<MaterialPropertyName>("dL_DD_muC_epsilon", 0.0, "Onsager mobility of comp D in epsilon phase"); 
  //*************************Derivatives of L w.r.t muD********************************************************//  
   //dL_BB_muD
  params.addParam<MaterialPropertyName>("dL_BB_muD_gamma", 0.0, "Onsager mobility of comp B in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BB_muD_delta", 0.0, "Onsager mobility of comp B in delta phase");
  params.addParam<MaterialPropertyName>("dL_BB_muD_epsilon", 0.0, "TOnsager mobility of comp B in epsilon phase");  
  //dL_BC_muD
  params.addParam<MaterialPropertyName>("dL_BC_muD_gamma", 0.0, "Onsager mobility of comp BC in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BC_muD_delta", 0.0, "Onsager mobility of comp BC in delta phase");
  params.addParam<MaterialPropertyName>("dL_BC_muD_epsilon", 0.0, "Onsager mobility of comp BC in epsilon phase");  
  //dL_BD_muD
  params.addParam<MaterialPropertyName>("dL_BD_muD_gamma", 0.0, "Onsager mobility of comp BD in gamma phase");
  params.addParam<MaterialPropertyName>("dL_BD_muD_delta", 0.0, "Onsager mobility of comp BD in delta phase");
  params.addParam<MaterialPropertyName>("dL_BD_muD_epsilon", 0.0, "TOnsager mobility of comp BD in epsilon phase");  
  //dL_CC_muD
  params.addParam<MaterialPropertyName>("dL_CC_muD_gamma", 0.0, "Onsager mobility of comp C in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CC_muD_delta", 0.0, "Onsager mobility of comp C in delta phase");
  params.addParam<MaterialPropertyName>("dL_CC_muD_epsilon", 0.0, "Onsager mobility of comp C in epsilon phase");  
  //dL_CD_muD
  params.addParam<MaterialPropertyName>("dL_CD_muD_gamma", 0.0, "Onsager mobility of comp CD in gamma phase");
  params.addParam<MaterialPropertyName>("dL_CD_muD_delta", 0.0, "Onsager mobility of comp CD in delta phase");
  params.addParam<MaterialPropertyName>("dL_CD_muD_epsilon", 0.0, "TOnsager mobility of comp CD in epsilon phase");  
  //dL_DD_muD
  params.addParam<MaterialPropertyName>("dL_DD_muD_gamma", 0.0, "Onsager mobility of comp D in gamma phase");
  params.addParam<MaterialPropertyName>("dL_DD_muD_delta", 0.0, "Onsager mobility of comp D in delta phase");
  params.addParam<MaterialPropertyName>("dL_DD_muD_epsilon", 0.0, "Onsager mobility of comp D in epsilon phase");  
  params.addRequiredParam<MaterialPropertyName>("h_alpha", "interpolation");
  return params; 
}

MultiCompMultiPhaseBase::MultiCompMultiPhaseBase(const InputParameters & parameters)
  : MultiPhaseBase(parameters),
  //inverse of thermodynamic factor
  _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
  _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
  _inv_B_tf_gamma(getMaterialProperty<Real>("inv_B_tf_gamma")),
  _inv_B_tf_delta(getMaterialProperty<Real>("inv_B_tf_delta")),
  _inv_B_tf_epsilon(getMaterialProperty<Real>("inv_B_tf_epsilon")),
  //inverse of thermodynamic factor
  _inv_BC_tf_alpha(getMaterialProperty<Real>("inv_BC_tf_alpha")),
  _inv_BC_tf_beta(getMaterialProperty<Real>("inv_BC_tf_beta")),
  _inv_BC_tf_gamma(getMaterialProperty<Real>("inv_BC_tf_gamma")),
  _inv_BC_tf_delta(getMaterialProperty<Real>("inv_BC_tf_delta")),
  _inv_BC_tf_epsilon(getMaterialProperty<Real>("inv_BC_tf_epsilon")),
  //inverse of thermodynamic factor
  _inv_BD_tf_alpha(getMaterialProperty<Real>("inv_BD_tf_alpha")),
  _inv_BD_tf_beta(getMaterialProperty<Real>("inv_BD_tf_beta")),
  _inv_BD_tf_gamma(getMaterialProperty<Real>("inv_BD_tf_gamma")),
  _inv_BD_tf_delta(getMaterialProperty<Real>("inv_BD_tf_delta")),
  _inv_BD_tf_epsilon(getMaterialProperty<Real>("inv_BD_tf_epsilon")),
   //inverse of thermodynamic factor
  _inv_C_tf_alpha(getMaterialProperty<Real>("inv_C_tf_alpha")),
  _inv_C_tf_beta(getMaterialProperty<Real>("inv_C_tf_beta")),
  _inv_C_tf_gamma(getMaterialProperty<Real>("inv_C_tf_gamma")),
  _inv_C_tf_delta(getMaterialProperty<Real>("inv_C_tf_delta")),
  _inv_C_tf_epsilon(getMaterialProperty<Real>("inv_C_tf_epsilon")),
   //inverse of thermodynamic factor
  _inv_CD_tf_alpha(getMaterialProperty<Real>("inv_CD_tf_alpha")),
  _inv_CD_tf_beta(getMaterialProperty<Real>("inv_CD_tf_beta")),
  _inv_CD_tf_gamma(getMaterialProperty<Real>("inv_CD_tf_gamma")),
  _inv_CD_tf_delta(getMaterialProperty<Real>("inv_CD_tf_delta")),
  _inv_CD_tf_epsilon(getMaterialProperty<Real>("inv_CD_tf_epsilon")),
   //inverse of thermodynamic factor
  _inv_D_tf_alpha(getMaterialProperty<Real>("inv_D_tf_alpha")),
  _inv_D_tf_beta(getMaterialProperty<Real>("inv_D_tf_beta")),
  _inv_D_tf_gamma(getMaterialProperty<Real>("inv_D_tf_gamma")),
  _inv_D_tf_delta(getMaterialProperty<Real>("inv_D_tf_delta")),
  _inv_D_tf_epsilon(getMaterialProperty<Real>("inv_D_tf_epsilon")),
  //Material properties to be interpolated
  //1) Onsager Mobility, 2) the derivative of OM 3) third derivative
  //L_BB within each phase
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  //assumed to be zero
  _L_BB_gamma(getMaterialProperty<Real>("L_BB_gamma")),
  _L_BB_delta(getMaterialProperty<Real>("L_BB_delta")),
  _L_BB_epsilon(getMaterialProperty<Real>("L_BB_epsilon")),
    //L_BC within each phase
  _L_BC_alpha(getMaterialProperty<Real>("L_BC_alpha")),
  _L_BC_beta(getMaterialProperty<Real>("L_BC_beta")),
  _L_BC_gamma(getMaterialProperty<Real>("L_BC_gamma")),
  _L_BC_delta(getMaterialProperty<Real>("L_BC_delta")),
  _L_BC_epsilon(getMaterialProperty<Real>("L_BC_epsilon")),
   //L_BD within each phase
  _L_BD_alpha(getMaterialProperty<Real>("L_BD_alpha")),
  _L_BD_beta(getMaterialProperty<Real>("L_BD_beta")),
  _L_BD_gamma(getMaterialProperty<Real>("L_BD_gamma")),
  _L_BD_delta(getMaterialProperty<Real>("L_BD_delta")),
  _L_BD_epsilon(getMaterialProperty<Real>("L_BD_epsilon")),
  //L_CC within each phase
  _L_CC_alpha(getMaterialProperty<Real>("L_CC_alpha")),
  _L_CC_beta(getMaterialProperty<Real>("L_CC_beta")),
  _L_CC_gamma(getMaterialProperty<Real>("L_CC_gamma")),
  _L_CC_delta(getMaterialProperty<Real>("L_CC_delta")),
  _L_CC_epsilon(getMaterialProperty<Real>("L_CC_epsilon")),
   //L_CD within each phase
  _L_CD_alpha(getMaterialProperty<Real>("L_CD_alpha")),
  _L_CD_beta(getMaterialProperty<Real>("L_CD_beta")),
  _L_CD_gamma(getMaterialProperty<Real>("L_CD_gamma")),
  _L_CD_delta(getMaterialProperty<Real>("L_CD_delta")),
  _L_CD_epsilon(getMaterialProperty<Real>("L_CD_epsilon")),
   //L_D within each phase
  _L_DD_alpha(getMaterialProperty<Real>("L_DD_alpha")),
  _L_DD_beta(getMaterialProperty<Real>("L_DD_beta")),
  _L_DD_gamma(getMaterialProperty<Real>("L_DD_gamma")),
  _L_DD_delta(getMaterialProperty<Real>("L_DD_delta")),
  _L_DD_epsilon(getMaterialProperty<Real>("L_DD_epsilon")),
  //dL_BB_muB
  _dL_BB_muB_alpha(getMaterialProperty<Real>("dL_BB_muB_alpha")),
  _dL_BB_muB_beta(getMaterialProperty<Real>("dL_BB_muB_beta")),
  _dL_BB_muB_gamma(getMaterialProperty<Real>("dL_BB_muB_gamma")),
  _dL_BB_muB_delta(getMaterialProperty<Real>("dL_BB_muB_delta")),
  _dL_BB_muB_epsilon(getMaterialProperty<Real>("dL_BB_muB_epsilon")),
  //dL_BC_muB
  _dL_BC_muB_alpha(getMaterialProperty<Real>("dL_BC_muB_alpha")),
  _dL_BC_muB_beta(getMaterialProperty<Real>("dL_BC_muB_beta")),
  _dL_BC_muB_gamma(getMaterialProperty<Real>("dL_BC_muB_gamma")),
  _dL_BC_muB_delta(getMaterialProperty<Real>("dL_BC_muB_delta")),
  _dL_BC_muB_epsilon(getMaterialProperty<Real>("dL_BC_muB_epsilon")),
  //dL_BD_muB
  _dL_BD_muB_alpha(getMaterialProperty<Real>("dL_BD_muB_alpha")),
  _dL_BD_muB_beta(getMaterialProperty<Real>("dL_BD_muB_beta")),
  _dL_BD_muB_gamma(getMaterialProperty<Real>("dL_BD_muB_gamma")),
  _dL_BD_muB_delta(getMaterialProperty<Real>("dL_BD_muB_delta")),
  _dL_BD_muB_epsilon(getMaterialProperty<Real>("dL_BD_muB_epsilon")),
  //dL_CC_muB
  _dL_CC_muB_alpha(getMaterialProperty<Real>("dL_CC_muB_alpha")),
  _dL_CC_muB_beta(getMaterialProperty<Real>("dL_CC_muB_beta")),
  _dL_CC_muB_gamma(getMaterialProperty<Real>("dL_CC_muB_gamma")),
  _dL_CC_muB_delta(getMaterialProperty<Real>("dL_CC_muB_delta")),
  _dL_CC_muB_epsilon(getMaterialProperty<Real>("dL_CC_muB_epsilon")),
  //dL_CD_muB
  _dL_CD_muB_alpha(getMaterialProperty<Real>("dL_CD_muB_alpha")),
  _dL_CD_muB_beta(getMaterialProperty<Real>("dL_CD_muB_beta")),
  _dL_CD_muB_gamma(getMaterialProperty<Real>("dL_CD_muB_gamma")),
  _dL_CD_muB_delta(getMaterialProperty<Real>("dL_CD_muB_delta")),
  _dL_CD_muB_epsilon(getMaterialProperty<Real>("dL_CD_muB_epsilon")),
  //dL_DD_muB
  _dL_DD_muB_alpha(getMaterialProperty<Real>("dL_DD_muB_alpha")),
  _dL_DD_muB_beta(getMaterialProperty<Real>("dL_DD_muB_beta")),
  _dL_DD_muB_gamma(getMaterialProperty<Real>("dL_DD_muB_gamma")),
  _dL_DD_muB_delta(getMaterialProperty<Real>("dL_DD_muB_delta")),
  _dL_DD_muB_epsilon(getMaterialProperty<Real>("dL_DD_muB_epsilon")),
   //dL_BB_muC
  _dL_BB_muC_alpha(getMaterialProperty<Real>("dL_BB_muC_alpha")),
  _dL_BB_muC_beta(getMaterialProperty<Real>("dL_BB_muC_beta")),
  _dL_BB_muC_gamma(getMaterialProperty<Real>("dL_BB_muC_gamma")),
  _dL_BB_muC_delta(getMaterialProperty<Real>("dL_BB_muC_delta")),
  _dL_BB_muC_epsilon(getMaterialProperty<Real>("dL_BB_muC_epsilon")),
  //dL_BC_muC
  _dL_BC_muC_alpha(getMaterialProperty<Real>("dL_BC_muC_alpha")),
  _dL_BC_muC_beta(getMaterialProperty<Real>("dL_BC_muC_beta")),
  _dL_BC_muC_gamma(getMaterialProperty<Real>("dL_BC_muC_gamma")),
  _dL_BC_muC_delta(getMaterialProperty<Real>("dL_BC_muC_delta")),
  _dL_BC_muC_epsilon(getMaterialProperty<Real>("dL_BC_muC_epsilon")),
  //dL_BD_muC
  _dL_BD_muC_alpha(getMaterialProperty<Real>("dL_BD_muC_alpha")),
  _dL_BD_muC_beta(getMaterialProperty<Real>("dL_BD_muC_beta")),
  _dL_BD_muC_gamma(getMaterialProperty<Real>("dL_BD_muC_gamma")),
  _dL_BD_muC_delta(getMaterialProperty<Real>("dL_BD_muC_delta")),
  _dL_BD_muC_epsilon(getMaterialProperty<Real>("dL_BD_muC_epsilon")),
  //dL_CC_muC
  _dL_CC_muC_alpha(getMaterialProperty<Real>("dL_CC_muC_alpha")),
  _dL_CC_muC_beta(getMaterialProperty<Real>("dL_CC_muC_beta")),
  _dL_CC_muC_gamma(getMaterialProperty<Real>("dL_CC_muC_gamma")),
  _dL_CC_muC_delta(getMaterialProperty<Real>("dL_CC_muC_delta")),
  _dL_CC_muC_epsilon(getMaterialProperty<Real>("dL_CC_muC_epsilon")),
  //dL_CD_muC
  _dL_CD_muC_alpha(getMaterialProperty<Real>("dL_CD_muC_alpha")),
  _dL_CD_muC_beta(getMaterialProperty<Real>("dL_CD_muC_beta")),
  _dL_CD_muC_gamma(getMaterialProperty<Real>("dL_CD_muC_gamma")),
  _dL_CD_muC_delta(getMaterialProperty<Real>("dL_CD_muC_delta")),
  _dL_CD_muC_epsilon(getMaterialProperty<Real>("dL_CD_muC_epsilon")),
  //dL_DD_muC
  _dL_DD_muC_alpha(getMaterialProperty<Real>("dL_DD_muC_alpha")),
  _dL_DD_muC_beta(getMaterialProperty<Real>("dL_DD_muC_beta")),
  _dL_DD_muC_gamma(getMaterialProperty<Real>("dL_DD_muC_gamma")),
  _dL_DD_muC_delta(getMaterialProperty<Real>("dL_DD_muC_delta")),
  _dL_DD_muC_epsilon(getMaterialProperty<Real>("dL_DD_muC_epsilon")),
   //dL_BB_muD
  _dL_BB_muD_alpha(getMaterialProperty<Real>("dL_BB_muD_alpha")),
  _dL_BB_muD_beta(getMaterialProperty<Real>("dL_BB_muD_beta")),
  _dL_BB_muD_gamma(getMaterialProperty<Real>("dL_BB_muD_gamma")),
  _dL_BB_muD_delta(getMaterialProperty<Real>("dL_BB_muD_delta")),
  _dL_BB_muD_epsilon(getMaterialProperty<Real>("dL_BB_muD_epsilon")),
  //dL_BC_muC
  _dL_BC_muD_alpha(getMaterialProperty<Real>("dL_BC_muD_alpha")),
  _dL_BC_muD_beta(getMaterialProperty<Real>("dL_BC_muD_beta")),
  _dL_BC_muD_gamma(getMaterialProperty<Real>("dL_BC_muD_gamma")),
  _dL_BC_muD_delta(getMaterialProperty<Real>("dL_BC_muD_delta")),
  _dL_BC_muD_epsilon(getMaterialProperty<Real>("dL_BC_muD_epsilon")),
  //dL_BD_muD
  _dL_BD_muD_alpha(getMaterialProperty<Real>("dL_BD_muD_alpha")),
  _dL_BD_muD_beta(getMaterialProperty<Real>("dL_BD_muD_beta")),
  _dL_BD_muD_gamma(getMaterialProperty<Real>("dL_BD_muD_gamma")),
  _dL_BD_muD_delta(getMaterialProperty<Real>("dL_BD_muD_delta")),
  _dL_BD_muD_epsilon(getMaterialProperty<Real>("dL_BD_muD_epsilon")),
  //dL_CC_muD
  _dL_CC_muD_alpha(getMaterialProperty<Real>("dL_CC_muD_alpha")),
  _dL_CC_muD_beta(getMaterialProperty<Real>("dL_CC_muD_beta")),
  _dL_CC_muD_gamma(getMaterialProperty<Real>("dL_CC_muD_gamma")),
  _dL_CC_muD_delta(getMaterialProperty<Real>("dL_CC_muD_delta")),
  _dL_CC_muD_epsilon(getMaterialProperty<Real>("dL_CC_muD_epsilon")),
  //dL_CD_muD
  _dL_CD_muD_alpha(getMaterialProperty<Real>("dL_CD_muD_alpha")),
  _dL_CD_muD_beta(getMaterialProperty<Real>("dL_CD_muD_beta")),
  _dL_CD_muD_gamma(getMaterialProperty<Real>("dL_CD_muD_gamma")),
  _dL_CD_muD_delta(getMaterialProperty<Real>("dL_CD_muD_delta")),
  _dL_CD_muD_epsilon(getMaterialProperty<Real>("dL_CD_muD_epsilon")),
  //dL_DD_muD
  _dL_DD_muD_alpha(getMaterialProperty<Real>("dL_DD_muD_alpha")),
  _dL_DD_muD_beta(getMaterialProperty<Real>("dL_DD_muD_beta")),
  _dL_DD_muD_gamma(getMaterialProperty<Real>("dL_DD_muD_gamma")),
  _dL_DD_muD_delta(getMaterialProperty<Real>("dL_DD_muD_delta")),
  _dL_DD_muD_epsilon(getMaterialProperty<Real>("dL_DD_muD_epsilon"))
{
}

//****************************************************************************//
//The following steps is followed here:
//1) Create the coefficents of the overall susceptibility matrix
//2) Generate the determinant of the matrix
//3) Create the coefficients of the overall thermodynamic factor matrix
//4) Interpolate the coeff. of the Onsager mobilities
//5) Perform a matrix product to obtain the Overall diffusioon coeff. 
//****************************************************************************//

//****************************************************************************//
//Step 1) Susceptibility matrix coefficients
//      overall_chi = [chi_BB chi_BC chi_BD
//                     chi_BC chi_CC chi_CD
//                     chi_BD chi_CD chi_DD]
//****************************************************************************//

Real
MultiCompMultiPhaseBase::chi_BB() const
{
   return (_h_alpha[_qp]  * _inv_B_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_B_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_B_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_B_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_B_tf_epsilon[_qp]);
 
}

Real
MultiCompMultiPhaseBase::chi_BC() const
{
   return (_h_alpha[_qp]  * _inv_BC_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_BC_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_BC_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_BC_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_BC_tf_epsilon[_qp]);
 
}

Real
MultiCompMultiPhaseBase::chi_BD() const
{
   return (_h_alpha[_qp]  * _inv_BD_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_BD_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_BD_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_BD_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_BD_tf_epsilon[_qp]);
 
}

Real
MultiCompMultiPhaseBase::chi_CC() const
{
   return (_h_alpha[_qp]  * _inv_C_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_C_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_C_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_C_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_C_tf_epsilon[_qp]);
 
}

Real
MultiCompMultiPhaseBase::chi_CD() const
{
   return (_h_alpha[_qp]  * _inv_CD_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_CD_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_CD_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_CD_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_CD_tf_epsilon[_qp]);
 
}

Real
MultiCompMultiPhaseBase::chi_DD() const
{
   return (_h_alpha[_qp]  * _inv_D_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_D_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_D_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_D_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_D_tf_epsilon[_qp]);
 
}
//****************************************************************************//
//Step 2) Generate the determinant of the overall susceptibility matrix
//      determinanr_chi = [chi_BB chi_BC chi_BD
//                     chi_BC chi_CC chi_CD
//                     chi_BD chi_CD chi_DD]
//****************************************************************************//

Real
MultiCompMultiPhaseBase::det_chi() const
{
  return (MultiCompMultiPhaseBase::chi_BB()* 
            (MultiCompMultiPhaseBase::chi_CC()*MultiCompMultiPhaseBase::chi_DD()
            -MultiCompMultiPhaseBase::chi_CD()*MultiCompMultiPhaseBase::chi_CD())
         -MultiCompMultiPhaseBase::chi_BC()*
            (MultiCompMultiPhaseBase::chi_BC()*MultiCompMultiPhaseBase::chi_DD()
            -MultiCompMultiPhaseBase::chi_CD()*MultiCompMultiPhaseBase::chi_BD())
         +MultiCompMultiPhaseBase::chi_BD()*
            (MultiCompMultiPhaseBase::chi_BC()*MultiCompMultiPhaseBase::chi_CD()
            -MultiCompMultiPhaseBase::chi_BD()*MultiCompMultiPhaseBase::chi_CC()));  
}

//****************************************************************************//
//Step 3) Coefficents of the overall thermodynamic factor matrix
//      theta_overall = [theta_BB theta_BC theta_BD
//                      theta_BC  theta_CC theta_CD
//                      theta_BD  theta_CD theta_DD]
// All coeffecients of the overall thermodynamic factor matrix
// In this case, we will first construct the diagonal elements
// and then the off-diagonal components
//****************************************************************************//

Real
MultiCompMultiPhaseBase::thermodynamic_factorBB() const
{
   //The coeff BB of the overall TF matrix is:
   //1/det(overall_chi)*(chi_CC*chi_DD - chi_CD*chi_CD);
   return ( (1.0/MultiCompMultiPhaseBase::det_chi()) 
              * (MultiCompMultiPhaseBase::chi_CC()*MultiCompMultiPhaseBase::chi_DD() 
                -MultiCompMultiPhaseBase::chi_CD()*MultiCompMultiPhaseBase::chi_CD())); 
}

Real
MultiCompMultiPhaseBase::thermodynamic_factorCC() const
{                     
   //The coeff CC of the overall TF matrix is:
   //1/det(overall_chi)*(chi_BB*chi_DD - chi_BD*chi_BD);
   return ( (1.0/MultiCompMultiPhaseBase::det_chi()) 
              * (MultiCompMultiPhaseBase::chi_BB()*MultiCompMultiPhaseBase::chi_DD() 
                -MultiCompMultiPhaseBase::chi_BD()*MultiCompMultiPhaseBase::chi_BD())); 
}

Real
MultiCompMultiPhaseBase::thermodynamic_factorDD() const
{   
     //The coeff CC of the overall TF matrix is:
    //1/det(overall_chi)*(chi_BB*chi_CC - chi_BC*chi_BC);              
    return ( (1.0/MultiCompMultiPhaseBase::det_chi()) 
              * (MultiCompMultiPhaseBase::chi_BB()*MultiCompMultiPhaseBase::chi_CC() 
                -MultiCompMultiPhaseBase::chi_BC()*MultiCompMultiPhaseBase::chi_BC()));                
}

//Now, we construct the off-diagonal components

Real
MultiCompMultiPhaseBase::thermodynamic_factorBC() const
{
   //The coeff BC of the overall TF matrix is:
   //-1/det(overall_chi)*(chi_BC*chi_BD - chi_CD*chi_BD);
   return -( (1.0/MultiCompMultiPhaseBase::det_chi()) 
              * (MultiCompMultiPhaseBase::chi_BC()*MultiCompMultiPhaseBase::chi_BD()
                -MultiCompMultiPhaseBase::chi_CD()*MultiCompMultiPhaseBase::chi_BD()));
}

Real
MultiCompMultiPhaseBase::thermodynamic_factorBD() const
{
    //The coeff BD of the overall TF matrix is:
    //1/det(overall_chi)*(chi_BC*chi_CD - chi_CC*chi_BD);
    return ( (1.0/MultiCompMultiPhaseBase::det_chi()) 
              * (MultiCompMultiPhaseBase::chi_BC()*MultiCompMultiPhaseBase::chi_CD()
                -MultiCompMultiPhaseBase::chi_CC()*MultiCompMultiPhaseBase::chi_BD()));
}

Real
MultiCompMultiPhaseBase::thermodynamic_factorCD() const
{
    //The coeff BD of the overall TF matrix is:
    //1/det(overall_chi)*(chi_BB*chi_CD - chi_BC*chi_BD);
    return ( (1.0/MultiCompMultiPhaseBase::det_chi()) 
              * (MultiCompMultiPhaseBase::chi_BB()*MultiCompMultiPhaseBase::chi_CD()
                -MultiCompMultiPhaseBase::chi_BC()*MultiCompMultiPhaseBase::chi_BD()));
                      
}

//****************************************************************************//
//Step 4) Interpolate the coefficients of the overall Onsager mobility matrix
//      L_overall = [L_BB L_BC L_BD
//                   L_BC L_CC L_CD
//                   L_BD L_CD L_DD]
// All coeffecients of the overall thermodynamic factor matrix
//****************************************************************************//


Real
MultiCompMultiPhaseBase::L_BB_interp() const
{
  return (_L_BB_alpha[_qp] *_h_alpha[_qp] 
        + _L_BB_beta[_qp]  * _h_beta[_qp] 
        + _L_BB_gamma[_qp] * _h_gamma[_qp] 
        + _L_BB_delta[_qp] * _h_delta[_qp]
        + _L_BB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::L_BC_interp() const
{
  return (_L_BC_alpha[_qp] * _h_alpha[_qp] 
        + _L_BC_beta[_qp]  * _h_beta[_qp] 
        + _L_BC_gamma[_qp] * _h_gamma[_qp] 
        + _L_BC_delta[_qp] * _h_delta[_qp]
        + _L_BC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::L_BD_interp() const
{
  return (_L_BD_alpha[_qp]  * _h_alpha[_qp] 
        + _L_BD_beta[_qp]   * _h_beta[_qp] 
        + _L_BD_gamma[_qp]  * _h_gamma[_qp] 
        + _L_BD_delta[_qp]  * _h_delta[_qp]
        + _L_BD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::L_CC_interp() const
{
  return (_L_CC_alpha[_qp]  * _h_alpha[_qp] 
        + _L_CC_beta[_qp]   * _h_beta[_qp] 
        + _L_CC_gamma[_qp]  * _h_gamma[_qp] 
        + _L_CC_delta[_qp]  * _h_delta[_qp]
        + _L_CC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::L_CD_interp() const
{
  return (_L_CD_alpha[_qp]  * _h_alpha[_qp] 
        + _L_CD_beta[_qp]   * _h_beta[_qp] 
        + _L_CD_gamma[_qp]  * _h_gamma[_qp] 
        + _L_CD_delta[_qp]  * _h_delta[_qp]
        + _L_CD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::L_DD_interp() const
{
  return (_L_DD_alpha[_qp]  * _h_alpha[_qp] 
        + _L_DD_beta[_qp]   * _h_beta[_qp] 
        + _L_DD_gamma[_qp]  * _h_gamma[_qp] 
        + _L_DD_delta[_qp]  * _h_delta[_qp]
        + _L_DD_epsilon[_qp]* _h_epsilon[_qp]);
}

//From the Onsager mobility and the thermodynamic factor, we can derive the 
//interdiffusion coefficents, note that these need not be symmetric
//First is for the flux of component of B

//****************************************************************************//
//Step 5) Matrix product of the overall thermodynamic factor matrix with 
//      DC_overall =   [DC_BB DC_BC DC_BD
//                      DC_CB DC_CC DC_CD
//                      DC_DB DC_DC DC_DD]
// All coeffecients of the overall thermodynamic factor matrix
//****************************************************************************//

Real
MultiCompMultiPhaseBase::DC_BB_interp() const
{ 
  return (MultiCompMultiPhaseBase::L_BB_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBB()
         +MultiCompMultiPhaseBase::L_BC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBC()
         +MultiCompMultiPhaseBase::L_BD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBD());
}

Real
MultiCompMultiPhaseBase::DC_BC_interp() const
{
  //Note thermodynamic factor CD = thermodynamic factor DC
  return (MultiCompMultiPhaseBase::L_BB_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBC()
         +MultiCompMultiPhaseBase::L_BC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCC()
         +MultiCompMultiPhaseBase::L_BD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCD());
}

Real
MultiCompMultiPhaseBase::DC_BD_interp() const
{  
  return (MultiCompMultiPhaseBase::L_BB_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBD()
         +MultiCompMultiPhaseBase::L_BC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCD()
         +MultiCompMultiPhaseBase::L_BD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorDD());
}
//First is for the flux of component of C
Real
MultiCompMultiPhaseBase::DC_CB_interp() const
{
  //Note Onsager Mobility BC = Onsager Mobility BC
  return (MultiCompMultiPhaseBase::L_BC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBB()
         +MultiCompMultiPhaseBase::L_CC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBC()
         +MultiCompMultiPhaseBase::L_CD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBD());
}

Real
MultiCompMultiPhaseBase::DC_CC_interp() const
{
  //Note Onsager Mobility BC = Onsager Mobility BC
  return (MultiCompMultiPhaseBase::L_BC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBC()
         +MultiCompMultiPhaseBase::L_CC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCC()
         +MultiCompMultiPhaseBase::L_CD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCD());
}

Real
MultiCompMultiPhaseBase::DC_CD_interp() const
{
  //Note Onsager Mobility CB = Onsager Mobility BC
  return (MultiCompMultiPhaseBase::L_BC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBD()
         +MultiCompMultiPhaseBase::L_CC_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCD()
         +MultiCompMultiPhaseBase::L_CD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorDD());
}

//First is for the flux of component of D
//Note that DC_BD is not equal to DC_DB
Real
MultiCompMultiPhaseBase::DC_DB_interp() const
{
  //Note Onsager Mobility DB = Onsager Mobility BD
  // Thermodynamic factor DB = Thermodynamic factor BD 
  return (MultiCompMultiPhaseBase::L_BD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBB()
         +MultiCompMultiPhaseBase::L_CD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBC()
         +MultiCompMultiPhaseBase::L_DD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBD());
}

//Note that DC_CD is not equal to DC_DC
Real
MultiCompMultiPhaseBase::DC_DC_interp() const
{
  //Note Onsager Mobility DB = Onsager Mobility BD
  return (MultiCompMultiPhaseBase::L_BD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBC()
         +MultiCompMultiPhaseBase::L_CD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCC()
         +MultiCompMultiPhaseBase::L_DD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCD());
}

Real
MultiCompMultiPhaseBase::DC_DD_interp() const
{
  //Note Onsager Mobility BC = Onsager Mobility BC
  return (MultiCompMultiPhaseBase::L_BD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorBD()
         +MultiCompMultiPhaseBase::L_CD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorCD()
         +MultiCompMultiPhaseBase::L_DD_interp()*MultiCompMultiPhaseBase::thermodynamic_factorDD());
}

// To model the composition dependence
//of the Onsager mobilities this is required
Real
MultiCompMultiPhaseBase::dL_BB_muB_interp() const
{
         
  return ( _dL_BB_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BB_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_BB_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BB_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_BB_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_BC_muB_interp() const
{
  return ( _dL_BC_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BC_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_BC_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BC_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_BC_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_BD_muB_interp() const
{
  return ( _dL_BD_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BD_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_BD_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BD_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_BD_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_CC_muB_interp() const
{
  return ( _dL_CC_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CC_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_CC_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CC_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_CC_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_CD_muB_interp() const
{
  return ( _dL_CD_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CD_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_CD_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CD_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_CD_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_DD_muB_interp() const
{
  return ( _dL_DD_muB_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_DD_muB_beta[_qp]   * _h_beta[_qp] 
         + _dL_DD_muB_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_DD_muB_delta[_qp]  * _h_delta[_qp]
         + _dL_DD_muB_epsilon[_qp]* _h_epsilon[_qp]);
}


/// w.r.t component C
Real
MultiCompMultiPhaseBase::dL_BB_muC_interp() const
{
  return ( _dL_BB_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BB_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_BB_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BB_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_BB_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_BC_muC_interp() const
{
  return ( _dL_BC_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BC_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_BC_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BC_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_BC_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_BD_muC_interp() const
{
  return ( _dL_BD_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BD_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_BD_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BD_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_BD_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_CC_muC_interp() const
{
  return ( _dL_CC_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CC_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_CC_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CC_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_CC_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_CD_muC_interp() const
{
  return ( _dL_CD_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CD_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_CD_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CD_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_CD_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_DD_muC_interp() const
{
  return ( _dL_DD_muC_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_DD_muC_beta[_qp]   * _h_beta[_qp] 
         + _dL_DD_muC_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_DD_muC_delta[_qp]  * _h_delta[_qp]
         + _dL_DD_muC_epsilon[_qp]* _h_epsilon[_qp]);
}

/// w.r.t component D

Real
MultiCompMultiPhaseBase::dL_BB_muD_interp() const
{
  return ( _dL_BB_muD_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BB_muD_beta[_qp]   * _h_beta[_qp] 
         + _dL_BB_muD_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BB_muD_delta[_qp]  * _h_delta[_qp]
         + _dL_BB_muD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_BC_muD_interp() const
{
  return ( _dL_BC_muD_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BC_muD_beta[_qp]   * _h_beta[_qp] 
         + _dL_BC_muD_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BC_muD_delta[_qp]  * _h_delta[_qp]
         + _dL_BC_muD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_BD_muD_interp() const
{
  return ( _dL_BD_muD_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_BD_muD_beta[_qp]   * _h_beta[_qp] 
         + _dL_BD_muD_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_BD_muD_delta[_qp]  * _h_delta[_qp]
         + _dL_BD_muD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_CC_muD_interp() const
{
  return ( _dL_CC_muD_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CC_muD_beta[_qp]   * _h_beta[_qp] 
         + _dL_CC_muD_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CC_muD_delta[_qp]  * _h_delta[_qp]
         + _dL_CC_muD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_CD_muD_interp() const
{
  return ( _dL_CD_muD_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_CD_muD_beta[_qp]   * _h_beta[_qp] 
         + _dL_CD_muD_gamma[_qp]  * _h_gamma[_qp] 
         + _dL_CD_muD_delta[_qp]  * _h_delta[_qp]
         + _dL_CD_muD_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
MultiCompMultiPhaseBase::dL_DD_muD_interp() const
{
  return (_dL_DD_muD_alpha[_qp]  * _h_alpha[_qp] 
         + _dL_DD_muD_beta[_qp]  * _h_beta[_qp] 
         + _dL_DD_muD_gamma[_qp] * _h_gamma[_qp] 
         + _dL_DD_muD_delta[_qp] * _h_delta[_qp]
         + _dL_DD_muD_epsilon[_qp]* _h_epsilon[_qp]);
}

//First derivatives  of each material property sum over all interpolation function
//This will be (# of property*# of phases)
//Onsger mobility L_BB

Real
MultiCompMultiPhaseBase::sum_dh_L_BB_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_BB_beta[_qp]  - _L_BB_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_BB_gamma[_qp] - _L_BB_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_BB_delta[_qp] - _L_BB_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_BB_epsilon[_qp]- _L_BB_alpha[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BB_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_BB_alpha[_qp] - _L_BB_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_BB_gamma[_qp] - _L_BB_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_BB_delta[_qp] - _L_BB_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_BB_epsilon[_qp]- _L_BB_beta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BB_gamma() const
{
    return (_dhalpha_dphigamma[_qp]     * (_L_BB_alpha[_qp] - _L_BB_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_BB_beta[_qp]  - _L_BB_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_BB_delta[_qp] - _L_BB_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_BB_epsilon[_qp]-_L_BB_gamma[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BB_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_BB_alpha[_qp] - _L_BB_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_BB_beta[_qp]  - _L_BB_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_BB_gamma[_qp] - _L_BB_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_BB_epsilon[_qp]-_L_BB_delta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BB_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_BB_alpha[_qp] - _L_BB_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_BB_beta[_qp]  - _L_BB_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_BB_gamma[_qp] - _L_BB_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]*   (_L_BB_delta[_qp] - _L_BB_epsilon[_qp]));
}

//Onsger mobility L_BC

Real
MultiCompMultiPhaseBase::sum_dh_L_BC_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_BC_beta[_qp]  - _L_BC_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_BC_gamma[_qp] - _L_BC_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_BC_delta[_qp] - _L_BC_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_BC_epsilon[_qp]- _L_BC_alpha[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BC_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_BC_alpha[_qp] - _L_BC_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_BC_gamma[_qp] - _L_BC_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_BC_delta[_qp] - _L_BC_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_BC_epsilon[_qp]- _L_BC_beta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BC_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_BC_alpha[_qp] - _L_BC_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_BC_beta[_qp]  - _L_BC_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_BC_delta[_qp] - _L_BC_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_BC_epsilon[_qp]-_L_BC_gamma[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BC_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_BC_alpha[_qp] - _L_BC_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_BC_beta[_qp]  - _L_BC_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_BC_gamma[_qp] - _L_BC_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_BC_epsilon[_qp]-_L_BC_delta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BC_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_BC_alpha[_qp] - _L_BC_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_BC_beta[_qp]  - _L_BC_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_BC_gamma[_qp] - _L_BC_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]*   (_L_BC_delta[_qp] - _L_BC_epsilon[_qp]));
}

//Property::sum_dh_L_BD_phases()

Real
MultiCompMultiPhaseBase::sum_dh_L_BD_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_BD_beta[_qp]  - _L_BD_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_BD_gamma[_qp] - _L_BD_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_BD_delta[_qp] - _L_BD_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_BD_epsilon[_qp]- _L_BD_alpha[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BD_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_BD_alpha[_qp] - _L_BD_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_BD_gamma[_qp] - _L_BD_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_BD_delta[_qp] - _L_BD_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_BD_epsilon[_qp]- _L_BD_beta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BD_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_BD_alpha[_qp] - _L_BD_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_BD_beta[_qp]  - _L_BD_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_BD_delta[_qp] - _L_BD_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_BD_epsilon[_qp]-_L_BD_gamma[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BD_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_BD_alpha[_qp] - _L_BD_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_BD_beta[_qp]  - _L_BD_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_BD_gamma[_qp] - _L_BD_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_BD_epsilon[_qp]-_L_BD_delta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_BD_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_BD_alpha[_qp] - _L_BD_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_BD_beta[_qp]  - _L_BD_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_BD_gamma[_qp] - _L_BD_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]*   (_L_BD_delta[_qp] - _L_BD_epsilon[_qp]));
}


//Property::sum_dh_L_CC_phases()

Real
MultiCompMultiPhaseBase::sum_dh_L_CC_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_CC_beta[_qp]  - _L_CC_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_CC_gamma[_qp] - _L_CC_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_CC_delta[_qp] - _L_CC_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_CC_epsilon[_qp]- _L_CC_alpha[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CC_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_CC_alpha[_qp] - _L_CC_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_CC_gamma[_qp] - _L_CC_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_CC_delta[_qp] - _L_CC_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_CC_epsilon[_qp]- _L_CC_beta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CC_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_CC_alpha[_qp] - _L_CC_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_CC_beta[_qp]  - _L_CC_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_CC_delta[_qp] - _L_CC_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_CC_epsilon[_qp]-_L_CC_gamma[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CC_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_CC_alpha[_qp] - _L_CC_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_CC_beta[_qp]  - _L_CC_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_CC_gamma[_qp] - _L_CC_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_CC_epsilon[_qp]-_L_CC_delta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CC_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_CC_alpha[_qp] - _L_CC_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_CC_beta[_qp]  - _L_CC_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_CC_gamma[_qp] - _L_CC_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]*   (_L_CC_delta[_qp] - _L_CC_epsilon[_qp]));
}

//Property::sum_dh_L_CD_phases()

Real
MultiCompMultiPhaseBase::sum_dh_L_CD_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_CD_beta[_qp]  - _L_CD_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_CD_gamma[_qp] - _L_CD_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_CD_delta[_qp] - _L_CD_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_CD_epsilon[_qp]- _L_CD_alpha[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CD_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_CD_alpha[_qp] - _L_CD_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_CD_gamma[_qp] - _L_CD_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_CD_delta[_qp] - _L_CD_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_CD_epsilon[_qp]- _L_CD_beta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CD_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_CD_alpha[_qp] - _L_CD_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_CD_beta[_qp]  - _L_CD_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_CD_delta[_qp] - _L_CD_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_CD_epsilon[_qp]-_L_CD_gamma[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CD_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_CD_alpha[_qp] - _L_CD_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_CD_beta[_qp]  - _L_CD_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_CD_gamma[_qp] - _L_CD_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_CD_epsilon[_qp]-_L_CD_delta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_CD_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_CD_alpha[_qp] - _L_CD_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_CD_beta[_qp]  - _L_CD_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_CD_gamma[_qp] - _L_CD_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]  * (_L_CD_delta[_qp] - _L_CD_epsilon[_qp]));
}

//Property::sum_dh_L_DD_phases()

Real
MultiCompMultiPhaseBase::sum_dh_L_DD_alpha() const
{
    return (_dhbeta_dphialpha[_qp]      * (_L_DD_beta[_qp]  - _L_DD_alpha[_qp])
             + _dhgamma_dphialpha[_qp]  * (_L_DD_gamma[_qp] - _L_DD_alpha[_qp])
             + _dhdelta_dphialpha[_qp]  * (_L_DD_delta[_qp] - _L_DD_alpha[_qp])
             + _dhepsilon_dphialpha[_qp]* (_L_DD_epsilon[_qp]- _L_DD_alpha[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_DD_beta() const
{
    return (_dhalpha_dphibeta[_qp]     * (_L_DD_alpha[_qp] - _L_DD_beta[_qp])
             + _dhgamma_dphibeta[_qp]  * (_L_DD_gamma[_qp] - _L_DD_beta[_qp])
             + _dhdelta_dphibeta[_qp]  * (_L_DD_delta[_qp] - _L_DD_beta[_qp])
             + _dhepsilon_dphibeta[_qp]*(_L_DD_epsilon[_qp]- _L_DD_beta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_DD_gamma() const
{
    return  (_dhalpha_dphigamma[_qp]    * (_L_DD_alpha[_qp] - _L_DD_gamma[_qp])
             + _dhbeta_dphigamma[_qp]   * (_L_DD_beta[_qp]  - _L_DD_gamma[_qp])
             + _dhdelta_dphigamma[_qp]  * (_L_DD_delta[_qp] - _L_DD_gamma[_qp])
             + _dhepsilon_dphigamma[_qp]* (_L_DD_epsilon[_qp]-_L_DD_gamma[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_DD_delta() const
{
    return (_dhalpha_dphidelta[_qp]     * (_L_DD_alpha[_qp] - _L_DD_delta[_qp])
             + _dhbeta_dphidelta[_qp]   * (_L_DD_beta[_qp]  - _L_DD_delta[_qp])
             + _dhgamma_dphidelta[_qp]  * (_L_DD_gamma[_qp] - _L_DD_delta[_qp])
             + _dhepsilon_dphidelta[_qp]* (_L_DD_epsilon[_qp]-_L_DD_delta[_qp]));
}

Real
MultiCompMultiPhaseBase::sum_dh_L_DD_epsilon() const
{
    return (_dhalpha_dphiepsilon[_qp]     * (_L_DD_alpha[_qp] - _L_DD_epsilon[_qp])
             + _dhbeta_dphiepsilon[_qp]   * (_L_DD_beta[_qp]  - _L_DD_epsilon[_qp])
             + _dhgamma_dphiepsilon[_qp]  * (_L_DD_gamma[_qp] - _L_DD_epsilon[_qp])
             + _dhdelta_dphiepsilon[_qp]  * (_L_DD_delta[_qp] - _L_DD_epsilon[_qp]));
}

Real
MultiCompMultiPhaseBase::computeQpResidual(){
  return 0;
}

Real
MultiCompMultiPhaseBase::computeQpJacobian(){
  return 0;
}

Real
MultiCompMultiPhaseBase::computeQpOffDiagJacobian(unsigned int /*jvar*/){
  return 0;
}
