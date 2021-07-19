//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "MultiPhaseBase.h"

class MultiCompMultiPhaseBase;

template <>
InputParameters validParams<MultiCompMultiPhaseBase>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * variable on which this kernel operates: X
 */
class MultiCompMultiPhaseBase: public MultiPhaseBase
{
public:
  MultiCompMultiPhaseBase(const InputParameters & parameters);

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  //All Onsager Mobility to be interpolated
  Real L_BB_interp() const;
  Real L_BC_interp() const;
  Real L_BD_interp() const;
  Real L_CC_interp() const;
  Real L_CD_interp() const;
  Real L_DD_interp() const;
  
  Real dL_BB_muB_interp() const;
  Real dL_BC_muB_interp() const;
  Real dL_BD_muB_interp() const;
  Real dL_CC_muB_interp() const;
  Real dL_CD_muB_interp() const;
  Real dL_DD_muB_interp() const;
  
  Real dL_BB_muC_interp() const;
  Real dL_BC_muC_interp() const;
  Real dL_BD_muC_interp() const;
  Real dL_CC_muC_interp() const;
  Real dL_CD_muC_interp() const;
  Real dL_DD_muC_interp() const;
  
  Real dL_BB_muD_interp() const;
  Real dL_BC_muD_interp() const;
  Real dL_BD_muD_interp() const;
  Real dL_CC_muD_interp() const;
  Real dL_CD_muD_interp() const;
  Real dL_DD_muD_interp() const;
    
  //Overall chi matrix
  Real chi_BB() const;
  Real chi_BC() const;
  Real chi_BD() const;
  Real chi_CC() const;
  Real chi_CD() const;
  Real chi_DD() const;
  
  //Determinant of the overall susceptibility matrix
  Real det_chi() const;
  
  //Overall thermodynamic factors for a quartenary alloy
  Real thermodynamic_factorBB() const;
  Real thermodynamic_factorBC() const;
  Real thermodynamic_factorBD() const;
  Real thermodynamic_factorCC() const;
  Real thermodynamic_factorCD() const;
  Real thermodynamic_factorDD() const;
  
  //Chemical diffusivity matrix not symmetric
  Real DC_BB_interp() const;
  Real DC_BC_interp() const;
  Real DC_BD_interp() const;
  Real DC_CB_interp() const;
  Real DC_CC_interp() const;
  Real DC_CD_interp() const;
  Real DC_DB_interp() const;
  Real DC_DC_interp() const;
  Real DC_DD_interp() const;
  
  //first derivative of the interpolation function
  //and sum over all phases
  Real sum_dh_L_BB_alpha() const;
  Real sum_dh_L_BB_beta() const;
  Real sum_dh_L_BB_gamma() const;
  Real sum_dh_L_BB_delta() const;
  Real sum_dh_L_BB_epsilon() const;
  
  Real sum_dh_L_BC_alpha() const;
  Real sum_dh_L_BC_beta() const;
  Real sum_dh_L_BC_gamma() const;
  Real sum_dh_L_BC_delta() const;
  Real sum_dh_L_BC_epsilon() const;
  
  Real sum_dh_L_BD_alpha() const;
  Real sum_dh_L_BD_beta() const;
  Real sum_dh_L_BD_gamma() const;
  Real sum_dh_L_BD_delta() const;
  Real sum_dh_L_BD_epsilon() const;
  
  Real sum_dh_L_CC_alpha() const;
  Real sum_dh_L_CC_beta() const;
  Real sum_dh_L_CC_gamma() const;
  Real sum_dh_L_CC_delta() const;
  Real sum_dh_L_CC_epsilon() const;
  
  Real sum_dh_L_CD_alpha() const;
  Real sum_dh_L_CD_beta() const;
  Real sum_dh_L_CD_gamma() const;
  Real sum_dh_L_CD_delta() const;
  Real sum_dh_L_CD_epsilon() const;
   
  Real sum_dh_L_DD_alpha() const;
  Real sum_dh_L_DD_beta() const;
  Real sum_dh_L_DD_gamma() const;
  Real sum_dh_L_DD_delta() const;
  Real sum_dh_L_DD_epsilon() const;
  
  //Consequently each material property should depend on phase
  //First derivative xB with respect to comp B diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  const MaterialProperty<Real> & _inv_B_tf_beta;
  const MaterialProperty<Real> & _inv_B_tf_gamma;
  const MaterialProperty<Real> & _inv_B_tf_delta; 
  const MaterialProperty<Real> & _inv_B_tf_epsilon; 
  
  //First derivative xB with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_BC_tf_alpha;
  const MaterialProperty<Real> & _inv_BC_tf_beta;
  const MaterialProperty<Real> & _inv_BC_tf_gamma;
  const MaterialProperty<Real> & _inv_BC_tf_delta; 
  const MaterialProperty<Real> & _inv_BC_tf_epsilon; 
  
  //First derivative xB with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_BD_tf_alpha;
  const MaterialProperty<Real> & _inv_BD_tf_beta;
  const MaterialProperty<Real> & _inv_BD_tf_gamma;
  const MaterialProperty<Real> & _inv_BD_tf_delta; 
  const MaterialProperty<Real> & _inv_BD_tf_epsilon; 
  
  //First derivative xC with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_C_tf_alpha;
  const MaterialProperty<Real> & _inv_C_tf_beta;
  const MaterialProperty<Real> & _inv_C_tf_gamma;
  const MaterialProperty<Real> & _inv_C_tf_delta; 
  const MaterialProperty<Real> & _inv_C_tf_epsilon; 
  
  //First derivative xC with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_CD_tf_alpha;
  const MaterialProperty<Real> & _inv_CD_tf_beta;
  const MaterialProperty<Real> & _inv_CD_tf_gamma;
  const MaterialProperty<Real> & _inv_CD_tf_delta; 
  const MaterialProperty<Real> & _inv_CD_tf_epsilon; 
  
   //First derivative xD with respect to comp D diffusion potential
  const MaterialProperty<Real> & _inv_D_tf_alpha;
  const MaterialProperty<Real> & _inv_D_tf_beta;
  const MaterialProperty<Real> & _inv_D_tf_gamma;
  const MaterialProperty<Real> & _inv_D_tf_delta; 
  const MaterialProperty<Real> & _inv_D_tf_epsilon;
  
  // All phase dependent material property L_BB
  const MaterialProperty<Real>& _L_BB_alpha; 
  const MaterialProperty<Real>& _L_BB_beta; 
  const MaterialProperty<Real>& _L_BB_gamma;
  const MaterialProperty<Real>& _L_BB_delta;
  const MaterialProperty<Real>& _L_BB_epsilon;
  
  // All phase dependent material property L_BC
  const MaterialProperty<Real>& _L_BC_alpha; 
  const MaterialProperty<Real>& _L_BC_beta; 
  const MaterialProperty<Real>& _L_BC_gamma;
  const MaterialProperty<Real>& _L_BC_delta;
  const MaterialProperty<Real>& _L_BC_epsilon;
  
  // All phase dependent material property L_BD
  const MaterialProperty<Real>& _L_BD_alpha; 
  const MaterialProperty<Real>& _L_BD_beta; 
  const MaterialProperty<Real>& _L_BD_gamma;
  const MaterialProperty<Real>& _L_BD_delta;
  const MaterialProperty<Real>& _L_BD_epsilon;
  
  // All phase dependent material property L_CC
  const MaterialProperty<Real>& _L_CC_alpha; 
  const MaterialProperty<Real>& _L_CC_beta; 
  const MaterialProperty<Real>& _L_CC_gamma;
  const MaterialProperty<Real>& _L_CC_delta;
  const MaterialProperty<Real>& _L_CC_epsilon;
  
  // All phase dependent material property L_CD
  const MaterialProperty<Real>& _L_CD_alpha; 
  const MaterialProperty<Real>& _L_CD_beta; 
  const MaterialProperty<Real>& _L_CD_gamma;
  const MaterialProperty<Real>& _L_CD_delta;
  const MaterialProperty<Real>& _L_CD_epsilon;
  
   // All phase dependent material property L_DD
  const MaterialProperty<Real>& _L_DD_alpha; 
  const MaterialProperty<Real>& _L_DD_beta; 
  const MaterialProperty<Real>& _L_DD_gamma;
  const MaterialProperty<Real>& _L_DD_delta;
  const MaterialProperty<Real>& _L_DD_epsilon;
 
  //***********************All derivatives with respect to B*****************/
  //dL_BB_muB
  const MaterialProperty<Real>& _dL_BB_muB_alpha;
  const MaterialProperty<Real>& _dL_BB_muB_beta; 
  const MaterialProperty<Real>& _dL_BB_muB_gamma;
  const MaterialProperty<Real>& _dL_BB_muB_delta;
  const MaterialProperty<Real>& _dL_BB_muB_epsilon;
  
  //dL_BC_muB
  const MaterialProperty<Real>& _dL_BC_muB_alpha;
  const MaterialProperty<Real>& _dL_BC_muB_beta; 
  const MaterialProperty<Real>& _dL_BC_muB_gamma;
  const MaterialProperty<Real>& _dL_BC_muB_delta;
  const MaterialProperty<Real>& _dL_BC_muB_epsilon;  
  
  //dL_BD_muB
  const MaterialProperty<Real>& _dL_BD_muB_alpha;
  const MaterialProperty<Real>& _dL_BD_muB_beta; 
  const MaterialProperty<Real>& _dL_BD_muB_gamma;
  const MaterialProperty<Real>& _dL_BD_muB_delta;
  const MaterialProperty<Real>& _dL_BD_muB_epsilon; 
  
  //dL_CC_muB
  const MaterialProperty<Real>& _dL_CC_muB_alpha;
  const MaterialProperty<Real>& _dL_CC_muB_beta; 
  const MaterialProperty<Real>& _dL_CC_muB_gamma;
  const MaterialProperty<Real>& _dL_CC_muB_delta;
  const MaterialProperty<Real>& _dL_CC_muB_epsilon;
  
  //dL_BC_muB
  const MaterialProperty<Real>& _dL_CD_muB_alpha;
  const MaterialProperty<Real>& _dL_CD_muB_beta; 
  const MaterialProperty<Real>& _dL_CD_muB_gamma;
  const MaterialProperty<Real>& _dL_CD_muB_delta;
  const MaterialProperty<Real>& _dL_CD_muB_epsilon;  
  
  //dL_BD_muB
  const MaterialProperty<Real>& _dL_DD_muB_alpha;
  const MaterialProperty<Real>& _dL_DD_muB_beta; 
  const MaterialProperty<Real>& _dL_DD_muB_gamma;
  const MaterialProperty<Real>& _dL_DD_muB_delta;
  const MaterialProperty<Real>& _dL_DD_muB_epsilon;   
  
  //***********************All derivatives with respect to C*****************/
  //dL_BB_muC
  const MaterialProperty<Real>& _dL_BB_muC_alpha;
  const MaterialProperty<Real>& _dL_BB_muC_beta; 
  const MaterialProperty<Real>& _dL_BB_muC_gamma;
  const MaterialProperty<Real>& _dL_BB_muC_delta;
  const MaterialProperty<Real>& _dL_BB_muC_epsilon;
  
  //dL_BC_muC
  const MaterialProperty<Real>& _dL_BC_muC_alpha;
  const MaterialProperty<Real>& _dL_BC_muC_beta; 
  const MaterialProperty<Real>& _dL_BC_muC_gamma;
  const MaterialProperty<Real>& _dL_BC_muC_delta;
  const MaterialProperty<Real>& _dL_BC_muC_epsilon;  
  
  //dL_BD_muC
  const MaterialProperty<Real>& _dL_BD_muC_alpha;
  const MaterialProperty<Real>& _dL_BD_muC_beta; 
  const MaterialProperty<Real>& _dL_BD_muC_gamma;
  const MaterialProperty<Real>& _dL_BD_muC_delta;
  const MaterialProperty<Real>& _dL_BD_muC_epsilon; 
  
  //dL_CC_muC
  const MaterialProperty<Real>& _dL_CC_muC_alpha;
  const MaterialProperty<Real>& _dL_CC_muC_beta; 
  const MaterialProperty<Real>& _dL_CC_muC_gamma;
  const MaterialProperty<Real>& _dL_CC_muC_delta;
  const MaterialProperty<Real>& _dL_CC_muC_epsilon;
  
  //dL_CD_muC
  const MaterialProperty<Real>& _dL_CD_muC_alpha;
  const MaterialProperty<Real>& _dL_CD_muC_beta; 
  const MaterialProperty<Real>& _dL_CD_muC_gamma;
  const MaterialProperty<Real>& _dL_CD_muC_delta;
  const MaterialProperty<Real>& _dL_CD_muC_epsilon;  
  
  //dL_DD_muC
  const MaterialProperty<Real>& _dL_DD_muC_alpha;
  const MaterialProperty<Real>& _dL_DD_muC_beta; 
  const MaterialProperty<Real>& _dL_DD_muC_gamma;
  const MaterialProperty<Real>& _dL_DD_muC_delta;
  const MaterialProperty<Real>& _dL_DD_muC_epsilon;  
  
  //***********************All derivatives with respect to D*****************/
  //dL_BB_muD
  const MaterialProperty<Real>& _dL_BB_muD_alpha;
  const MaterialProperty<Real>& _dL_BB_muD_beta; 
  const MaterialProperty<Real>& _dL_BB_muD_gamma;
  const MaterialProperty<Real>& _dL_BB_muD_delta;
  const MaterialProperty<Real>& _dL_BB_muD_epsilon;
  
  //dL_BC_muD
  const MaterialProperty<Real>& _dL_BC_muD_alpha;
  const MaterialProperty<Real>& _dL_BC_muD_beta; 
  const MaterialProperty<Real>& _dL_BC_muD_gamma;
  const MaterialProperty<Real>& _dL_BC_muD_delta;
  const MaterialProperty<Real>& _dL_BC_muD_epsilon;  
  
  //dL_BD_muD
  const MaterialProperty<Real>& _dL_BD_muD_alpha;
  const MaterialProperty<Real>& _dL_BD_muD_beta; 
  const MaterialProperty<Real>& _dL_BD_muD_gamma;
  const MaterialProperty<Real>& _dL_BD_muD_delta;
  const MaterialProperty<Real>& _dL_BD_muD_epsilon; 
  
  //dL_CC_muD
  const MaterialProperty<Real>& _dL_CC_muD_alpha;
  const MaterialProperty<Real>& _dL_CC_muD_beta; 
  const MaterialProperty<Real>& _dL_CC_muD_gamma;
  const MaterialProperty<Real>& _dL_CC_muD_delta;
  const MaterialProperty<Real>& _dL_CC_muD_epsilon;
  
  //dL_CD_muD
  const MaterialProperty<Real>& _dL_CD_muD_alpha;
  const MaterialProperty<Real>& _dL_CD_muD_beta; 
  const MaterialProperty<Real>& _dL_CD_muD_gamma;
  const MaterialProperty<Real>& _dL_CD_muD_delta;
  const MaterialProperty<Real>& _dL_CD_muD_epsilon;  
  
  //dL_DD_muD
  const MaterialProperty<Real>& _dL_DD_muD_alpha;
  const MaterialProperty<Real>& _dL_DD_muD_beta; 
  const MaterialProperty<Real>& _dL_DD_muD_gamma;
  const MaterialProperty<Real>& _dL_DD_muD_delta;
  const MaterialProperty<Real>& _dL_DD_muD_epsilon;    
  
};
