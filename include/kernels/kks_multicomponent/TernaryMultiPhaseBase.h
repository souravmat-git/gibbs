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

class TernaryMultiPhaseBase;

template <>
InputParameters validParams<TernaryMultiPhaseBase>();

/**
 * Kernel to precompute
 * the a)overall diffusivity matrix
 * and b)overall Onsager mobility and its derivatives
 * for a ternary A-B-C alloy
 */
class TernaryMultiPhaseBase: public MultiPhaseBase
{
public:
  TernaryMultiPhaseBase(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  //Overall chi matrix
  Real chi_BB() const;
  Real chi_BC() const;
  Real chi_CC() const;
  
  //Determinant of the overall susceptibility matrix
  Real det_chi() const;
  
  //Overall thermodynamic factors for a quartenary alloy
  Real thermodynamic_factorBB() const;
  Real thermodynamic_factorBC() const;
  Real thermodynamic_factorCC() const;

  //All Onsager Mobility to be interpolated
  Real L_BB_interp() const;
  Real L_BC_interp() const;
  Real L_CC_interp() const;
  
  //Chemical diffusivity matrix not symmetric
  Real DC_BB_interp() const;
  Real DC_BC_interp() const;
  Real DC_CB_interp() const;
  Real DC_CC_interp() const;
  
  //Its first derivatives
  Real dL_BB_muB_interp() const;
  Real dL_BC_muB_interp() const;
  Real dL_CC_muB_interp() const;
  
  Real dL_BB_muC_interp() const;
  Real dL_BC_muC_interp() const;
  Real dL_CC_muC_interp() const; 
  
  //first derivative of the interpolation function
  //and sum over all phases
  Real sum_dh_L_BB_alpha() const;
  Real sum_dh_L_BB_beta()  const;
  Real sum_dh_L_BB_gamma() const;
  Real sum_dh_L_BB_delta() const;
  Real sum_dh_L_BB_epsilon() const;
  
  Real sum_dh_L_BC_alpha() const;
  Real sum_dh_L_BC_beta() const;
  Real sum_dh_L_BC_gamma() const;
  Real sum_dh_L_BC_delta() const;
  Real sum_dh_L_BC_epsilon() const;
    
  Real sum_dh_L_CC_alpha() const;
  Real sum_dh_L_CC_beta() const;
  Real sum_dh_L_CC_gamma() const;
  Real sum_dh_L_CC_delta() const;
  Real sum_dh_L_CC_epsilon() const;
  
  //This base class requires three material property
  //1) The chi matrix of each phase of the system
  //2) The Onsager matrix of each phase
  //Based on this it constructs two derived material property
  //1) The interpolated thermodynamic factor matrix of the overall system
  //2) The interpolated diffusivity matrix of the overall system
  
  //Note that the coeffecints of the susceptibility matrix are here referred
  // to as the inverse of the tf matrix.
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
  
  //First derivative xC with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_C_tf_alpha;
  const MaterialProperty<Real> & _inv_C_tf_beta;
  const MaterialProperty<Real> & _inv_C_tf_gamma;
  const MaterialProperty<Real> & _inv_C_tf_delta; 
  const MaterialProperty<Real> & _inv_C_tf_epsilon; 
  
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
  
  // All phase dependent material property L_CC
  const MaterialProperty<Real>& _L_CC_alpha; 
  const MaterialProperty<Real>& _L_CC_beta; 
  const MaterialProperty<Real>& _L_CC_gamma;
  const MaterialProperty<Real>& _L_CC_delta;
  const MaterialProperty<Real>& _L_CC_epsilon;
  
  //To model the diffusion potential dependence of the Onsager mobility 
  //matrix this is required

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
 
  //dL_CC_muB
  const MaterialProperty<Real>& _dL_CC_muB_alpha;
  const MaterialProperty<Real>& _dL_CC_muB_beta; 
  const MaterialProperty<Real>& _dL_CC_muB_gamma;
  const MaterialProperty<Real>& _dL_CC_muB_delta;
  const MaterialProperty<Real>& _dL_CC_muB_epsilon;
   
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
  
  //dL_CC_muC
  const MaterialProperty<Real>& _dL_CC_muC_alpha;
  const MaterialProperty<Real>& _dL_CC_muC_beta; 
  const MaterialProperty<Real>& _dL_CC_muC_gamma;
  const MaterialProperty<Real>& _dL_CC_muC_delta;
  const MaterialProperty<Real>& _dL_CC_muC_epsilon;
 
};
