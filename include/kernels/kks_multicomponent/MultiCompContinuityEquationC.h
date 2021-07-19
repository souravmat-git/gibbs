//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MULTICOMPCONTINUITYEQUATIONC_H
#define MULTICOMPCONTINUITYEQUATIONC_H

#include "Kernel.h"

class MultiCompContinuityEquationC;

template <>
InputParameters validParams<MultiCompContinuityEquationC>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * The kernel assumes data for quartenary A-B-C-D alloy
 */
class MultiCompContinuityEquationC : public Kernel
{
public:
  MultiCompContinuityEquationC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:  
  const VariableGradient & _grad_B_diff_pot;
  unsigned int _grad_B_diff_pot_var;
  
  const VariableValue & _eta;
  unsigned int _eta_var; 
  
  //Composition variables
  const VariableValue & _xB_beta;
  unsigned int _xB_beta_var;
  
  const VariableValue & _xC_beta;
  unsigned int _xC_beta_var;
  
  const VariableValue & _xB_alpha;
  unsigned int _xB_alpha_var;
  
  const VariableValue & _xC_alpha;
  unsigned int _xC_alpha_var;
  
  //Diffusion potential variables
  const VariableGradient & _grad_D_diff_pot;
  unsigned int _grad_D_diff_pot_var;

/// Isotropic onsager mobility depdendent on concentration
  const MaterialProperty<Real> & _L_CC_beta;
  const MaterialProperty<Real> & _L_BC_beta;
  const MaterialProperty<Real> & _L_CD_beta;
  
  const MaterialProperty<Real> & _L_CC_alpha;
  const MaterialProperty<Real> & _L_BC_alpha;
  const MaterialProperty<Real> & _L_CD_alpha;
  
   //The derivatives of the Onsgare mobility with
  //respect to xB_beta &  xC_beta
  const MaterialProperty<Real> & _dL_CC_xB_beta;
  const MaterialProperty<Real> & _dL_BC_xB_beta;
  
  const MaterialProperty<Real> & _dL_CC_xC_beta;
  const MaterialProperty<Real> & _dL_BC_xC_beta;
  
  //Derivative of L with resco xB_alpha & xC_alpha
  const MaterialProperty<Real> & _dL_CC_xB_alpha;
  const MaterialProperty<Real> & _dL_BC_xB_alpha;
  
  const MaterialProperty<Real> & _dL_CC_xC_alpha;
  const MaterialProperty<Real> & _dL_BC_xC_alpha;
  
  //Interpolation function
  
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  //Real _L_CC_interp, _L_BC_interp, _L_CD_interp;
};
#endif // MULTICOMPONENTCONTINUITYEQUATIONC_H
