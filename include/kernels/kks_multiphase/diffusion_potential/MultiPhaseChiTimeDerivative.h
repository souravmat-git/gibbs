//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeDerivative.h"
// Forward Declaration
class MultiPhaseChiTimeDerivative;

template <>
InputParameters validParams<MultiPhaseChiTimeDerivative>();
/**
 * This calculates the time derivative for a variable multiplied by a generalized susceptibility
 **/
class MultiPhaseChiTimeDerivative : public TimeDerivative
{
public:
  MultiPhaseChiTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  Real interpolated_chi() const;
  Real third_deriv() const;

  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;
  
  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;
  
  const VariableValue & _phase_gamma;
  unsigned int _phase_gamma_var;
  
  
  //First derivatives with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  const MaterialProperty<Real> &  _inv_B_tf_beta;
  const MaterialProperty<Real> &  _inv_B_tf_gamma;
     
  //Second derivatives of composition with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_td_alpha;
  const MaterialProperty<Real> & _inv_B_td_beta;
  const MaterialProperty<Real> & _inv_B_td_gamma;
  
  //interpolation function
  const MaterialProperty<Real> & _h_alpha;
  const MaterialProperty<Real> & _h_beta;
  const MaterialProperty<Real> & _h_gamma;    
 
  //first derivative of interpolation function               
  const MaterialProperty<Real> & _dhbeta_dphialpha;
  const MaterialProperty<Real> & _dhgamma_dphialpha;
        
  const MaterialProperty<Real> & _dhalpha_dphibeta;
  const MaterialProperty<Real> & _dhgamma_dphibeta;
        
  const MaterialProperty<Real> & _dhalpha_dphigamma;
  const MaterialProperty<Real> & _dhbeta_dphigamma;
  
 
};
