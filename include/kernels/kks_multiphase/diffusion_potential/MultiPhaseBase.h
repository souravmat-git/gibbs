//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Kernel.h"

class MultiPhaseBase;

template <>
InputParameters validParams<MultiPhaseBase>();

/**
 * This kernel is the base class for multiphase system
 * Later this can be inherited by binary, ternary, and quarternary systems
 * equation for mass conservation
 * variable on which this kernel operates: X
 */
class MultiPhaseBase: public Kernel
{
public:
  MultiPhaseBase(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  //All phases for any multi-phase system
  //Here, we have assumed 5 phases;
  //It can be extended to N..phases;
 
  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;
  
  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;
  
  const VariableValue & _phase_gamma;
  unsigned int _phase_gamma_var;
  
  const VariableValue & _phase_delta;
  unsigned int _phase_delta_var;
  
  const VariableValue & _phase_epsilon;
  unsigned int _phase_epsilon_var;
 

  std::string _h_alpha_name, _h_beta_name, 
              _h_gamma_name, _h_delta_name, _h_epsilon_name;
  
  //_h_alpha, _h_beta....
  const MaterialProperty<Real>& _h_alpha;
  const MaterialProperty<Real>& _h_beta;
  const MaterialProperty<Real>& _h_gamma;
  const MaterialProperty<Real>& _h_delta;
  const MaterialProperty<Real>& _h_epsilon;
  
  //First derivatives of the interpolation function for phase alpha_var
  const MaterialProperty<Real>& _dhbeta_dphialpha;
  const MaterialProperty<Real>& _dhgamma_dphialpha;
  const MaterialProperty<Real>& _dhdelta_dphialpha;
  const MaterialProperty<Real>& _dhepsilon_dphialpha;
  
  //First derivatives of the interpolation function for phase beta_var
  const MaterialProperty<Real>& _dhalpha_dphibeta;
  const MaterialProperty<Real>& _dhgamma_dphibeta;
  const MaterialProperty<Real>& _dhdelta_dphibeta;
  const MaterialProperty<Real>& _dhepsilon_dphibeta;
  
  //First derivatives of the interpolation function for phase gamma_var
  const MaterialProperty<Real>& _dhalpha_dphigamma;
  const MaterialProperty<Real>& _dhbeta_dphigamma;
  const MaterialProperty<Real>& _dhdelta_dphigamma;
  const MaterialProperty<Real>& _dhepsilon_dphigamma;
  
  //First derivatives of the interpolation function for phase delta_var
  const MaterialProperty<Real>& _dhalpha_dphidelta;
  const MaterialProperty<Real>& _dhbeta_dphidelta;
  const MaterialProperty<Real>& _dhgamma_dphidelta;
  const MaterialProperty<Real>& _dhepsilon_dphidelta;
  
  //First derivatives of the interpolation function for phase epsilon_var
  const MaterialProperty<Real>& _dhalpha_dphiepsilon;
  const MaterialProperty<Real>& _dhbeta_dphiepsilon;
  const MaterialProperty<Real>& _dhgamma_dphiepsilon;
  const MaterialProperty<Real>& _dhdelta_dphiepsilon;
  
};
