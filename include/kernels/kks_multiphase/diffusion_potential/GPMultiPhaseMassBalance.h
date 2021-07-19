//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef GPMULTIPHASEMASSBALANCE_H
//#define GPMULTIPHASEMASSBALANCE_H

#pragma once

#include "Kernel.h"

class GPMultiPhaseMassBalance;

template <>
InputParameters validParams<GPMultiPhaseMassBalance>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * variable on which this kernel operates: mu
 */
class GPMultiPhaseMassBalance: public Kernel
{
public:
  GPMultiPhaseMassBalance(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:  

  //Declare constant memeber functon
  Real L_BB_interp() const;
  Real dL_BB_muB_interp() const;
 
  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;
  
  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;

  std::string _h_alpha_name, _h_beta_name, _h_gamma_name;
  
  const MaterialProperty<Real>& _h_alpha;
  const MaterialProperty<Real>& _h_beta;
  const MaterialProperty<Real>& _h_gamma;
  
  //First derivatives of the interpolation function
  const MaterialProperty<Real>& _dhbeta_dphialpha;
  
  const MaterialProperty<Real>& _dhalpha_dphibeta;  
  
  /// Isotropic diffusion mobility dependent on diffusion potential
  const MaterialProperty<Real> & _L_BB_beta;
  
  const MaterialProperty<Real> & _L_BB_alpha;
  
  const MaterialProperty<Real> & _dL_BB_muB_beta;
  
  const MaterialProperty<Real> & _dL_BB_muB_alpha;
};

//#endif // GPMULTIPHASEMASSBALANCE_H
