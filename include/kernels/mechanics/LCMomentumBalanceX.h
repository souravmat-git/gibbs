//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

//MOOSE includes
#include "Kernel.h"

class LCMomentumBalanceX;

template <>
InputParameters validParams<LCMomentumBalanceX>();
/**
 * Kernel to implement momentum balance in two phase
 * material
 **/
class LCMomentumBalanceX : public Kernel
{
public:
  LCMomentumBalanceX(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  Real compliance_inv() const;
  Real dcomplinv_dh()   const;
  Real d2complinv_dh2()  const;
  
  //Phase-field is a coupled variable to this kernel
  const VariableGradient & _grad_eta;
  unsigned int _eta_var;
  
  //The overall stress is also a coupled variable
  const VariableValue & _stress_x;
  unsigned int _stress_x_var;
  
  //Strains as a function of stress
  const MaterialProperty<Real> & _ex_beta;
  const MaterialProperty<Real> & _ex_alpha;
     
  //Compliance of beta and alpha phase
  const MaterialProperty<Real> & _comp_beta;
  const MaterialProperty<Real> & _comp_alpha;
  
  //Interpolation function
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;

};
