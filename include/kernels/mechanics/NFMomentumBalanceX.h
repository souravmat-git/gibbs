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

class NFMomentumBalanceX;

template <>
InputParameters validParams<NFMomentumBalanceX>();

/**
 * Kernel to implement momentum balance in two phase
 * material
 */
class NFMomentumBalanceX : public Kernel
{
public:
  NFMomentumBalanceX(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
private:

  const VariableValue & _ex_beta;
  unsigned int _ex_beta_var;

  //Stress in the beta phase and the stress in the alpha phase
  const MaterialProperty<Real> & _sx_beta;

 //material constant (lambda + 2mu) in the beta phase and the alpha phase
  const MaterialProperty<Real> & _mat_const_beta;

};
