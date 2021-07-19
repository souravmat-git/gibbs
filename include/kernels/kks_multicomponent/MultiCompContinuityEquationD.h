//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MULTICOMPCONTINUITYEQUATIOND_H
#define MULTICOMPCONTINUITYEQUATIOND_H

#include "Kernel.h"

class MultiCompContinuityEquationD;

template <>
InputParameters validParams<MultiCompContinuityEquationD>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * The kernel assumes data for quartenary A-B-C-D alloy
 */
class MultiCompContinuityEquationD : public Kernel
{
public:
  MultiCompContinuityEquationD(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:  
  const VariableGradient & _grad_C_diff_pot;
  unsigned int _grad_C_diff_pot_var;
  
  const VariableGradient & _grad_B_diff_pot;
  unsigned int _grad_B_diff_pot_var;
  
  const VariableValue & _eta;
  unsigned int _eta_var;

/// Isotropic onsager mobility dependent of concentration
  const MaterialProperty<Real> & _L_DD_beta;
  const MaterialProperty<Real> & _L_CD_beta;
  const MaterialProperty<Real> & _L_BD_beta;
  
  const MaterialProperty<Real> & _L_DD_alpha;
  const MaterialProperty<Real> & _L_CD_alpha;
  const MaterialProperty<Real> & _L_BD_alpha;
  
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  Real _L_DD_interp, _L_CD_interp, _L_BD_interp;
};
#endif // MULTICOMPONENTCONTINUITYEQUATIONC_H
