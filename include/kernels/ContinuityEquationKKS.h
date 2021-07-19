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

class ContinuityEquationKKS;

template <>
InputParameters validParams<ContinuityEquationKKS>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * based on the KKS model
 */
class ContinuityEquationKKS : public Kernel
{
public:
  ContinuityEquationKKS(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:

  //Declare constant member function
  //Real f_BB() const;  //See Eqn 29 of Kim's paper f_{cc}
  Real L_BB_interp() const;
  Real dL_BB_xB_beta_interp() const;
  
  const VariableValue & _eta;
  unsigned int _eta_var;
  
  //The overall mole fraction is also a coupled variable
  //const VariableGradient & _grad_xB;
  //unsigned int _xB_var;
  
  const VariableGradient &  _grad_xB_beta;
  unsigned int _xB_beta_var;
  
  //Thermodynamic factor of only one phase is required
  //This can be either beta or alpha phase but never both
  //Here, phase beta is assumed
  //const MaterialProperty<Real> & _B_tf_alpha;
  const MaterialProperty<Real> & _B_tf_beta;
  
  //Interpolation functions
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  /// Isotropic diffusion mobility dependent of concentration
  const MaterialProperty<Real> & _L_BB_alpha;
  const MaterialProperty<Real> & _L_BB_beta;
  
  //Note that the properries are dependent on phase compositions  
  const MaterialProperty<Real> & _dL_BB_xB_alpha;
  const MaterialProperty<Real> & _dL_BB_xB_beta;
};

