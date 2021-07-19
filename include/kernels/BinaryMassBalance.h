//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef BINARYMASSBALANCE_H
#define BINARYMASSBALANCE_H

#include "Kernel.h"

class BinaryMassBalance;

template <>
InputParameters validParams<BinaryMassBalance>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * variable on which this kernel operates: X
 */
class BinaryMassBalance: public Kernel
{
public:
  BinaryMassBalance(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:  

  //Declare constant memeber functon
  Real thermodynamic_factor() const;
  Real L_BB_interp() const;
  Real dL_BB_muB_interp() const;
 
  const VariableValue & _eta;
  unsigned int _eta_var;
  
  const VariableGradient & _grad_B_diff_pot;
  unsigned int _B_diff_pot_var;

  //First derivative of composition with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  const MaterialProperty<Real> & _inv_B_tf_beta;
  
  //_h is a reference variable that indicates the interpolation function
  const MaterialProperty<Real>& _h;
  const MaterialProperty<Real>& _dh; 
  
  /// Isotropic diffusion mobility dependent on diffusion potential
  const MaterialProperty<Real> & _L_BB_beta;
  const MaterialProperty<Real> & _L_BB_alpha;
  
  const MaterialProperty<Real> & _dL_BB_muB_beta;
  const MaterialProperty<Real> & _dL_BB_muB_alpha;
};

#endif // BINARYMASSBALANCE_H
