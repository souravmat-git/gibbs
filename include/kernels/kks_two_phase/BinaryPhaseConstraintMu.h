//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef BINARYPHASECONSTRAINTMU_H
#define BINARYPHASECONSTRAINTMU_H

//Include dependencies
#include "Kernel.h"

// Forward Declarations
class BinaryPhaseConstraintMu;

template <>
InputParameters validParams<BinaryPhaseConstraintMu>();
/**
 * Enforce sum of phase concentrations to be the real concentration.
 * The non-linear variable that this kernel acts on
 * is the mole fraction of component B
 * This kernel will be coupled to diff_comp_C, eta, mole_fraction
 **/
class BinaryPhaseConstraintMu : public Kernel
{
public:
  BinaryPhaseConstraintMu(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  const VariableValue & _eta;
  unsigned int _eta_var;
  
  const VariableValue & _xB;
  unsigned int _xB_var;
  
  //phase comp_alpha is a function of diffusion potential
  const MaterialProperty<Real> & _xB_alpha;
  
  //Phase comp_beta is a function of diffusion potential
  const MaterialProperty<Real> & _xB_beta;
  
  //First derivative with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  
  //First derivative with respect to diffusion potential
  const MaterialProperty<Real> &  _inv_B_tf_beta;
  
  //_h is a reference variable that indicates the interpolation function
  const MaterialProperty<Real>& _h;
  // _dh is a reference variable that indicates the first derivative of the 
  // interpolation function
  const MaterialProperty<Real>& _dh; 
};
#endif // BINARYPHASECONSTRAINT_H
