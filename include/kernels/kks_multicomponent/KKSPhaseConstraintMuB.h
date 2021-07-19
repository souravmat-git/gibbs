//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

//Include dependencies
#include "MultiCompMultiPhaseBase.h"

// Forward Declarations
class KKSPhaseConstraintMuB;

template <>
InputParameters validParams<KKSPhaseConstraintMuB>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 * The non-linear variable that this kernel acts on is
 * the diffusion potential of component B
 * This kernel is coupled to the mole fraction of B
 * through the continuity equation
 *
 */
class KKSPhaseConstraintMuB : public MultiCompMultiPhaseBase
{
public:
  KKSPhaseConstraintMuB(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  //Acts on the component B of alloy system A-B-C-D
  //Non-linear variable = B_diff_pot
  
  const VariableValue & _xB;
  unsigned int _xB_var;
  
  //xC_alpha is a function of diffusion potentials B,C,D
  const MaterialProperty<Real> & _xB_alpha;
  const MaterialProperty<Real> & _xB_beta;
  const MaterialProperty<Real> & _xB_gamma;
  const MaterialProperty<Real> & _xB_delta;
  const MaterialProperty<Real> & _xB_epsilon;
  
  //Diffusion potential of component C is a coupled 
  //variable to this kernel for a ternary alloy A-B-C
  const VariableValue & _C_diff_pot;
  unsigned int _C_diff_pot_var;
  
  //For a quaternary alloy A-B-C-D
  //Diffusion potential of component D is a coupled 
  //variable to this kernel
  const VariableValue & _D_diff_pot;
  unsigned int _D_diff_pot_var;
};
