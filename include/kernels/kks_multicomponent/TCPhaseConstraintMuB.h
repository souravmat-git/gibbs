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
#include "MultiPhaseBase.h"

// Forward Declarations
class TCPhaseConstraintMuB;

template <>
InputParameters validParams<TCPhaseConstraintMuB>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 * The non-linear variable that this kernel acts on is
 * the diffusion potential of component B
 * This kernel is coupled to the mole fraction of B
 * through the continuity equation
 *
 */
class TCPhaseConstraintMuB : public MultiPhaseBase
{
public:
  TCPhaseConstraintMuB(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real chi_BB() const;
  Real chi_BC() const;
  
  //Acts on the component B of alloy system A-B-C
  //Non-linear variable = B_diff_pot
  
  const VariableValue & _xB;
  unsigned int _xB_var;
  
  //xC_alpha is a function of diffusion potentials B,C,D
  const MaterialProperty<Real> & _xB_alpha;
  const MaterialProperty<Real> & _xB_beta;
  const MaterialProperty<Real> & _xB_gamma;
  const MaterialProperty<Real> & _xB_delta;
  const MaterialProperty<Real> & _xB_epsilon;
  
  //Consequently each material property should depend on phase
  //First derivative xB with respect to comp B diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  const MaterialProperty<Real> & _inv_B_tf_beta;
  const MaterialProperty<Real> & _inv_B_tf_gamma;
  const MaterialProperty<Real> & _inv_B_tf_delta; 
  const MaterialProperty<Real> & _inv_B_tf_epsilon; 
  
  //First derivative xB with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_BC_tf_alpha;
  const MaterialProperty<Real> & _inv_BC_tf_beta;
  const MaterialProperty<Real> & _inv_BC_tf_gamma;
  const MaterialProperty<Real> & _inv_BC_tf_delta; 
  const MaterialProperty<Real> & _inv_BC_tf_epsilon; 
  
  //Diffusion potential of component C is a coupled 
  //variable to this kernel for a ternary alloy A-B-C
  const VariableValue & _C_diff_pot;
  unsigned int _C_diff_pot_var;
  
  //For a quaternary alloy A-B-C-D
  //Diffusion potential of component D is a coupled 
  //variable to this kernel
  //const VariableValue & _D_diff_pot;
  //unsigned int _D_diff_pot_var;
};
