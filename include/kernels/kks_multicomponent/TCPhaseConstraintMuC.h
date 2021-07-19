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
class TCPhaseConstraintMuC;

template <>
InputParameters validParams<TCPhaseConstraintMuC>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 * The non-linear variable that this kernel acts on
 * is the mole fraction of comp C
 * This kernel will be coupled to diff_comp_B, eta, mole_fraction,
 * diff_comp_D
 *
 */
class TCPhaseConstraintMuC : public MultiPhaseBase
{
public:
  TCPhaseConstraintMuC(const InputParameters & parameters);

protected:
  // The override command ensures that a virtual member function
  // from a base class is overridden  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  Real chi_BC() const;
  Real chi_CC() const;

  //Acts on the component C of alloy system A-B-C
  //Non-linear variable = C_diff_pot
  
  const VariableValue & _xC;
  unsigned int _xC_var;
  
  //xC_alpha is a function of diffusion potentials B,C,D
  const MaterialProperty<Real> & _xC_alpha;
  const MaterialProperty<Real> & _xC_beta;
  const MaterialProperty<Real> & _xC_gamma;
  const MaterialProperty<Real> & _xC_delta;
  const MaterialProperty<Real> & _xC_epsilon;
  
  //First derivative xB with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_BC_tf_alpha;
  const MaterialProperty<Real> & _inv_BC_tf_beta;
  const MaterialProperty<Real> & _inv_BC_tf_gamma;
  const MaterialProperty<Real> & _inv_BC_tf_delta; 
  const MaterialProperty<Real> & _inv_BC_tf_epsilon; 
  
  //First derivative xC with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_C_tf_alpha;
  const MaterialProperty<Real> & _inv_C_tf_beta;
  const MaterialProperty<Real> & _inv_C_tf_gamma;
  const MaterialProperty<Real> & _inv_C_tf_delta; 
  const MaterialProperty<Real> & _inv_C_tf_epsilon; 
   
  //Diffusion potential of component B is a coupled 
  //variable to this kernel for a ternary alloy A-B-C
  const VariableValue & _B_diff_pot;
  unsigned int _B_diff_pot_var;
    
  //For a quaternary alloy A-B-C-D
  //Diffusion potential of component D is a coupled 
  //variable to this kernel
  //const VariableValue & _D_diff_pot;
  //unsigned int _D_diff_pot_var;
};
