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
class MCPhaseConstraintMuD;

template <>
InputParameters validParams<MCPhaseConstraintMuD>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 * The non-linear variable that this kernel acts on
 * is the diffusion potential of component D
 * This kernel will be coupled to diff_comp_B, eta, mole_fraction,
 * diff_comp_C
 *
 */
class MCPhaseConstraintMuD : public MultiPhaseBase
{
public:
  MCPhaseConstraintMuD(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  Real chi_CD() const;
  Real chi_BD() const;
  Real chi_DD() const;
  
  //For a quaternary alloy A-B-C-D
  //Mole fraction of component D is a coupled 
  //variable to this kernel
  const VariableValue & _xD;
  unsigned int _xD_var;
  
  //xD_alpha is a function of diffusion potentials B,C,D
  const MaterialProperty<Real> & _xD_alpha;
  const MaterialProperty<Real> & _xD_beta;
  const MaterialProperty<Real> & _xD_gamma;
  const MaterialProperty<Real> & _xD_delta;
  const MaterialProperty<Real> & _xD_epsilon;
    
  //First derivative xD with respect to comp C diffusion potential
  const MaterialProperty<Real> & _inv_CD_tf_alpha;
  const MaterialProperty<Real> & _inv_CD_tf_beta;
  const MaterialProperty<Real> & _inv_CD_tf_gamma;
  const MaterialProperty<Real> & _inv_CD_tf_delta; 
  const MaterialProperty<Real> & _inv_CD_tf_epsilon; 
  
  //First derivative xD with respect to comp B diffusion potential
  const MaterialProperty<Real> & _inv_BD_tf_alpha;
  const MaterialProperty<Real> & _inv_BD_tf_beta;
  const MaterialProperty<Real> & _inv_BD_tf_gamma;
  const MaterialProperty<Real> & _inv_BD_tf_delta; 
  const MaterialProperty<Real> & _inv_BD_tf_epsilon;
  
  //First derivative xD with respect to comp D diffusion potential
  const MaterialProperty<Real> & _inv_D_tf_alpha;
  const MaterialProperty<Real> & _inv_D_tf_beta;
  const MaterialProperty<Real> & _inv_D_tf_gamma;
  const MaterialProperty<Real> & _inv_D_tf_delta; 
  const MaterialProperty<Real> & _inv_D_tf_epsilon; 

  //Diffusion potential of component C is a coupled 
  //variable to this kernel for a quaternary alloy A-B-C-D
  const VariableValue & _C_diff_pot;
  unsigned int _C_diff_pot_var;
  
  //Diffusion potential of component B is a coupled 
  //variable to this kernel for a quaternary alloy A-B-C-D   
  const VariableValue & _B_diff_pot;
  unsigned int _B_diff_pot_var; 
};
