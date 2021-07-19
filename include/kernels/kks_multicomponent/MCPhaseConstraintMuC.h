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
#include "TCPhaseConstraintMuC.h"

// Forward Declarations
class MCPhaseConstraintMuC;

template <>
InputParameters validParams<MCPhaseConstraintMuC>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 * The non-linear variable that this kernel acts on is
 * the diffusion potential of component B for a A-B-C-D alloy
 * This kernel is coupled to the mole fraction of C
 * through the continuity equation
 *
 */
class MCPhaseConstraintMuC : public TCPhaseConstraintMuC
{
public:
  MCPhaseConstraintMuC(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  Real chi_CD() const;
  
  //First derivative xC with respect to comp D diffusion potential
  const MaterialProperty<Real> & _inv_CD_tf_alpha;
  const MaterialProperty<Real> & _inv_CD_tf_beta;
  const MaterialProperty<Real> & _inv_CD_tf_gamma;
  const MaterialProperty<Real> & _inv_CD_tf_delta; 
  const MaterialProperty<Real> & _inv_CD_tf_epsilon; 
  
  //For a quaternary alloy A-B-C-D
  //Diffusion potential of component D is also a coupled 
  //variable to this kernel
  const VariableValue & _D_diff_pot;
  unsigned int _D_diff_pot_var;
};
