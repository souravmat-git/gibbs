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
#include "Kernel.h"

// Forward Declarations
class TwoPhaseStrainConstraint;
template <>
InputParameters validParams<TwoPhaseStrainConstraint>();
/**
 * Enforce sum of phase strains to the total strain (e_x)
 * The non-linear variable that this kernel acts on
 * is the stress in the x-direction
 * This kernel will be coupled to stress tensor, eta function
 **/
class TwoPhaseStrainConstraint : public Kernel
{
public:
 TwoPhaseStrainConstraint(const InputParameters & parameters);

protected:  
  // The override command ensures that a virtual member function
  // from a base class is overridden  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  const VariableValue & _eta;
  unsigned int _eta_var;

  //Overall displacement gradient is a coupled variable
  const VariableGradient & _grad_ux;
  unsigned int _ux_var;
  
  //Strain_x is a funtion of stress through the constitutive equation
  const MaterialProperty<Real> & _ex_alpha;
  const MaterialProperty<Real> & _ex_beta;
  
  //First derivative of the phase strain with respect to stress in alpha
  const MaterialProperty<Real> & _compliance_alpha;
  const MaterialProperty<Real> & _compliance_beta;
  
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh; 
 
};
