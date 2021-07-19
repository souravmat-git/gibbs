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

class ConstraintSumPF;

template <>
InputParameters validParams< ConstraintSumPF>();

/**
 * Compute the Allen-Cahn interface term with constant Mobility and Interfacial parameter
 */
class  ConstraintSumPF : public Kernel
{
public:
   ConstraintSumPF(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override ;
  virtual Real computeQpJacobian() override ;
  virtual Real computeQpOffDiagJacobian (unsigned int jvar) override ;
  
  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;
  
  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;
  
  const VariableValue & _phase_gamma;
  unsigned int _phase_gamma_var;
 
};
