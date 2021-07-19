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

class LagrangeMultiplierPF;

template <>
InputParameters validParams<LagrangeMultiplierPF>();

/**
 * Compute the Allen-Cahn interface term with constant Mobility and Interfacial parameter
 */
class LagrangeMultiplierPF : public Kernel
{
public:
  LagrangeMultiplierPF(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override ;
  virtual Real computeQpJacobian() override ;
  virtual Real computeQpOffDiagJacobian (unsigned int jvar) override ;
  
  const VariableValue & _lambda;
  unsigned int _lambda_var;

  const MaterialProperty<Real> & _L;
};
