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

// Forward Declarations
class TwoDimDisplacementGradientConstraintYFun;
class Function;

template <>
InputParameters validParams<TwoDimDisplacementGradientConstraintYFun>();


class TwoDimDisplacementGradientConstraintYFun : public Kernel
{
public:
  TwoDimDisplacementGradientConstraintYFun(const InputParameters & parameters);

protected:  
  // This kernel acts on ex_beta; variable
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  Real _h()  const;
  Real _ay() const;
  
  const Function & _eta;
  //const Function & _ex;
  
  const VariableGradient & _grad_uy;
  unsigned int _uy_var;

  //ex_alpha is the compatible strain
  const VariableValue  & _ey_alpha;
  unsigned int _ey_alpha_var;

  //Transformation strains of each phase
  const MaterialProperty<Real> &  _eyT_alpha;
  const MaterialProperty<Real> &  _eyT_beta;
};
