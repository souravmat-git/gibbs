//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryPhaseDiffusionPotentialKKSFun.h"
#include "Function.h"

registerMooseObject("gibbsApp", BinaryPhaseDiffusionPotentialKKSFun);

template <>
InputParameters
validParams<BinaryPhaseDiffusionPotentialKKSFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("KKS model kernel to enforce the pointwise equality of phase chemical "
                             "potentials dF_{beta}/dc{beta} - dF_{alpha}/dx_{alpha}. The non-linear variable of this "
                             "kernel is xB_{alpha}.");
  params.addRequiredCoupledVar(
      "xB_beta", "Component B in beta phase"); // xB_{beta} is the coupled variable
  return params;
}

BinaryPhaseDiffusionPotentialKKSFun::BinaryPhaseDiffusionPotentialKKSFun(const InputParameters & parameters)
  : Kernel(parameters),
    _xB_beta(coupledValue("xB_beta")),
    _xB_beta_var(coupled("xB_beta")),
   //Material properties required for binary alloy A-B
    _B_diff_pot_alpha(getMaterialProperty<Real>("B_diff_pot_alpha")),
    _B_diff_pot_beta(getMaterialProperty<Real>("B_diff_pot_beta")),
    _B_therm_factor_alpha(getMaterialProperty<Real>("B_therm_factor_alpha")),
    _B_therm_factor_beta(getMaterialProperty<Real>("B_therm_factor_beta")) 
{
}

Real
BinaryPhaseDiffusionPotentialKKSFun::computeQpResidual()
{
  return (_test[_i][_qp] * (_B_diff_pot_beta[_qp] - _B_diff_pot_alpha[_qp]));
}

Real
BinaryPhaseDiffusionPotentialKKSFun::computeQpJacobian()
{
   // derivative with repect to xB_alpha
  return -(_test[_i][_qp] * _phi[_j][_qp] * _B_therm_factor_alpha[_qp]);
}

Real
BinaryPhaseDiffusionPotentialKKSFun::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xB_beta_var)
  {
    // derivative with respect to the coupled-variable xB_beta
    return (_test[_i][_qp] * _phi[_j][_qp] * _B_therm_factor_beta[_qp]);
  }   
  else
    return 0.0;     
}
