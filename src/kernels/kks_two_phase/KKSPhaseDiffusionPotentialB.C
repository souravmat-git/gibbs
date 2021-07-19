//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseDiffusionPotentialB.h"
registerMooseObject("gibbsApp", KKSPhaseDiffusionPotentialB);

template <>
InputParameters
validParams<KKSPhaseDiffusionPotentialB>()
{
  InputParameters params = validParams<BinaryPhaseDiffusionPotentialKKS>();
  params.addClassDescription("KKS model kernel to enforce the pointwise equality of phase chemical "
                             "potentials dF_{beta}/dc{beta} - dF_{alpha}/dx_{alpha}. The non-linear variable of this "
                             "kernel is xB_{alpha}.");
  params.addRequiredCoupledVar("xC_beta", "Component C in beta phase");
  params.addRequiredCoupledVar("xC_alpha","Component C in alpha phase");
  params.addCoupledVar("xD_beta",0.0, "Component D in beta phase");
  params.addCoupledVar("xD_alpha", 0.0,"Component D in alpha phase");
  params.addRequiredParam<MaterialPropertyName>("BC_therm_factor_beta","mu(B).x(C) - mu(A).x(C)");
  params.addRequiredParam<MaterialPropertyName>("BC_therm_factor_alpha","mu(B).x(C) -mu(A).x(C)");
  params.addParam<MaterialPropertyName>("BD_therm_factor_beta", 0.0,"mu(B).x(D) - mu(B).x(D)");
  params.addParam<MaterialPropertyName>("BD_therm_factor_alpha",0.0,"mu(B).x(D) - mu(B).x(D)");
  return params;
}

KKSPhaseDiffusionPotentialB::KKSPhaseDiffusionPotentialB(const InputParameters & parameters)
  : BinaryPhaseDiffusionPotentialKKS(parameters),
    //Component C is coupled with comp B
    _xC_beta(coupledValue("xC_beta")),
    _xC_beta_var(coupled("xC_beta")),
    _xC_alpha(coupledValue("xC_alpha")),
    _xC_alpha_var(coupled("xC_alpha")),
    //Component D is coupled with comp B
    _xD_beta(coupledValue("xD_beta")),
    _xD_beta_var(coupled("xD_beta")),
    _xD_alpha(coupledValue("xD_alpha")),
    _xD_alpha_var(coupled("xD_alpha")),
    //Material properties required for ternary alloy A-B-C
    _BC_therm_factor_beta(getMaterialProperty<Real>("BC_therm_factor_beta")),
    _BC_therm_factor_alpha(getMaterialProperty<Real>("BC_therm_factor_alpha")),
    //Material properties required for quartenary alloy A-B-C-D
    _BD_therm_factor_beta(getMaterialProperty<Real>("BD_therm_factor_beta")),
    _BD_therm_factor_alpha(getMaterialProperty<Real>("BD_therm_factor_alpha"))  
{
}

Real
KKSPhaseDiffusionPotentialB::computeQpResidual()
{
  return BinaryPhaseDiffusionPotentialKKS::computeQpResidual();
}

Real
KKSPhaseDiffusionPotentialB::computeQpJacobian()
{
  return BinaryPhaseDiffusionPotentialKKS::computeQpJacobian();
}

Real
KKSPhaseDiffusionPotentialB::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xB_beta_var)
  {
    return BinaryPhaseDiffusionPotentialKKS::computeQpOffDiagJacobian(jvar);
  }  
  else if (jvar == _xC_beta_var)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_BC_therm_factor_beta[_qp]);
  }
  else if (jvar == _xC_alpha_var)
  {
    return - (_test[_i][_qp] * _phi[_j][_qp] * (_BC_therm_factor_alpha[_qp]));
  }
  else if (jvar == _xD_beta_var)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_BD_therm_factor_beta[_qp]);
  }
  else if (jvar == _xD_alpha_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp] * _BD_therm_factor_alpha[_qp]);
  }
  else
    return 0.0;     
}
