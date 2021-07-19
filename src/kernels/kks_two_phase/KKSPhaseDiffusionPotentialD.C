//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseDiffusionPotentialD.h"
registerMooseObject("gibbsApp", KKSPhaseDiffusionPotentialD);

template <>
InputParameters
validParams<KKSPhaseDiffusionPotentialD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("KKS model kernel to enforce the pointwise equality of phase chemical "
                             "potentials dF_{beta}/dc{beta} - dF_{alpha}/dx_{alpha}. The non-linear variable of this "
                             "kernel is xD_{alpha}.");
  params.addRequiredCoupledVar(
      "xD_beta", "Component D in beta phase"); // xD_{beta} is the coupled variable
  params.addRequiredCoupledVar("xC_beta","Component C in beta phase");
  params.addRequiredCoupledVar("xC_alpha","Component C in alpha phase");
  params.addRequiredCoupledVar("xB_beta","Component B in beta phase");
  params.addRequiredCoupledVar("xB_alpha","Component B in alpha phase");
  return params;
}

KKSPhaseDiffusionPotentialD::KKSPhaseDiffusionPotentialD(const InputParameters & parameters)
  : Kernel(parameters),
    _xD_beta(coupledValue("xD_beta")),
    _xD_beta_var(coupled("xD_beta")),
    //Component D is coupled with comp C
    _xC_beta(coupledValue("xC_beta")),
    _xC_beta_var(coupled("xC_beta")),
    _xC_alpha(coupledValue("xC_alpha")),
    _xC_alpha_var(coupled("xC_alpha")),
    //Component D is coupled with comp B
    _xB_beta(coupledValue("xB_beta")),
    _xB_beta_var(coupled("xB_beta")),
    _xB_alpha(coupledValue("xB_alpha")),
    _xB_alpha_var(coupled("xB_alpha")),
    //Material properties required for a quartenary alloy A-B-C-D
    _D_diff_pot_alpha(getMaterialProperty<Real>("D_diff_pot_alpha")),
    _D_diff_pot_beta(getMaterialProperty<Real>("D_diff_pot_beta")),
    _D_therm_factor_alpha(getMaterialProperty<Real>("D_therm_factor_alpha")),
    _D_therm_factor_beta(getMaterialProperty<Real>("D_therm_factor_beta")),
    _CD_therm_factor_beta(getMaterialProperty<Real>("CD_therm_factor_beta")),
    _CD_therm_factor_alpha(getMaterialProperty<Real>("CD_therm_factor_alpha")),
    //Material properties required for quartenary alloy A-B-C-D
    _BD_therm_factor_beta(getMaterialProperty<Real>("BD_therm_factor_beta")),
    _BD_therm_factor_alpha(getMaterialProperty<Real>("BD_therm_factor_alpha"))  
{
}

Real
KKSPhaseDiffusionPotentialD::computeQpResidual()
{
  return _test[_i][_qp] * (_D_diff_pot_beta[_qp]- _D_diff_pot_alpha[_qp]);
}

Real
KKSPhaseDiffusionPotentialD::computeQpJacobian()
{
  // derivative with repect to xD_alpha
  return - _test[_i][_qp] * _phi[_j][_qp] * (_D_therm_factor_alpha[_qp]);
}

Real
KKSPhaseDiffusionPotentialD::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xD_beta_var)
  {
 // derivative with respect to the coupled-variable _xD_beta
    return  _test[_i][_qp] * _phi[_j][_qp] * (_D_therm_factor_beta[_qp]);
  }  
  else if (jvar == _xC_beta_var)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_CD_therm_factor_beta[_qp]);
  }
  else if (jvar == _xC_alpha_var)
  {
    return - (_test[_i][_qp] * _phi[_j][_qp] * (_CD_therm_factor_alpha[_qp]));
  }
  else if (jvar == _xB_beta_var)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * (_BD_therm_factor_beta[_qp]);
  }
  else if (jvar == _xB_alpha_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp] * _BD_therm_factor_alpha[_qp]);
  }
  else
    return 0.0;     
}
