//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseDiffusionPotentialC.h"
registerMooseObject("gibbsApp", KKSPhaseDiffusionPotentialC);

template <>
InputParameters
validParams<KKSPhaseDiffusionPotentialC>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("KKS model kernel to enforce the pointwise equality of phase chemical "
                             "potentials dF_{beta}/dc{beta} - dF_{alpha}/dx_{alpha}. The non-linear variable of this "
                             "kernel is xC_{alpha}.");
  params.addRequiredCoupledVar(
      "xC_beta", "Component C in beta phase"); // xC_{beta} is the coupled variable
  params.addRequiredCoupledVar("xB_beta", "Component B in beta phase");
  params.addRequiredCoupledVar("xB_alpha", "Component B in alpha phase");
  params.addCoupledVar("xD_beta", 0.0, "Component D in beta phase");
  params.addCoupledVar("xD_alpha",0.0,"Component D in alpha phase");
  params.addParam<MaterialPropertyName>("CD_therm_factor_beta", 0.0, "mu(C).x(D)_beta");
  params.addParam<MaterialPropertyName>("CD_therm_factor_alpha", 0.0, "mu(C).x(D)_alpha");
  return params;
}

KKSPhaseDiffusionPotentialC::KKSPhaseDiffusionPotentialC(const InputParameters & parameters)
  : Kernel(parameters),
    _xC_beta(coupledValue("xC_beta")),
    _xC_beta_var(coupled("xC_beta")),
    //Coupled to Component B
    _xB_beta(coupledValue("xB_beta")),
    _xB_beta_var(coupled("xB_beta")),
    _xB_alpha(coupledValue("xB_alpha")),
    _xB_alpha_var(coupled("xB_alpha")),
    //Coupled to Component D
    _xD_beta(coupledValue("xD_beta")),
    _xD_beta_var(coupled("xD_beta")),
    _xD_alpha(coupledValue("xD_alpha")),
    _xD_alpha_var(coupled("xD_alpha")),
    //Material properties required for ternary alloys A-B-C
    _C_diff_pot_alpha(getMaterialProperty<Real>("C_diff_pot_alpha")),
    _C_diff_pot_beta(getMaterialProperty<Real>("C_diff_pot_beta")),
    _C_therm_factor_alpha(getMaterialProperty<Real>("C_therm_factor_alpha")),
    _C_therm_factor_beta(getMaterialProperty<Real>("C_therm_factor_beta")),
    _BC_therm_factor_beta(getMaterialProperty<Real>("BC_therm_factor_beta")),
    _BC_therm_factor_alpha(getMaterialProperty<Real>("BC_therm_factor_alpha")),
    //Material properties required for quaternary alloys A-B-C-D
    _CD_therm_factor_beta(getMaterialProperty<Real>("CD_therm_factor_beta")),
    _CD_therm_factor_alpha(getMaterialProperty<Real>("CD_therm_factor_alpha"))
{
}

Real
KKSPhaseDiffusionPotentialC::computeQpResidual()
{
  return _test[_i][_qp] * (_C_diff_pot_beta[_qp]- _C_diff_pot_alpha[_qp]);
}

Real
KKSPhaseDiffusionPotentialC::computeQpJacobian()
{
  // derivative with repect to xB_alpha
  return - _test[_i][_qp] * _phi[_j][_qp] * (_C_therm_factor_alpha[_qp]);
}

Real
KKSPhaseDiffusionPotentialC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xC_beta_var)
  {
 // derivative with respect to the coupled-variable xB_beta
    return  _test[_i][_qp] * _phi[_j][_qp] * (_C_therm_factor_beta[_qp]);
  }
  else if (jvar == _xB_beta_var)
  {
     //Differentiate \tilde{mu}_{C}^{\beta} - \tilde{mu}_{C}^{\alpha} = 0
    //w.r.t x_B^{\beta}
    return (_test[_i][_qp] * _phi[_j][_qp] * _BC_therm_factor_beta[_qp]);
  }
  else if (jvar == _xB_alpha_var)
  {
    //Differentiate \tilde{mu}_{C}^{\beta} - \tilde{mu}_{C}^{\alpha} = 0
    //w.r.t x_B^{\alpha}
    return -(_test[_i][_qp] * _phi[_j][_qp] * _BC_therm_factor_alpha[_qp]);
  } 
  else if (jvar == _xD_beta_var)
  {
    return (_test[_i][_qp] * _phi[_j][_qp] * _CD_therm_factor_beta[_qp]);  
  }
  else if (jvar == _xD_alpha_var)
  {
    return (_test[_i][_qp] * _phi[_j][_qp] * _CD_therm_factor_alpha[_qp]);
  }
  else
    return 0.0;     
}
