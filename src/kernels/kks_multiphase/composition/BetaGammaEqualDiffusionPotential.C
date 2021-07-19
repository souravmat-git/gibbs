//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BetaGammaEqualDiffusionPotential.h"
registerMooseObject("gibbsApp", BetaGammaEqualDiffusionPotential);

template <>
InputParameters
validParams<BetaGammaEqualDiffusionPotential>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("KKS model kernel to enforce the pointwise equality of phase chemical "
                             "potentials dF_{beta}/dc{beta} - dF_{gamma}/dc_{gamma}. The non-linear variable of this "
                             "kernel is c_{beta}.");
  params.addRequiredCoupledVar(
      "phase_comp_gamma", "Phase composition in gamma phase"); // c_{alpha} is the coupled variable
  params.addRequiredParam<MaterialPropertyName>("df_dcgamma","diffusion potential of gamma phase");
  params.addRequiredParam<MaterialPropertyName>("df_dcbeta","diffusion potential of beta phase" );
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2gamma", "diffusion potential of gamma phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2beta", "Equilibrium composition in the beta phase");
  return params;
}

BetaGammaEqualDiffusionPotential::BetaGammaEqualDiffusionPotential(const InputParameters & parameters)
  : Kernel(parameters),
    _phase_comp_gamma(coupledValue("phase_comp_gamma")),
    _phase_comp_gamma_var(coupled("phase_comp_gamma")),
    _df_gamma(getMaterialProperty<Real>("df_dcgamma")),
    _df_beta(getMaterialProperty<Real>("df_dcbeta")),
    _d2f_gamma(getMaterialProperty<Real>("d2f_dc2gamma")),
    _d2f_beta(getMaterialProperty<Real>("d2f_dc2beta"))
{
}

Real
BetaGammaEqualDiffusionPotential::computeQpResidual()
{
  return _test[_i][_qp] * (_df_beta[_qp]- _df_gamma[_qp]);
}

Real
BetaGammaEqualDiffusionPotential::computeQpJacobian()
{
  // derivative with repect to the variable that this kernel operates on
  // for this case with respect to c_{beta}
  return  _test[_i][_qp] * _phi[_j][_qp] * (_d2f_beta[_qp]);
}

Real
BetaGammaEqualDiffusionPotential::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_comp_gamma_var)
  {
 // derivative with respect to the coupled-variable
    return - _test[_i][_qp] * _phi[_j][_qp] * (_d2f_gamma[_qp]);
  }
  else
    return 0.0;     
}
