//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AlphaGammaEqualDiffusionPotential.h"
registerMooseObject("gibbsApp", AlphaGammaEqualDiffusionPotential);

template <>
InputParameters
validParams<AlphaGammaEqualDiffusionPotential>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("KKS model kernel to enforce the pointwise equality of phase chemical "
                             "potentials dF_{alpha}/dc{alpha} - dF_{gamma}/dc_{gamma}. The non-linear variable of this "
                             "kernel is c_{alpha}.");
  params.addRequiredCoupledVar(
      "phase_comp_gamma", "Phase composition in gamma phase"); // c_{alpha} is the coupled variable
  params.addRequiredParam<MaterialPropertyName>("df_dcgamma","diffusion potential of gamma phase");
  params.addRequiredParam<MaterialPropertyName>("df_dcalpha","diffusion potential of beta phase" );
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2gamma", "diffusion potential of gamma phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2alpha", "Equilibrium composition in the beta phase");
  return params;
}

AlphaGammaEqualDiffusionPotential::AlphaGammaEqualDiffusionPotential(const InputParameters & parameters)
  : Kernel(parameters),
    _phase_comp_gamma(coupledValue("phase_comp_gamma")),
    _phase_comp_gamma_var(coupled("phase_comp_gamma")),
    _df_gamma(getMaterialProperty<Real>("df_dcgamma")),
    _df_alpha(getMaterialProperty<Real>("df_dcalpha")),
    _d2f_gamma(getMaterialProperty<Real>("d2f_dc2gamma")),
    _d2f_alpha(getMaterialProperty<Real>("d2f_dc2alpha"))
{
}

Real
AlphaGammaEqualDiffusionPotential::computeQpResidual()
{
  return _test[_i][_qp] * (_df_alpha[_qp]- _df_gamma[_qp]);
}

Real
AlphaGammaEqualDiffusionPotential::computeQpJacobian()
{
  // derivative with repect to the variable that this kernel operates on
  // for this case with respect to c_{beta}
  return  _test[_i][_qp] * _phi[_j][_qp] * (_d2f_alpha[_qp]);
}

Real
AlphaGammaEqualDiffusionPotential::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_comp_gamma_var)
  {
 // derivative with respect to the coupled-variable
    return - _test[_i][_qp] * _phi[_j][_qp] * (_d2f_gamma[_qp]);
  }
  else
    return 0.0;     
}
