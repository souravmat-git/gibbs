//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This was written by S.Chatterjee

#include "MultiCompContinuityEquationD.h"

registerMooseObject("gibbsApp", MultiCompContinuityEquationD);

template <>
InputParameters
validParams<MultiCompContinuityEquationD>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation on mu_{1}"
                              "Eqn:  nabla * (M_{11}(nabla mu_{1}) + nabla * (M_{12}(nabla mu_{2})) = 0");
  params.addRequiredCoupledVar("C_diff_pot", "diffusion potential of component C");
  params.addRequiredCoupledVar("B_diff_pot","Diffusion potential of comp B");
  params.addRequiredCoupledVar("eta", "Phase field variable");
  return params;
}

MultiCompContinuityEquationD::MultiCompContinuityEquationD(const InputParameters & parameters)
  : Kernel(parameters),
    _grad_C_diff_pot(coupledGradient("C_diff_pot")),
    _grad_C_diff_pot_var(coupled("C_diff_pot")),
    _grad_B_diff_pot(coupledGradient("B_diff_pot")),
    _grad_B_diff_pot_var(coupled("B_diff_pot")),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _L_DD_beta(getMaterialProperty<Real>("L_DD_beta")),
    _L_CD_beta(getMaterialProperty<Real>("L_CD_beta")),
    _L_BD_beta(getMaterialProperty<Real>("L_BD_beta")),
    _L_DD_alpha(getMaterialProperty<Real>("L_DD_alpha")),
    _L_CD_alpha(getMaterialProperty<Real>("L_CD_alpha")),
    _L_BD_alpha(getMaterialProperty<Real>("L_BD_alpha")),
    _h(getMaterialProperty<Real> ("h")),
    _dh(getMaterialProperty<Real>("dh"))
{
}

Real
MultiCompContinuityEquationD::computeQpResidual()
{
   _L_DD_interp = _L_DD_beta[_qp] * _h[_qp] +  _L_DD_alpha[_qp] * (1.0 - _h[_qp]);
   _L_CD_interp = _L_CD_beta[_qp] * _h[_qp] +  _L_CD_alpha[_qp] * (1.0 - _h[_qp]);
   _L_BD_interp = _L_BD_beta[_qp] * _h[_qp] +  _L_BD_alpha[_qp] * (1.0 - _h[_qp]);

  return (_grad_test[_i][_qp] *(_L_DD_interp * _grad_u[_qp] 
                              + _L_CD_interp * _grad_C_diff_pot[_qp]
                              + _L_BD_interp * _grad_B_diff_pot[_qp]));
}

Real
MultiCompContinuityEquationD::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * (_L_DD_beta[_qp] * _h[_qp]
                                +_L_DD_alpha[_qp] * (1.0 - _h[_qp]))* _grad_phi[_j][_qp]);
}

Real
MultiCompContinuityEquationD::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _grad_C_diff_pot_var)
  {
    return (_grad_test[_i][_qp] * (_L_CD_beta[_qp] * _h[_qp]
                                 + _L_CD_alpha[_qp] * (1.0 - _h[_qp]))* _grad_phi[_j][_qp]);
  }
  else if (jvar == _grad_B_diff_pot_var)
  {
    return (_grad_test[_i][_qp] * (_L_BD_beta[_qp] * _h[_qp]
                                 + _L_BD_alpha[_qp] * (1.0 - _h[_qp]))* _grad_phi[_j][_qp]);
  }
  else if (jvar == _eta_var)
  {
    return (_grad_test[_i][_qp] * _dh[_qp] * ((_L_DD_beta[_qp] - _L_DD_alpha[_qp]) * _grad_u[_qp]
                                             +(_L_CD_beta[_qp] - _L_CD_alpha[_qp]) * _grad_C_diff_pot[_qp]
                                             +(_L_BD_beta[_qp] - _L_BD_alpha[_qp]) * _grad_B_diff_pot[_qp])* _phi[_j][_qp]);
  }
  else
  {
    return 0.0; 
  }
}
