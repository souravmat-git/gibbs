//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This code modifies the KKSPhaseConcentration in MOOSE
//* by including the hand coded interpolation function
//* instead of using the Material property- S.Chatterjee

#include "BinaryPhaseConstraint.h"

registerMooseObject("gibbsApp", BinaryPhaseConstraint);

template <>
InputParameters
validParams<BinaryPhaseConstraint>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (1-h(eta))*xB_alpha + h(eta)*xB_beta - xB = 0."
                             "non-linear variable of this kernel is XB");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  params.addRequiredCoupledVar("B_diff_pot", "Component B diffusion potential");
  return params;
}

BinaryPhaseConstraint::BinaryPhaseConstraint(const InputParameters & parameters)
  :  Kernel(parameters),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _B_diff_pot(coupledValue("B_diff_pot")),
    _B_diff_pot_var(coupled("B_diff_pot")),
    _xB_alpha(getMaterialProperty<Real>("xB_alpha")),
    _xB_beta(getMaterialProperty<Real>("xB_beta")),
    _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
    _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
    //interpolation material
    _h(getMaterialProperty<Real>("h")),
    _dh(getMaterialProperty<Real>("dh"))
{
}

Real
BinaryPhaseConstraint::computeQpResidual()
{
  // Variable that this kernel operates on X_{B}
  return _test[_i][_qp] * (_h[_qp] *_xB_beta[_qp] + (1.0-_h[_qp]) * _xB_alpha[_qp] 
         - _u[_qp]);
}

Real
BinaryPhaseConstraint::computeQpJacobian()
{
 //Derivative with respect to XB
  return -(_test[_i][_qp] * _phi[_j][_qp]) ;
}

Real
BinaryPhaseConstraint::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _eta_var)
  {
    return _test[_i][_qp] * (_xB_beta[_qp] - _xB_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
  }
  else if (jvar == _B_diff_pot_var)
  {
    return (_test[_i][_qp] * (_h[_qp] *  _inv_B_tf_beta[_qp] + (1.0-_h[_qp]) * _inv_B_tf_alpha[_qp])
                         *  _phi[_j][_qp]);
  }
  else
    return 0.0;
}
