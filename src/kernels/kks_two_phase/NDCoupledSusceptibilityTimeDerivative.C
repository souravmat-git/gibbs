//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NDCoupledSusceptibilityTimeDerivative.h"

registerMooseObject("gibbsApp", NDCoupledSusceptibilityTimeDerivative);

template <>
InputParameters
validParams<NDCoupledSusceptibilityTimeDerivative>()
{
  InputParameters params = validParams<CoupledSusceptibilityTimeDerivative>();
  params.addClassDescription("A modified coupled time derivative Kernel that multiplies the time "
                             "derivative of a coupled variable by a generalized susceptibility");
  params.addCoupledVar("args", "Vector of arguments of the susceptibility");
  return params;
}

NDCoupledSusceptibilityTimeDerivative::NDCoupledSusceptibilityTimeDerivative(const InputParameters & parameters)
  :CoupledSusceptibilityTimeDerivative(parameters),
  _xB_alpha(getMaterialProperty<Real>("xB_alpha")),
  _xB_beta(getMaterialProperty<Real>("xB_beta")),
  _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
  _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
  //interpolation material
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h"))
{
}

Real
NDCoupledSusceptibilityTimeDerivative::computeQpResidual()
{
  return CoupledTimeDerivative::computeQpResidual() * (_xB_beta[_qp] - _xB_alpha[_qp])*_dh[_qp];
}

Real
NDCoupledSusceptibilityTimeDerivative::computeQpJacobian()
{
  return CoupledTimeDerivative::computeQpResidual() * (_inv_B_tf_beta[_qp] - _inv_B_tf_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
}

Real
NDCoupledSusceptibilityTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
 // get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);

  if (jvar == _v_var)
    return (CoupledTimeDerivative::computeQpOffDiagJacobian(jvar) * (_xB_beta[_qp] - _xB_alpha[_qp])*_dh[_qp] +
           CoupledTimeDerivative::computeQpResidual() * _phi[_j][_qp] * (_xB_beta[_qp] - _xB_alpha[_qp])*_d2h[_qp] );

  return CoupledTimeDerivative::computeQpResidual() * _phi[_j][_qp] * (*_dFdarg[cvar])[_qp];
}
