//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BinaryPhaseConstraintMu.h"

registerMooseObject("gibbsApp", BinaryPhaseConstraintMu);

template <>
InputParameters
validParams<BinaryPhaseConstraintMu>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (1-h(eta))*xB_alpha + h(eta)*xB_beta - xB = 0."
                             "non-linear variable of this kernel is XB");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  params.addRequiredCoupledVar("xB", "Component B mole fraction");
  return params;
}

BinaryPhaseConstraintMu::BinaryPhaseConstraintMu(const InputParameters & parameters)
  :  Kernel(parameters),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _xB(coupledValue("xB")),
    _xB_var(coupled("xB")),
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
BinaryPhaseConstraintMu::computeQpResidual()
{
  // Variable on which the kernel operates: B_diff_pot
  return _test[_i][_qp] * (_h[_qp] *_xB_beta[_qp] + (1.0-_h[_qp]) * _xB_alpha[_qp] 
         - _xB[_qp]);
}

Real
BinaryPhaseConstraintMu::computeQpJacobian()
{
 //Derivative with respect to B_diff_pot
  return (_test[_i][_qp] * (_h[_qp] *  _inv_B_tf_beta[_qp] + (1.0-_h[_qp]) * _inv_B_tf_alpha[_qp])
                         *  _phi[_j][_qp]);

}

Real
BinaryPhaseConstraintMu::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _eta_var)
  {
    return _test[_i][_qp] * (_xB_beta[_qp] - _xB_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
  }
  else if (jvar == _xB_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp]) ;
  }
  else
    return 0.0;
}
