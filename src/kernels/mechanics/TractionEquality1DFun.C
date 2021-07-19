//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TractionEquality1DFun.h"
#include "Function.h"
registerMooseObject("gibbsApp", TractionEquality1DFun);

template <>
InputParameters
validParams<TractionEquality1DFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Enforces equality of traction along x-direction");
  params.addRequiredParam<FunctionName>("eta", "Phase-field variable that distinguish phases");
  params.addRequiredCoupledVar("ux", "Displacement field in x-direction");
  return params;
}

TractionEquality1DFun::TractionEquality1DFun(const InputParameters & parameters)
  : Kernel(parameters),
   _eta(getFunction("eta")),
   _ux_var(coupled("ux")),
   _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
   _sx_beta(getMaterialProperty<Real>("sx_beta")),
   _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
   _mat_const_beta(getMaterialProperty<Real>("mat_const_beta")),
   _h(getMaterialProperty<Real>("h"))
{
}

Real
TractionEquality1DFun::rev_interp_mat_const() const {
  return (_mat_const_beta[_qp]  * (1.0- _h[_qp])
         +_mat_const_alpha[_qp] * _h[_qp]);
}         

Real
TractionEquality1DFun::computeQpResidual()
{//This kernel acts on a, i,e, the magnitude of strain jump
  return _test[_i][_qp] * (_sx_beta[_qp] - _sx_alpha[_qp]);
}

Real
TractionEquality1DFun::computeQpJacobian()
{                                
  return _test[_i][_qp] * rev_interp_mat_const() * _phi[_j][_qp];
}

Real
TractionEquality1DFun::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ux_var)
  {
    return _test[_i][_qp] * (_mat_const_beta[_qp] - _mat_const_alpha[_qp]) * _grad_phi[_j][_qp](0);
  }
  else
    return 0;
 }  
