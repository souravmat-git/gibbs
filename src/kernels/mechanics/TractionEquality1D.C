//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TractionEquality1D.h"
registerMooseObject("gibbsApp", TractionEquality1D);

template <>
InputParameters
validParams<TractionEquality1D>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Enforces equality of traction along x-direction");
  params.addRequiredCoupledVar("eta", "Phase-field variable that distinguish phases");
  params.addRequiredCoupledVar("ux", "Displacement field in x-direction");
  return params;
}

TractionEquality1D::TractionEquality1D(const InputParameters & parameters)
  : Kernel(parameters),
   _eta_var(coupled("eta")),
   _ux_var(coupled("ux")),
   _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
   _sx_beta(getMaterialProperty<Real>("sx_beta")),
   _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
   _mat_const_beta(getMaterialProperty<Real>("mat_const_beta")),
   _h(getMaterialProperty<Real>("h")),
   _dh(getMaterialProperty<Real>("dh"))
{
}

Real
TractionEquality1D::rev_interp_mat_const() const {
  return (_mat_const_beta[_qp]  * (1.0- _h[_qp])
         +_mat_const_alpha[_qp] * _h[_qp]);
}         

Real
TractionEquality1D::computeQpResidual()
{  //This kernel acts on a, i,e, the magnitude of strain jump
  return _test[_i][_qp] * (_sx_beta[_qp] - _sx_alpha[_qp]);
}

Real
TractionEquality1D::computeQpJacobian()
{                                
  return (_test[_i][_qp] * rev_interp_mat_const() * _phi[_j][_qp]);
}

Real
TractionEquality1D::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _eta_var)
  {
    return -_test[_i][_qp] *_u[_qp] * _dh[_qp]* (_mat_const_beta[_qp] - _mat_const_alpha[_qp]) * _phi[_j][_qp];     
    
  }
  else if (jvar == _ux_var)
  {
    return _test[_i][_qp] * (_mat_const_beta[_qp] - _mat_const_alpha[_qp]) * _grad_phi[_j][_qp](0);
  }
  else
    return 0;
 }  
