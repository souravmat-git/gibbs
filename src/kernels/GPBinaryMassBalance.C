//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "GPBinaryMassBalance.h"

registerMooseObject("gibbsApp", GPBinaryMassBalance);

template <>
InputParameters
validParams<GPBinaryMassBalance>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  return params;
}

GPBinaryMassBalance::GPBinaryMassBalance(const InputParameters & parameters)
  :Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  //Interpolation Material
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  //Kinetic material
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _dL_BB_muB_beta(getMaterialProperty<Real>("dL_BB_muB_beta")),
  _dL_BB_muB_alpha(getMaterialProperty<Real>("dL_BB_muB_alpha"))
{
}

Real
GPBinaryMassBalance::L_BB_interp() const
{
  return (_L_BB_beta[_qp]* _h[_qp] + _L_BB_alpha[_qp]*(1.0-_h[_qp]));
}

Real
GPBinaryMassBalance::dL_BB_muB_interp() const
{
  return (_dL_BB_muB_beta[_qp]* _h[_qp] + _dL_BB_muB_alpha[_qp]*(1.0-_h[_qp]));
}

Real
GPBinaryMassBalance::computeQpResidual()
{   
  //Variable on which the kernel operates: mu_B
  return (_grad_test[_i][_qp] * GPBinaryMassBalance::L_BB_interp() * _grad_u[_qp]);
}


Real
GPBinaryMassBalance::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * ( (GPBinaryMassBalance::L_BB_interp() * _grad_phi[_j][_qp])
                               +  (GPBinaryMassBalance::dL_BB_muB_interp() * _grad_u[_qp]) * _phi[_j][_qp]) );
}


Real
GPBinaryMassBalance::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _eta_var)
  {
    return (_grad_test[_i][_qp] * _dh[_qp]* (_L_BB_beta[_qp] - _L_BB_alpha[_qp]) * _grad_u[_qp]* _phi[_j][_qp] );
  }
  else
  {
    return 0;
  }
}
