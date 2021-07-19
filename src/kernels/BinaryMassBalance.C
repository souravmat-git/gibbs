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

#include "BinaryMassBalance.h"

registerMooseObject("gibbsApp", BinaryMassBalance);

template <>
InputParameters
validParams<BinaryMassBalance>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  params.addRequiredCoupledVar("B_diff_pot", "Component B diffusion potential");
  return params;
}

BinaryMassBalance::BinaryMassBalance(const InputParameters & parameters)
  : Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  _grad_B_diff_pot(coupledGradient("B_diff_pot")),
  _B_diff_pot_var(coupled("B_diff_pot")),
  //Phase Material
  _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
  _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
  //_inv_B_td_alpha(getMaterialProperty<Real>("inv_B_td_alpha")),
  //_inv_B_td_beta(getMaterialProperty<Real>("inv_B_td_beta")),
  //interpolation material
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
BinaryMassBalance::thermodynamic_factor() const
{
 //the thermodynamic factor is same as: BinaryPhaseConstraint::computeQpDiagJac
  return (1.0/(_h[_qp] * _inv_B_tf_beta[_qp] + (1.0-_h[_qp]) * _inv_B_tf_alpha[_qp]));
}

Real
BinaryMassBalance::L_BB_interp() const
{
  return (_L_BB_beta[_qp]* _h[_qp] + _L_BB_alpha[_qp]*(1.0-_h[_qp]));
}

Real
BinaryMassBalance::dL_BB_muB_interp() const
{
  return (_dL_BB_muB_beta[_qp]* _h[_qp] + _dL_BB_muB_alpha[_qp]*(1.0-_h[_qp]));
}

Real
BinaryMassBalance::computeQpResidual()
{     //Variable on which the kernel operates: X_B
  return (_grad_test[_i][_qp] * BinaryMassBalance::L_BB_interp() * _grad_B_diff_pot[_qp]);
}


Real
BinaryMassBalance::computeQpJacobian(){  
  //if the derivative of the residual wrt to the comp is taken
  return (_grad_test[_i][_qp] *(BinaryMassBalance::L_BB_interp() * BinaryMassBalance::thermodynamic_factor() * _grad_phi[_j][_qp]));
}


Real
BinaryMassBalance::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _B_diff_pot_var)
  { 
      return (_grad_test[_i][_qp] * ( (BinaryMassBalance::L_BB_interp() * _grad_phi[_j][_qp])
                               +  (BinaryMassBalance::dL_BB_muB_interp()* _grad_B_diff_pot[_qp]) * _phi[_j][_qp]) );
  }
  else if (jvar == _eta_var)
  {
    return (_grad_test[_i][_qp] * _dh[_qp]* (_L_BB_beta[_qp] - _L_BB_alpha[_qp]) * _grad_B_diff_pot[_qp])* _phi[_j][_qp];
  }
  else
  {
    return 0;
  }
}
