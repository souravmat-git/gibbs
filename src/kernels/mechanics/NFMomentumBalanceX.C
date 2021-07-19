//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* This was written by S.Chatterjee

#include "NFMomentumBalanceX.h"
registerMooseObject("gibbsApp", NFMomentumBalanceX);

template <>
InputParameters
validParams<NFMomentumBalanceX>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  params.addRequiredCoupledVar("ex_beta", "Strain in the beta phase");
  return params;
}

NFMomentumBalanceX::NFMomentumBalanceX(const InputParameters & parameters)
  :Kernel(parameters),
  _ex_beta(coupledValue("ex_beta")),
  _ex_beta_var(coupled("ex_beta")),
  //Fetch the material properties
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _mat_const_beta(getMaterialProperty<Real>("mat_const_beta"))
{
}

Real
NFMomentumBalanceX::computeQpResidual(){
  //Non-linear variable that this equation is coupled to is 
  //the displacement field variable
    return -_grad_test[_i][_qp](0) * _sx_beta[_qp];          
}

Real
NFMomentumBalanceX::computeQpJacobian(){
 //derivative wrt ux is zero.
  return 0;
}

Real
NFMomentumBalanceX::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _ex_beta_var){
    return (-_grad_test[_i][_qp](0)  * _mat_const_beta[_qp] * _phi[_j][_qp]);
  }
  else 
    return 0;
}           
