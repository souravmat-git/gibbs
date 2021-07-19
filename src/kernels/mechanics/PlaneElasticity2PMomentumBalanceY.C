//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* This was written by S.Chatterjee

#include "PlaneElasticity2PMomentumBalanceY.h"
registerMooseObject("gibbsApp", PlaneElasticity2PMomentumBalanceY);

template <>
InputParameters
validParams<PlaneElasticity2PMomentumBalanceY>()
{
  InputParameters params = validParams<PlaneElasticityTwoPhaseBase>();
  params.addRequiredCoupledVar("disp_x","Displacement in the x-direction");
  params.addClassDescription("This kernel implements the momentum balance in y-dir");
  return params;
}

PlaneElasticity2PMomentumBalanceY::PlaneElasticity2PMomentumBalanceY(const InputParameters & parameters)
  :PlaneElasticityTwoPhaseBase(parameters),
  //_grad_ux(coupledGradient("disp_x")),
  _ux_var(coupled("disp_x"))
{
}
        

Real
PlaneElasticity2PMomentumBalanceY::computeQpResidual(){ 
  //The non-linear variable that this kernel acts on is the displacement in the x-direction;
  //The variable is coupled to the displacement in the y-direction
  //_grad_u[_qp](0) - x-component of the gradient of the displacement in x-direction
  //_grad_u[_qp](1) - y-component of the gradient of the displacement in x-direction
  
  Real _interpolated_sy  = (_h[_qp] * _sy_beta[_qp]  + (1.0-_h[_qp]) * _sy_alpha[_qp]);   
  Real _interpolated_sxy = (_h[_qp] * _sxy_beta[_qp] + (1.0-_h[_qp]) * _sxy_alpha[_qp]);
  
  return - _grad_test[_i][_qp](0) * _interpolated_sxy 
         - _grad_test[_i][_qp](1) * _interpolated_sy;
}

Real
PlaneElasticity2PMomentumBalanceY::computeQpJacobian(){
  //Derivative with respect to displacement in the y-direction
  return - _grad_test[_i][_qp](0) * PlaneElasticityTwoPhaseBase::C66_interp() * _grad_phi[_j][_qp](0)
         - _grad_test[_i][_qp](1) * PlaneElasticityTwoPhaseBase::C22_interp() * _grad_phi[_j][_qp](1);
}

Real
PlaneElasticity2PMomentumBalanceY::computeQpOffDiagJacobian(unsigned int jvar){ 
  if (jvar == _ux_var){
    //Derivative wrt to displacement in the x-direction
    return -_grad_test[_i][_qp](0)* PlaneElasticityTwoPhaseBase::C66_interp() * _grad_phi[_j][_qp](1)
           -_grad_test[_i][_qp](1)* PlaneElasticityTwoPhaseBase::C12_interp() * _grad_phi[_j][_qp](0);
 }
 else if (jvar == _eta_var){ 
  return -_grad_test[_i][_qp](0) * (_sxy_beta[_qp] - _sxy_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp]
         -_grad_test[_i][_qp](1) * (_sy_beta[_qp]  - _sy_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
 }          
 else 
   return 0;
}
