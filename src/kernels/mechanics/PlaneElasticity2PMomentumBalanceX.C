//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* This was written by S.Chatterjee

#include "PlaneElasticity2PMomentumBalanceX.h"
registerMooseObject("gibbsApp", PlaneElasticity2PMomentumBalanceX);

template <>
InputParameters
validParams<PlaneElasticity2PMomentumBalanceX>()
{
  InputParameters params = validParams<PlaneElasticityTwoPhaseBase>();
  params.addRequiredCoupledVar("disp_y","Displacement in the y-direction");
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  return params;
}

PlaneElasticity2PMomentumBalanceX::PlaneElasticity2PMomentumBalanceX(const InputParameters & parameters)
  :PlaneElasticityTwoPhaseBase(parameters),
  //_grad_uy(coupledGradient("disp_y")),
  _uy_var(coupled("disp_y"))
{
}
         
Real
PlaneElasticity2PMomentumBalanceX::computeQpResidual(){ 
  //The non-linear variable that this kernel acts on is the displacement in the x-direction;
  //The variable is coupled to the displacement in the y-direction
  //_grad_u[_qp](0) - x-component of the gradient of the displacement in x-direction
  //_grad_u[_qp](1) - y-component of the gradient of the displacement in x-direction
  
  Real _interpolated_sx  = (_h[_qp]* _sx_beta[_qp]  + (1.0-_h[_qp]) * _sx_alpha[_qp]);   
  Real _interpolated_sxy = (_h[_qp]* _sxy_beta[_qp] + (1.0-_h[_qp]) * _sxy_alpha[_qp]);
  
  return (- _grad_test[_i][_qp](0) * _interpolated_sx
          - _grad_test[_i][_qp](1) * _interpolated_sxy);                         ;
}

Real
PlaneElasticity2PMomentumBalanceX::computeQpJacobian(){
  //Derivative with respect to displacement in the x-direction
  return -_grad_test[_i][_qp](0) * PlaneElasticityTwoPhaseBase::C11_interp() * _grad_phi[_j][_qp](0)
         -_grad_test[_i][_qp](1) * PlaneElasticityTwoPhaseBase::C66_interp() * _grad_phi[_j][_qp](1);
}

Real
PlaneElasticity2PMomentumBalanceX::computeQpOffDiagJacobian(unsigned int jvar){ 
  if (jvar == _uy_var){
    //Derivative wrt to displacement in the y-direction
    return -_grad_test[_i][_qp](0)* PlaneElasticityTwoPhaseBase::C12_interp() * _grad_phi[_j][_qp](1)
           -_grad_test[_i][_qp](1)* PlaneElasticityTwoPhaseBase::C66_interp() * _grad_phi[_j][_qp](0);
 }
 else if (jvar == _eta_var){ 
  return -_grad_test[_i][_qp](0) * (_sx_beta[_qp]  - _sx_alpha[_qp])  * _dh[_qp] * _phi[_j][_qp]
         -_grad_test[_i][_qp](1) * (_sxy_beta[_qp] - _sxy_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
 }          
 else 
   return 0;
}
