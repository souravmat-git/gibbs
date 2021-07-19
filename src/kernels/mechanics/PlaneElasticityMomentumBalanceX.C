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

#include "PlaneElasticityMomentumBalanceX.h"

registerMooseObject("gibbsApp", PlaneElasticityMomentumBalanceX);

template <>
InputParameters
validParams<PlaneElasticityMomentumBalanceX>()
{
  InputParameters params = validParams<PlaneElasticityBase>();
  params.addRequiredCoupledVar("disp_y","Displacement in the y-direction");
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  params.addRequiredParam<MaterialPropertyName>("force_x", "Body force in x-dir");

  return params;
}

PlaneElasticityMomentumBalanceX::PlaneElasticityMomentumBalanceX(const InputParameters & parameters)
  :PlaneElasticityBase(parameters),
  _grad_uy(coupledGradient("disp_y")),
  _uy_var(coupled("disp_y")),
  _fx(getMaterialProperty<Real>(getParam<MaterialPropertyName>("force_x")))
{
}

Real
PlaneElasticityMomentumBalanceX::computeQpResidual()
{ 
  //The non-linear variable that this kernel acts on is the displacement in the x-direction;
  //The variable is coupled to the displacement in the y-direction
  //_grad_u[_qp](0) - x-component of the gradient of the displacement in x-direction
  //_grad_u[_qp](1) - y-component of the gradient of the displacement in x-direction
  
  return _grad_test[_i][_qp](0) * (_C11[_qp] * _grad_u[_qp](0) + _C12[_qp] * _grad_uy[_qp](1)) 
       + _grad_test[_i][_qp](1) * (_C66[_qp] * _grad_u[_qp](1) + _C66[_qp] * _grad_uy[_qp](0))
       - _test[_i][_qp]*_fx[_qp] ;
}

Real
PlaneElasticityMomentumBalanceX::computeQpJacobian()
{
  //Derivative with respect to displacement in the x-direction
  return _grad_test[_i][_qp](0) * _C11[_qp] * _grad_phi[_j][_qp](0) +
         _grad_test[_i][_qp](1) * _C66[_qp] * _grad_phi[_j][_qp](1);
}

Real
PlaneElasticityMomentumBalanceX::computeQpOffDiagJacobian(unsigned int jvar)
{ 
  if (jvar == _uy_var){
    return _grad_test[_i][_qp](0)* _C12[_qp] * _grad_phi[_j][_qp](1) +
           _grad_test[_i][_qp](1)* _C66[_qp] * _grad_phi[_j][_qp](0);
 }          
 else
   return 0;
}
