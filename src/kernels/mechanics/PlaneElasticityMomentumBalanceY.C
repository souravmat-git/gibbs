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

#include "PlaneElasticityMomentumBalanceY.h"

registerMooseObject("gibbsApp", PlaneElasticityMomentumBalanceY);

template <>
InputParameters
validParams<PlaneElasticityMomentumBalanceY>()
{
  InputParameters params = validParams<PlaneElasticityBase>();
  params.addRequiredCoupledVar("disp_x","Displacement in the x-direction");
  params.addClassDescription("This kernel implements the momentum balance in y-dir");
  params.addRequiredParam<MaterialPropertyName>("force_y", "Body force in y-dir");

  return params;
}

PlaneElasticityMomentumBalanceY::PlaneElasticityMomentumBalanceY(const InputParameters & parameters)
  :PlaneElasticityBase(parameters),
  _grad_ux(coupledGradient("disp_x")),
  _ux_var(coupled("disp_x")),
  _fy(getMaterialProperty<Real>(getParam<MaterialPropertyName>("force_y")))
{
}

Real
PlaneElasticityMomentumBalanceY::computeQpResidual()
{ 
  //The non-linear variable that this kernel acts on is the displacement in the x-direction;
  //The variable is coupled to the displacement in the y-direction
  //_grad_u[_qp](0) - x-component of the gradient of the displacement in y-direction
  //_grad_u[_qp](1) - y-component of the gradient of the displacement in y-direction
  
  return _grad_test[_i][_qp](0) * (_C66[_qp] * _grad_ux[_qp](1) + _C66[_qp] * _grad_u[_qp](0)) 
      +  _grad_test[_i][_qp](1) * (_C12[_qp] * _grad_ux[_qp](0) + _C22[_qp] * _grad_u[_qp](1))
      -  _test[_i][_qp]*_fy[_qp]  ;
}

Real
PlaneElasticityMomentumBalanceY::computeQpJacobian()
{
  //Derivative with respect to displacement in the y-direction i.e. uy
  return _grad_test[_i][_qp](0) * _C66[_qp] * _grad_phi[_j][_qp](0)
      +  _grad_test[_i][_qp](1) * _C22[_qp] * _grad_phi[_j][_qp](1);
}

Real
PlaneElasticityMomentumBalanceY::computeQpOffDiagJacobian(unsigned int jvar)
{ 
  if (jvar == _ux_var){
    return _grad_test[_i][_qp](0) * _C66[_qp] * _grad_phi[_j][_qp](1)
        +  _grad_test[_i][_qp](1) * _C12[_qp] * _grad_phi[_j][_qp](0);
 }          
 else
   return 0;
}
