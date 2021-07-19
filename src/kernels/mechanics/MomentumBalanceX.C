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

#include "MomentumBalanceX.h"

registerMooseObject("gibbsApp", MomentumBalanceX);

template <>
InputParameters
validParams<MomentumBalanceX>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  params.addRequiredParam<MaterialPropertyName>("youngs_modulus", "Youngs Modulus");
  params.addRequiredParam<MaterialPropertyName>("body_force_x", "Body force in x-dir");
  params.addRequiredParam<MaterialPropertyName>("lc", "Characteristic length scale of the problem");
  return params;
}

MomentumBalanceX::MomentumBalanceX(const InputParameters & parameters)
  :Kernel(parameters),
  _E(getMaterialProperty<Real>(getParam<MaterialPropertyName>("youngs_modulus"))),
  _bx(getMaterialProperty<Real>(getParam<MaterialPropertyName>("body_force_x"))),
  _lc(getMaterialProperty<Real>(getParam<MaterialPropertyName>("lc")))
{
}

Real
MomentumBalanceX::computeQpResidual()
{
  //The non-linear variable that this kernel acts on is u_x;
  
  return (-_E[_qp]*_lc[_qp]*(_grad_test[_i][_qp] * _grad_u[_qp]) 
                       +(_test[_i][_qp] * _lc[_qp]* _lc[_qp]* _bx[_qp])) ;
}

Real
MomentumBalanceX::computeQpJacobian()
{
  //Derivative with respect to the non-linear variable
  return (-_E[_qp]*_lc[_qp]* _grad_test[_i][_qp] * _grad_phi[_j][_qp]);
}

Real
MomentumBalanceX::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{ 
    return 0;
}
