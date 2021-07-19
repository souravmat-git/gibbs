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

#include "QMomentumBalance1DFun.h"
#include "Function.h"
registerMooseObject("gibbsApp", QMomentumBalance1DFun);

template <>
InputParameters
validParams<QMomentumBalance1DFun>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  params.addRequiredParam<FunctionName>("eta", "Phase field variable");
  params.addRequiredParam<MaterialPropertyName>("da_de", "First derivative of jump in strain");
  return params;
}

QMomentumBalance1DFun::QMomentumBalance1DFun(const InputParameters & parameters)
  :Kernel(parameters),
  _eta(getFunction("eta")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _mat_const_beta(getMaterialProperty<Real>("mat_const_beta")),
  _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
  _da_de(getMaterialProperty<Real>(getParam<MaterialPropertyName>("da_de"))),
  _h(getMaterialProperty<Real>("h"))
{
}

Real
QMomentumBalance1DFun::interp_modulus() const{
  return _h[_qp]* _mat_const_beta[_qp] + (1.0-_h[_qp]) * _mat_const_alpha[_qp];
}

Real
QMomentumBalance1DFun::computeQpResidual()
{//The non-linear variable that this kernel acts on is ux
    return -_grad_test[_i][_qp](0) *(_h[_qp] * _sx_beta[_qp] + (1.0-_h[_qp]) * _sx_alpha[_qp]) ;
}

Real
QMomentumBalance1DFun::computeQpJacobian()
{  //The derivative w.r.t displacement field
  return -_grad_test[_i][_qp](0) * (QMomentumBalance1DFun::interp_modulus()
         + _da_de[_qp]*_h[_qp]*(1-_h[_qp])*(_mat_const_beta[_qp] - _mat_const_alpha[_qp]))*_grad_phi[_j][_qp](0);
}

Real
QMomentumBalance1DFun::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
    return 0.0;
}
