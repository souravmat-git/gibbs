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

#include "QMomentumBalance1D_3P.h"
registerMooseObject("gibbsApp", QMomentumBalance1D_3P);

template <>
InputParameters
validParams<QMomentumBalance1D_3P>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  //params.addRequiredCoupledVar("phase_alpha", "Phase field variable");
  //params.addRequiredParam<MaterialPropertyName>("a", "Jump in strain");
  //params.addRequiredParam<MaterialPropertyName>("da_de", "First derivative of jump in strain");
  //params.addRequiredParam<MaterialPropertyName>("da_dh", "First derivative of jump in strain");
  return params;
}

QMomentumBalance1D_3P::QMomentumBalance1D_3P(const InputParameters & parameters)
  :Kernel(parameters),
  _sx_alpha(getMaterialProperty<Real>("alpha_elastic_stress")),
  _mat_const_alpha(getMaterialProperty<Real>("alpha_mat_const"))
{
}

//Real
//QMomentumBalance1D_3P::interp_modulus() const{
//  return (_h[_qp]* _mat_const_beta[_qp] + (1.0-_h[_qp]) * _mat_const_alpha[_qp]);
//}


Real
QMomentumBalance1D_3P::computeQpResidual()
{//The non-linear variable that this kernel acts on is ux
    return -_grad_test[_i][_qp](0) * _sx_alpha[_qp];
}

Real
QMomentumBalance1D_3P::computeQpJacobian()
{  //The derivative w.r.t displacement field
  return - _grad_test[_i][_qp](0) * _mat_const_alpha[_qp] * _grad_phi[_j][_qp](0);
}

Real
QMomentumBalance1D_3P::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
 //f (jvar == _eta_var)
 //{ //Derivative w.r.t phi
//  return -_grad_test[_i][_qp](0) * _dh[_qp] * ( (_sx_beta[_qp]- _sx_alpha[_qp]) -_a[_qp]* QMomentumBalance1D_3P::interp_modulus()
                                // + _da_dh[_qp] * _h[_qp] * (1.0- _h[_qp]) * (_mat_const_beta[_qp] - _mat_const_alpha[_qp]) ) * _phi[_j][_qp];
 //}
 //else
    return 0.0;
}
