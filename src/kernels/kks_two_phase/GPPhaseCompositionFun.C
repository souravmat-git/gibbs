//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This code modifies the KKSPhaseConcentration in MOOSE
//* by including the hand coded interpolation function
//* instead of using the Material property- S.Chatterjee

#include "GPPhaseCompositionFun.h"
#include "Function.h"
registerMooseObject("gibbsApp", GPPhaseCompositionFun);

template <>
InputParameters
validParams<GPPhaseCompositionFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (1-h(eta))*c_alpha + h(eta)*c_beta - c = 0."
                             "non-linear variable of this kernel is c_beta.");
  params.addRequiredParam<FunctionName>("xB", "Mole fraction of comp B");
  params.addRequiredParam<FunctionName>("eta", "To distinguish phases");
  return params;
}

GPPhaseCompositionFun::GPPhaseCompositionFun(const InputParameters & parameters)
  : Kernel(parameters),
   _xB_alpha(getMaterialProperty<Real>("xB_alpha")),
   _xB_beta(getMaterialProperty<Real>("xB_beta")),
   _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
   _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
   _xB(getFunction("xB")),
   _eta(getFunction("eta")),
   //interpolation material
   _h(getMaterialProperty<Real>("h"))
{
}

Real
GPPhaseCompositionFun::computeQpResidual(){
  // R: the non-linear varible that this kernel acts on is muB
  return _test[_i][_qp] * (_h[_qp] *_xB_beta[_qp] + (1.0 - _h[_qp]) * _xB_alpha[_qp]
                           - _xB.value(_t, _q_point[_qp]));
}

Real
GPPhaseCompositionFun::computeQpJacobian(){
  return _test[_i][_qp] * (_h[_qp] * _inv_B_tf_beta[_qp]  + (1.0- _h[_qp])* _inv_B_tf_beta[_qp])
                        * _phi[_j][_qp];
}

Real
GPPhaseCompositionFun::computeQpOffDiagJacobian(unsigned int /*jvar*/){

    return 0.0;
}
