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

#include "KKSPhaseCompositionFun.h"
#include "Function.h"
registerMooseObject("gibbsApp", KKSPhaseCompositionFun);

template <>
InputParameters
validParams<KKSPhaseCompositionFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (1-h(eta))*c_alpha + h(eta)*c_beta - c = 0."
                             "non-linear variable of this kernel is c_beta.");
  params.addRequiredCoupledVar("xB_alpha", "Phase concentration in alpha phase");
  params.addRequiredParam<FunctionName>("mole_fraction", "Concentration function");
  params.addRequiredParam<FunctionName>("eta", "To distinguish phases");
  return params;
}

KKSPhaseCompositionFun::KKSPhaseCompositionFun(const InputParameters & parameters)
  :  Kernel(parameters),
    _xB_alpha(coupledValue("xB_alpha")),
    _xB_alpha_var(coupled("xB_alpha")),
    _mole_fraction(getFunction("mole_fraction")),
    //_mole_fraction_var(coupled("mole_fraction")),
    _eta(getFunction("eta"))
    //_eta_var(coupled("eta")),
    //_h(getMaterialProperty<Real>("h")),
    //_dh(getMaterialProperty<Real>("dh"))
{
}

Real 
KKSPhaseCompositionFun::_h() const{

  Real _phi_val = _eta.value(_t, _q_point[_qp]);
  return std::pow(_phi_val,3.0)*(6.0 * std::pow(_phi_val,2.0) - 15.0 * _phi_val + 10);
}

Real
KKSPhaseCompositionFun::computeQpResidual(){
  // R: w*((1-h(eta))*ca + h(eta)*cb - c)
  return _test[_i][_qp] * (_h() *_u[_qp] + (1.0 - _h()) * _xB_alpha[_qp] 
         - _mole_fraction.value(_t, _q_point[_qp]));
}

Real
KKSPhaseCompositionFun::computeQpJacobian(){
  return _test[_i][_qp] * _h() * _phi[_j][_qp];
}

Real
KKSPhaseCompositionFun::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xB_alpha_var){
    return _test[_i][_qp] * (1.0 - _h()) * _phi[_j][_qp];
  }
  //else if (jvar == _mole_fraction_var)
  //{
  //  return -_test[_i][_qp] * _phi[_j][_qp];
  //}
  //else if (jvar == _eta_var)
  //{
  //  return _test[_i][_qp] * (_u[_qp] - _xB_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
  //}
  else
    return 0.0;
}
