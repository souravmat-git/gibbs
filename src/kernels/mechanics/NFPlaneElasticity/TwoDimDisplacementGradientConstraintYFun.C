//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoDimDisplacementGradientConstraintYFun.h"
#include "Function.h" 
registerMooseObject("gibbsApp", TwoDimDisplacementGradientConstraintYFun);

template <>

InputParameters
validParams<TwoDimDisplacementGradientConstraintYFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel acts on ux_beta");
  params.addRequiredCoupledVar("ex_alpha","Strain in the alpha phase");
  params.addRequiredParam<FunctionName>("eta", "Phase-field variable that distinguish phases");
  params.addRequiredCoupledVar("ux", "Displacement in the x-direction");
  return params;
}

TwoDimDisplacementGradientConstraintYFun::TwoDimDisplacementGradientConstraintYFun(const InputParameters & parameters)
  : Kernel(parameters),
   _eta(getFunction("eta")),
   _grad_uy(coupledGradient("uy")),
   _uy_var(coupled("uy")),
   _ey_alpha(coupledValue("ey_alpha")),
   _ey_alpha_var(coupled("ey_alpha")),
   _eyT_alpha(getMaterialProperty<Real>("eyT_alpha")),
   _eyT_beta(getMaterialProperty<Real>("eyT_beta"))
{
}

Real 
TwoDimDisplacementGradientConstraintYFun::_h() const{

  Real _phi_val = _eta.value(_t, _q_point[_qp]);
  return std::pow(_phi_val,3.0)*(6.0 * std::pow(_phi_val,2.0) - 15.0 * _phi_val + 10.0);
}

Real
TwoDimDisplacementGradientConstraintYFun::_ay() const{

  //compatible strain = elastic_strain + transformation strain
  Real _eyc_beta = _u[_qp] + _eyT_beta[_qp];
  Real _eyc_alpha = _ey_alpha[_qp] + _eyT_alpha[_qp];

  return     (_h() * _eyc_beta
      + (1.0-_h()) * _eyc_alpha - _grad_uy[_qp](1));
}

Real
TwoDimDisplacementGradientConstraintYFun::computeQpResidual(){  
  return _test[_i][_qp]* _ay();
}

Real
TwoDimDisplacementGradientConstraintYFun::computeQpJacobian(){
 //Derivative wrt compataibl strain of the beta phase
  return _test[_i][_qp] * _h() * _phi[_j][_qp];
}

Real
TwoDimDisplacementGradientConstraintYFun::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _ey_alpha_var){
    return _test[_i][_qp] * (1.0-_h()) * _phi[_j][_qp];
  }
  else if (jvar == _uy_var){
    return -_test[_i][_qp] * _grad_phi[_j][_qp](1);
  }
  else
    return 0.0;
}
