//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacementGradientConstraintXFun.h"
#include "Function.h" 
registerMooseObject("gibbsApp", DisplacementGradientConstraintXFun);

template <>

InputParameters
validParams<DisplacementGradientConstraintXFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel acts on ux_beta");
  params.addRequiredCoupledVar("ex_alpha","Strain in the alpha phase");
  params.addRequiredParam<FunctionName>("eta", "Phase-field variable that distinguish phases");
  params.addRequiredCoupledVar("ux", "Displacement in the x-direction");
  return params;
}

DisplacementGradientConstraintXFun::DisplacementGradientConstraintXFun(const InputParameters & parameters)
  : Kernel(parameters),
   _eta(getFunction("eta")),
   _grad_ux(coupledGradient("ux")),
   _ux_var(coupled("ux")),
   _ex_alpha(coupledValue("ex_alpha")),
   _ex_alpha_var(coupled("ex_alpha")),
   _eT_alpha(getMaterialProperty<Real>("eT_alpha")),
   _eT_beta(getMaterialProperty<Real>("eT_beta"))
{
}

Real 
DisplacementGradientConstraintXFun::_h() const{

  Real _phi_val = _eta.value(_t, _q_point[_qp]);
  return std::pow(_phi_val,3.0)*(6.0 * std::pow(_phi_val,2.0) - 15.0 * _phi_val + 10.0);
}

Real
DisplacementGradientConstraintXFun::_ax() const{

  //compatible strain = elastic_strain + transformation strain
  Real _ec_beta = _u[_qp] + _eT_beta[_qp];
  Real _ec_alpha = _ex_alpha[_qp] + _eT_alpha[_qp];

  return     (_h() * _ec_beta
      + (1.0-_h()) * _ec_alpha - _grad_ux[_qp](0));
}

Real
DisplacementGradientConstraintXFun::computeQpResidual(){  
  return _test[_i][_qp]* _ax();
}

Real
DisplacementGradientConstraintXFun::computeQpJacobian(){
 //Derivative wrt compataibl strain of the beta phase
  return _test[_i][_qp] * _h() * _phi[_j][_qp];
}

Real
DisplacementGradientConstraintXFun::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _ex_alpha_var){
    return _test[_i][_qp] * (1.0-_h()) * _phi[_j][_qp];
  }
  else if (jvar == _ux_var){
    return -_test[_i][_qp] * _grad_phi[_j][_qp](0);
  }
  else
    return 0.0;
}
