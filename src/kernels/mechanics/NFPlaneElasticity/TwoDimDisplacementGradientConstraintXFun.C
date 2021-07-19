//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoDimDisplacementGradientConstraintXFun.h"
#include "Function.h" 
registerMooseObject("gibbsApp", TwoDimDisplacementGradientConstraintXFun);

template <>

InputParameters
validParams<TwoDimDisplacementGradientConstraintXFun>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel acts on ux_beta");
  params.addRequiredCoupledVar("exy_alpha","Strain in the alpha phase");
  params.addRequiredParam<FunctionName>("eta", "Phase-field variable that distinguish phases");
  params.addRequiredCoupledVar("ux", "Displacement in the x-direction");
  params.addRequiredCoupledVar("uy", "Displacement in the y-direction");
  return params;
}

TwoDimDisplacementGradientConstraintXFun::TwoDimDisplacementGradientConstraintXFun(const InputParameters & parameters)
  : Kernel(parameters),
   _eta(getFunction("eta")),
   _grad_ux(coupledGradient("ux")),
   _ux_var(coupled("ux")),
   _grad_uy(coupledGradient("uy")),
   _uy_var(coupled("uy")),
   _exy_alpha(coupledValue("exy_alpha")),
   _exy_alpha_var(coupled("exy_alpha")),
   _exyT_alpha(getMaterialProperty<Real>("exyT_alpha")),
   _exyT_beta(getMaterialProperty<Real>("exyT_beta"))
{
}

Real 
TwoDimDisplacementGradientConstraintXFun::_h() const{

  Real _phi_val = _eta.value(_t, _q_point[_qp]);
  return std::pow(_phi_val,3.0)*(6.0 * std::pow(_phi_val,2.0) - 15.0 * _phi_val + 10.0);
}

Real
TwoDimDisplacementGradientConstraintXFun::_axy() const{

  //compatible strain = elastic_strain + transformation strain
  Real _exyc_beta = _u[_qp] + _exyT_beta[_qp];
  Real _exyc_alpha = _exy_alpha[_qp] + _exyT_alpha[_qp];
  
  Real _exy = 0.5*(_grad_ux[_qp](1) + _grad_uy[_qp](0));

  return (_h() * _exyc_beta
       + (1.0-_h()) * _exyc_alpha - _exy);
}

Real
TwoDimDisplacementGradientConstraintXFun::computeQpResidual(){  
  return _test[_i][_qp]* _axy();
}

Real
TwoDimDisplacementGradientConstraintXFun::computeQpJacobian(){
 //Derivative wrt compataibl strain of the beta phase
  return _test[_i][_qp] * _h() * _phi[_j][_qp];
}

Real
TwoDimDisplacementGradientConstraintXFun::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _exy_alpha_var){
    return _test[_i][_qp] * (1.0-_h()) * _phi[_j][_qp];
  }
  else if (jvar == _ux_var){
    return -_test[_i][_qp] * _grad_phi[_j][_qp](1);
  }
  else if (jvar == _uy_var){
    return - _test[_i][_qp] * _grad_phi[_j][_qp](0);
  }
  else
    return 0.0;
}
