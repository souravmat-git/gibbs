//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacementGradientConstraintX.h"
registerMooseObject("gibbsApp", DisplacementGradientConstraintX);

template <>

InputParameters
validParams<DisplacementGradientConstraintX>()
{
  InputParameters params = validParams<ObtainUnitNormalBase>();
  params.addClassDescription("This kernel acts on ux_beta");
  params.addRequiredCoupledVar("ex_alpha","Strain in the alpha phase");
  params.addRequiredCoupledVar("ux", "Overall displacement");
  params.addRequiredCoupledVar("eta", "Phase-field variable that distinguish phases");
  return params;
}

DisplacementGradientConstraintX::DisplacementGradientConstraintX(const InputParameters & parameters)
  : ObtainUnitNormalBase(parameters),
   _eta(coupledValue("eta")),
   _eta_var(coupled("eta")),
   _ex_alpha(coupledValue("ex_alpha")),
   _ex_alpha_var(coupled("ex_alpha")),
   _grad_ux(coupledGradient("ux")),
   _ux_var(coupled("ux")),
   _eT_alpha(getMaterialProperty<Real>("eT_alpha")),
   _eT_beta(getMaterialProperty<Real>("eT_beta")),
   _h(getMaterialProperty<Real>("h")),
   _dh(getMaterialProperty<Real>("dh"))
{
}

Real
DisplacementGradientConstraintX::axx() const{
  
  Real _ec_beta = _u[_qp] + _eT_beta[_qp];  
  Real _ec_alpha = _ex_alpha[_qp] + _eT_alpha[_qp];
  
  return     (_h[_qp]  * _ec_beta
      + (1.0-_h[_qp])  * _ec_alpha - _grad_ux[_qp](0));
}

Real
DisplacementGradientConstraintX::computeQpResidual(){  
  return _test[_i][_qp]*axx()*ObtainUnitNormalBase::nx();
}

Real
DisplacementGradientConstraintX::computeQpJacobian(){
  return _test[_i][_qp] * _h[_qp] * _phi[_j][_qp]* ObtainUnitNormalBase::nx();
}

Real
DisplacementGradientConstraintX::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _ex_alpha_var){
    return _test[_i][_qp] * (1.0-_h[_qp]) * _phi[_j][_qp]* ObtainUnitNormalBase::nx();
  }
  else if (jvar == _ux_var){
    return - _test[_i][_qp] * _grad_phi[_j][_qp](0)* ObtainUnitNormalBase::nx();
  }
  else if (jvar == _eta_var){

    Real _ec_beta  = _u[_qp] +  _eT_beta[_qp];
    Real _ec_alpha = _ex_alpha[_qp] + _eT_alpha[_qp];
      
    return _test[_i][_qp] * (_ec_beta - _ec_alpha) * _dh[_qp] * _phi[_j][_qp]* ObtainUnitNormalBase::nx();
  }
  else
    return 0.0;
}
