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

#include "ContinuityEquationWBM.h"
registerMooseObject("gibbsApp", ContinuityEquationWBM);

template <>
InputParameters
validParams<ContinuityEquationWBM>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("eta", "Phase field variable"),
  params.addRequiredParam<MaterialPropertyName>("diff_mobility", "Diffusion mobility");
  return params;
}

ContinuityEquationWBM::ContinuityEquationWBM(const InputParameters & parameters)
  : Kernel(parameters),
  //_eta(coupledValue("eta")),
  _grad_eta(coupledGradient("eta")),
  _eta_var(coupled("eta")),
  //_diff_pot(coupledValue("B_diff_pot")),
  //_diff_pot_var(coupled("B_diff_pot")),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  //Material properties
  _dfalpha_dc(getMaterialProperty<Real>("B_diff_pot_alpha")),
  _d2falpha_dc2(getMaterialProperty<Real>("B_therm_factor_alpha")),
  _dfbeta_dc(getMaterialProperty<Real>("B_diff_pot_beta")),
  _d2fbeta_dc2(getMaterialProperty<Real>("B_therm_factor_beta")),
  _M(getMaterialProperty<Real>("diff_mobility"))
{
}

Real
ContinuityEquationWBM::interpolated_tf() const {
  return (_d2fbeta_dc2[_qp]*_h[_qp] + _d2falpha_dc2[_qp]* (1.0- _h[_qp]));
}

RealGradient
ContinuityEquationWBM::grad_mu() const{
   //Here, u is the mole fraction variable,
   //Note that mu is a function of both c (or u) and phi
   //The first term is the derivative of mu w.r.t c
   //and, the second term is due to derivative of phi
  return  ContinuityEquationWBM::interpolated_tf() * _grad_u[_qp]
        + (_dfbeta_dc[_qp] - _dfalpha_dc[_qp]) * _dh[_qp] * _grad_eta[_qp];
}

Real
ContinuityEquationWBM::computeQpResidual(){
  return (_grad_test[_i][_qp] *  _M[_qp] * ContinuityEquationWBM::grad_mu());
}

Real
ContinuityEquationWBM::computeQpJacobian(){
  //Here, the derivative is with respect to the mole fraction variable
  Real jac1_xB = _grad_test[_i][_qp] * _M[_qp] * ContinuityEquationWBM::interpolated_tf() * _grad_phi[_j][_qp];
  Real jac2_xB = _grad_test[_i][_qp] * _M[_qp] * (_d2fbeta_dc2[_qp] - _d2falpha_dc2[_qp]) * _dh[_qp] * _grad_eta[_qp] * _phi[_j][_qp];
  
  return (jac1_xB  + jac2_xB);
}

Real
ContinuityEquationWBM::computeQpOffDiagJacobian(unsigned int jvar){
 
  if (jvar == _eta_var){
  //Derivative of the residual wrt phi (assuming constant mobility) will have three terms 
  //The first term is due to ContinuityEquationWBM::interpolated_tf(), the second due to dh in second residual
  Real jac1_phi =  _grad_test[_i][_qp] * _M[_qp] * (_d2fbeta_dc2[_qp] - _d2falpha_dc2[_qp])* _dh[_qp]  * _grad_u[_qp]   * _phi[_j][_qp];
  Real jac2_phi =  _grad_test[_i][_qp] * _M[_qp] * (_dfbeta_dc[_qp] - _dfalpha_dc[_qp])    * _d2h[_qp] * _grad_eta[_qp] * _phi[_j][_qp];
  Real jac3_phi =  _grad_test[_i][_qp] * _M[_qp] * (_dfbeta_dc[_qp] - _dfalpha_dc[_qp])    * _dh[_qp]  * _grad_phi[_j][_qp];
    
    return (jac1_phi + jac2_phi + jac3_phi);
  }
  else
    return 0;  
}
