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

#include "ContinuityEquationKKS.h"

registerMooseObject("gibbsApp", ContinuityEquationKKS);

template <>
InputParameters
validParams<ContinuityEquationKKS>(){
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  //params.addRequiredCoupledVar("xB",  "Overall mole fraction");
  params.addRequiredCoupledVar("xB_beta",  " Phase Mole fraction of component B");
  return params;
}

ContinuityEquationKKS::ContinuityEquationKKS(const InputParameters & parameters)
 : Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  //_grad_xB(coupledGradient("xB")),
  //_xB_var(coupled("xB")),
  _grad_xB_beta(coupledGradient("xB_beta")),
  _xB_beta_var(coupled("xB_beta")),
  //thermodynamic_factor of both phases
  //_B_tf_alpha(getMaterialProperty<Real>("B_therm_factor_alpha")),
  _B_tf_beta(getMaterialProperty<Real>("B_therm_factor_beta")),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  //Kinetic materials
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  //The mobilities are functions of phase compositions
  _dL_BB_xB_alpha(getMaterialProperty<Real>("dL_BB_xB_alpha")),
  _dL_BB_xB_beta(getMaterialProperty<Real>("dL_BB_xB_beta"))
{
}

//To check the implementation of f_BB see Eqn(29) in Kim's paper
//This term can be interpreted as the overall thermodyanmic factor
//Real
//ContinuityEquationKKS::f_BB() const {
//  return _B_tf_alpha[_qp]*_B_tf_beta[_qp]
//         *(std::pow(_B_tf_alpha[_qp] * _h[_qp] + _B_tf_beta[_qp] * (1-_h[_qp]),-1));
//}

Real
ContinuityEquationKKS::L_BB_interp() const {
  return (_L_BB_beta[_qp] * _h[_qp] + _L_BB_alpha[_qp] * (1.0 - _h[_qp]));
}

Real
ContinuityEquationKKS::dL_BB_xB_beta_interp() const{
  return (_dL_BB_xB_beta[_qp] * _h[_qp] + _dL_BB_xB_alpha[_qp] * (1.0 - _h[_qp])); 
}

Real
ContinuityEquationKKS::computeQpResidual(){
   //The nonlinear variable that this kerne acts on the is the overall 
   //diffusion potential which is equal to phase diffusion potential
   // _B_diff_pot = _B_diff_pot_beta = _B_diff_pot_alpha
   return (_grad_test[_i][_qp] * ContinuityEquationKKS::L_BB_interp() * _grad_u[_qp]);
}

Real
ContinuityEquationKKS::computeQpJacobian(){
  //Take the derivative wrt diffusion potential of comp B
  return (_grad_test[_i][_qp] * ContinuityEquationKKS::L_BB_interp() * _grad_phi[_j][_qp]);
}

Real
ContinuityEquationKKS::computeQpOffDiagJacobian(unsigned int jvar){ 

  if (jvar == _eta_var)
  {
    return (_grad_test[_i][_qp]* _dh[_qp]*(_L_BB_beta[_qp] -_L_BB_alpha[_qp])*_grad_u[_qp]*_phi[_j][_qp]);
  }
  else if (jvar == _xB_beta_var)
  {
    //Note that the derivative is wrt phase compositions and not wrt the overall composition
    return (_grad_test[_i][_qp]* ContinuityEquationKKS::L_BB_interp()*_B_tf_beta[_qp]*_grad_xB_beta[_qp]
          + _grad_test[_i][_qp]* ContinuityEquationKKS::dL_BB_xB_beta_interp()*_grad_u[_qp]*_phi[_j][_qp]);
  }
  //else if (jvar == _xB_var)
  //{
  //  return _grad_test[_i][_qp] * ContinuityEquationKKS::L_BB_interp()* ContinuityEquationKKS::f_BB()*_grad_phi[_j][_qp];
  //}
  else
  {
    return 0;
  }    
}
