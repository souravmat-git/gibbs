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

#include "LCMomentumBalanceX.h"
registerMooseObject("gibbsApp", LCMomentumBalanceX);

template <>
InputParameters
validParams<LCMomentumBalanceX>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the momentum balance in x-dir");
  params.addRequiredCoupledVar("eta", "Phase field variable");
  params.addRequiredCoupledVar("stress_x", "Stress in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("ex_beta", "Strain in the beta phase" );
  params.addRequiredParam<MaterialPropertyName>("ex_alpha", "Strain in the alpha phase");
  params.addRequiredParam<MaterialPropertyName>("compliance_beta", "Youngs Modulus of beta phase");
  params.addRequiredParam<MaterialPropertyName>("compliance_alpha", "Youngs Modulus in alpha phase");
  return params;
}

LCMomentumBalanceX::LCMomentumBalanceX(const InputParameters & parameters)
  :Kernel(parameters),
  _grad_eta(coupledGradient("eta")),
  _eta_var(coupled("eta")),
  _stress_x(coupledValue("stress_x")),
  _stress_x_var(coupled("stress_x")),
  //Material properties
  _ex_beta(getMaterialProperty<Real>(getParam<MaterialPropertyName>("ex_beta"))),
  _ex_alpha(getMaterialProperty<Real>(getParam<MaterialPropertyName>("ex_beta"))),
  _comp_beta(getMaterialProperty<Real>(getParam<MaterialPropertyName>("compliance_beta"))),
  _comp_alpha(getMaterialProperty<Real>(getParam<MaterialPropertyName>("compliance_alpha"))),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h"))
{
}

Real
LCMomentumBalanceX::compliance_inv() const{
  return std::pow((_h[_qp]*_comp_beta[_qp] + (1.0-_h[_qp])*_comp_alpha[_qp]),-1);
}

Real
LCMomentumBalanceX::dcomplinv_dh() const{
  return -(_comp_beta[_qp] - _comp_alpha[_qp]) * std::pow((_h[_qp]* _comp_beta[_qp] + (1.0-_h[_qp])*_comp_alpha[_qp]) ,-2);
}

Real
LCMomentumBalanceX::d2complinv_dh2() const{
  return  2.0*std::pow((_comp_beta[_qp] - _comp_alpha[_qp]),2)* 
            std::pow((_h[_qp]* _comp_beta[_qp] + (1.0-_h[_qp])*_comp_alpha[_qp]) ,-3);
}


Real
LCMomentumBalanceX::computeQpResidual(){   
   //This kernel acts on the displacement variable    
   //res1  for homogeneous elasticty
   //res2 and res3 for inhomohegenous case
   //Calculate the residuals
    
   Real res1 = -_grad_test[_i][_qp](0) * _grad_u[_qp](0);
   Real res2 = -_test[_i][_qp]*(_ex_beta[_qp]- _ex_alpha[_qp]) * _dh[_qp] * _grad_eta[_qp](0); 
   
    return (res1 + res2);         
}

Real
LCMomentumBalanceX::computeQpJacobian(){
  //Derivative w.r.t displacement variable
  
    Real jac1_ux = -_grad_test[_i][_qp](0) * _grad_phi[_j][_qp](0);
  
  return -(jac1_ux);
}

Real
LCMomentumBalanceX::computeQpOffDiagJacobian(unsigned int jvar){
 if (jvar == _eta_var){   
  Real jac_phi = -_test[_i][_qp] * (_ex_beta[_qp] - _ex_alpha[_qp])* 
                    (_grad_eta[_qp](0)* _d2h[_qp] * _phi[_j][_qp] + _dh[_qp] * _grad_phi[_j][_qp](0));
 
 return (jac_phi);
        
  }             
 else if (jvar == _stress_x_var){
   return -_grad_test[_i][_qp](0) * _phi[_j][_qp];
 }
 else
  return 0;         
}
