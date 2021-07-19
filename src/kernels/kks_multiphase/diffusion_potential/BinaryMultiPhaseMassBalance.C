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

#include "BinaryMultiPhaseMassBalance.h"

registerMooseObject("gibbsApp", BinaryMultiPhaseMassBalance);

template <>
InputParameters
validParams<BinaryMultiPhaseMassBalance>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("phase_alpha", "phase 1");
  params.addRequiredCoupledVar("phase_beta", "Phase 2");
  params.addCoupledVar("phase_gamma", 0.0, "Phase 3");
  params.addCoupledVar("phase_delta", 0.0, "Phase 4");
  params.addCoupledVar("phase_epsilon", 0.0, "Phase 5");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");    
  //params.addParam<MaterialPropertyName>("inv_B_td_gamma", 0.0, "Third derivative in gamma phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_delta", 0.0, "Thermodynamic factor in delta phase");
  //params.addParam<MaterialPropertyName>("inv_B_td_delta", 0.0, "Third derivative in delta phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_epsilon", 0.0, "Thermodynamic factor in epsilon phase");
  //params.addParam<MaterialPropertyName>("inv_B_td_epsilon", 0.0, "Third derivative in epsilon phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_gamma", 0.0, "Thermodynamic factor in gamma phase");
  //params.addParam<MaterialPropertyName>("inv_B_td_gamma", 0.0, "Third derivative in gamma phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_delta", 0.0, "Thermodynamic factor in delta phase");
  //params.addParam<MaterialPropertyName>("inv_B_td_delta", 0.0, "Third derivative in delta phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_epsilon", 0.0, "Thermodynamic factor in epsilon phase");
  //params.addParam<MaterialPropertyName>("inv_B_td_epsilon", 0.0, "Third derivative in epsilon phase");
  params.addRequiredParam<MaterialPropertyName>("h_alpha", "interpolation");
  params.addRequiredParam<MaterialPropertyName>("h_beta", "interpolation");
  params.addParam<MaterialPropertyName>("h_gamma",0.0, "interpolation");
  params.addParam<MaterialPropertyName>("h_delta",0.0, "interpolation");
  params.addParam<MaterialPropertyName>("h_epsilon",0.0, "interpolation");
  return params; 
}

BinaryMultiPhaseMassBalance::BinaryMultiPhaseMassBalance(const InputParameters & parameters)
  : Kernel(parameters),
  _phase_alpha(coupledValue("phase_alpha")),
  _phase_alpha_var(coupled("phase_alpha")),
  _phase_beta(coupledValue("phase_beta")),
  _phase_beta_var(coupled("phase_beta")),
  _phase_gamma(coupledValue("phase_gamma")),
  _phase_gamma_var(coupled("phase_gamma")),
  _phase_delta(coupledValue("phase_delta")),
  _phase_delta_var(coupled("phase_delta")),
  _phase_epsilon(coupledValue("phase_epsilon")),
  _phase_epsilon_var(coupled("phase_epsilon")),
  //component B
  _grad_B_diff_pot(coupledGradient("B_diff_pot")),
  _B_diff_pot_var(coupled("B_diff_pot")),  
  //inverse of thermodynamic factor
  _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
  _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
  _inv_B_tf_gamma(getMaterialProperty<Real>("inv_B_tf_gamma")),
  _inv_B_tf_delta(getMaterialProperty<Real>("inv_B_tf_delta")),
  _inv_B_tf_epsilon(getMaterialProperty<Real>("inv_B_tf_epsilon")),
  //interpolation material
  _h_alpha_name(getParam<MaterialPropertyName>("h_alpha")),
  _h_beta_name(getParam<MaterialPropertyName>("h_beta")),
  _h_gamma_name(getParam<MaterialPropertyName>("h_gamma")),
  _h_delta_name(getParam<MaterialPropertyName>("h_delta")),
  _h_epsilon_name(getParam<MaterialPropertyName>("h_epsilon")),
  _h_alpha(getMaterialProperty<Real>(_h_alpha_name)),
  _h_beta(getMaterialProperty<Real>(_h_beta_name)),
  _h_gamma(getMaterialProperty<Real>(_h_gamma_name)),
  _h_delta(getMaterialProperty<Real>(_h_delta_name)),
  _h_epsilon(getMaterialProperty<Real>(_h_epsilon_name)),
  // for coupled variable phase_alpha
  _dhbeta_dphialpha(getMaterialProperty<Real>("dhbeta_dphialpha")),
  _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),
  _dhdelta_dphialpha(getMaterialProperty<Real>("dhdelta_dphialpha")),
  _dhepsilon_dphialpha(getMaterialProperty<Real>("dhepsilon_dphialpha")),
      
  // for coupled variable phase_beta
  _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
  _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")),
  _dhdelta_dphibeta(getMaterialProperty<Real>("dhdelta_dphibeta")),
  _dhepsilon_dphibeta(getMaterialProperty<Real>("dhepsilon_dphibeta")),
  
  // for coupled variable phase_gamma
  _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
  _dhbeta_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
  _dhdelta_dphigamma(getMaterialProperty<Real>("dhdelta_dphigamma")),
  _dhepsilon_dphigamma(getMaterialProperty<Real>("dhepsilon_dphigamma")),
  
  // for coupled variable phase_delta
  _dhalpha_dphidelta(getMaterialProperty<Real>("dhalpha_dphidelta")),
  _dhbeta_dphidelta(getMaterialProperty<Real>("dhbeta_dphidelta")),
  _dhgamma_dphidelta(getMaterialProperty<Real>("dhgamma_dphidelta")),
  _dhepsilon_dphidelta(getMaterialProperty<Real>("dhepsilon_dphidelta")),
  
   // for coupled variable phase_epsilon
  _dhalpha_dphiepsilon(getMaterialProperty<Real>("dhalpha_dphiepsilon")),
  _dhbeta_dphiepsilon(getMaterialProperty<Real>("dhbeta_dphiepsilon")),
  _dhgamma_dphiepsilon(getMaterialProperty<Real>("dhgamma_dphiepsilon")),
  _dhdelta_dphiepsilon(getMaterialProperty<Real>("dhdelta_dphiepsilon")),
  
  //Mobility within each phase
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  _L_BB_gamma(getMaterialProperty<Real>("L_BB_gamma")),
  _L_BB_delta(getMaterialProperty<Real>("L_BB_delta")),
  _L_BB_epsilon(getMaterialProperty<Real>("L_BB_epsilon")),
  _dL_BB_muB_alpha(getMaterialProperty<Real>("dL_BB_muB_alpha")),
  _dL_BB_muB_beta(getMaterialProperty<Real>("dL_BB_muB_beta")),
  _dL_BB_muB_gamma(getMaterialProperty<Real>("dL_BB_muB_gamma")),
  _dL_BB_muB_delta(getMaterialProperty<Real>("dL_BB_muB_delta")),
  _dL_BB_muB_epsilon(getMaterialProperty<Real>("dL_BB_muB_epsilon"))  
  //Only required if thermodynamic factor varies with diffusion potential
  //_inv_B_td_alpha(getMaterialProperty<Real>("inv_B_td_alpha")),
  //_inv_B_td_beta(getMaterialProperty<Real>("inv_B_td_beta")),
  //_inv_B_td_gamma(getMaterialProperty<Real>("inv_B_td_gamma")),
  //_inv_B_td_delta(getMaterialProperty<Real>("inv_B_td_delta")),
  //_inv_B_td_epsilon(getMaterialProperty<Real>("inv_B_td_epsilon"))
{
}

Real
BinaryMultiPhaseMassBalance::thermodynamic_factor() const
{
 //the thermodynamic factor is same as: BinaryPhaseConstraint::computeQpDiagJac
  return (std::pow((  _h_alpha[_qp] * _inv_B_tf_alpha[_qp] + _h_beta[_qp] * _inv_B_tf_beta[_qp]
                    + _h_gamma[_qp] * _inv_B_tf_gamma[_qp]
                    + _h_delta[_qp] * _inv_B_tf_delta[_qp]
                    + _h_epsilon[_qp] * _inv_B_tf_epsilon[_qp]), -1.0));
}

//Real
//BinaryMultiPhaseMassBalance::third_deriv() const
//{
 //this is rate of change of the curvature of the free energy curve with composition
  //return (std::pow((  _h_alpha[_qp] * _inv_B_td_alpha[_qp] + _h_beta[_qp] * _inv_B_td_beta[_qp]
  //                  + _h_gamma[_qp] * _inv_B_td_gamma[_qp]
  //                  + _h_delta[_qp] * _inv_B_td_delta[_qp]
  //                  + _h_epsilon[_qp] * _inv_B_td_epsilon[_qp]), -1.0 ));
//}

Real
BinaryMultiPhaseMassBalance::L_BB_interp() const
{
  return (_L_BB_alpha[_qp]*_h_alpha[_qp] + _L_BB_beta[_qp]* _h_beta[_qp] 
        + _L_BB_gamma[_qp]* _h_gamma[_qp] +  _L_BB_delta[_qp]* _h_delta[_qp]
        + _L_BB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
BinaryMultiPhaseMassBalance::dL_BB_muB_interp() const
{
  return ( _dL_BB_muB_alpha[_qp]*_h_alpha[_qp] + _dL_BB_muB_beta[_qp]* _h_beta[_qp] 
          + _dL_BB_muB_gamma[_qp]* _h_gamma[_qp] + _dL_BB_muB_delta[_qp] * _h_delta[_qp]
          + _dL_BB_muB_epsilon[_qp]* _h_epsilon[_qp]);
}

Real
BinaryMultiPhaseMassBalance::computeQpResidual()
{   
  //Variable on which the kernel operates: X_B
  return (_grad_test[_i][_qp] * BinaryMultiPhaseMassBalance::L_BB_interp() * _grad_B_diff_pot[_qp]);
}

Real
BinaryMultiPhaseMassBalance::computeQpJacobian()
{
  return _grad_test[_i][_qp] * (BinaryMultiPhaseMassBalance::L_BB_interp() * 
                                 BinaryMultiPhaseMassBalance::thermodynamic_factor() * _grad_phi[_j][_qp]);
  
}

Real
BinaryMultiPhaseMassBalance::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _B_diff_pot_var)
  { 
    return (_grad_test[_i][_qp] * ( (BinaryMultiPhaseMassBalance::L_BB_interp() * _grad_phi[_j][_qp])
                                +  (BinaryMultiPhaseMassBalance::dL_BB_muB_interp()* _grad_B_diff_pot[_qp]) * _phi[_j][_qp]) );
  }
 else if (jvar == _phase_alpha_var)
 {
    return (_grad_test[_i][_qp] *(_dhbeta_dphialpha[_qp] * (_L_BB_beta[_qp] - _L_BB_alpha[_qp])
                                + _dhgamma_dphialpha[_qp]* (_L_BB_gamma[_qp]- _L_BB_alpha[_qp])
                                + _dhdelta_dphialpha[_qp]* (_L_BB_delta[_qp]- _L_BB_alpha[_qp])
                              + _dhepsilon_dphialpha[_qp]*(_L_BB_epsilon[_qp]- _L_BB_alpha[_qp]))* _grad_B_diff_pot[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _phase_beta_var)
 {
   return (_grad_test[_i][_qp] * (_dhalpha_dphibeta[_qp] * (_L_BB_alpha[_qp] - _L_BB_beta[_qp]) 
                                 +_dhgamma_dphibeta[_qp] * (_L_BB_gamma[_qp] - _L_BB_beta[_qp])
                                 +_dhdelta_dphibeta[_qp] * (_L_BB_delta[_qp] - _L_BB_beta[_qp])
                              + _dhepsilon_dphibeta[_qp]* (_L_BB_epsilon[_qp] - _L_BB_beta[_qp]))* _grad_B_diff_pot[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _phase_gamma_var)
 {
   return (_grad_test[_i][_qp] *(_dhalpha_dphigamma[_qp] * (_L_BB_alpha[_qp] - _L_BB_gamma[_qp]) 
                                +_dhbeta_dphigamma[_qp] *  (_L_BB_beta[_qp]  - _L_BB_gamma[_qp])
                                +_dhdelta_dphigamma[_qp] * (_L_BB_delta[_qp] - _L_BB_gamma[_qp])
                              + _dhepsilon_dphigamma[_qp]* (_L_BB_epsilon[_qp] - _L_BB_gamma[_qp]))* _grad_B_diff_pot[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _phase_delta_var)
 {
   return (_grad_test[_i][_qp] *(_dhalpha_dphidelta[_qp] * (_L_BB_alpha[_qp] - _L_BB_delta[_qp]) 
                                +_dhbeta_dphidelta[_qp] *  (_L_BB_beta[_qp]  - _L_BB_delta[_qp])
                                +_dhgamma_dphidelta[_qp] * (_L_BB_gamma[_qp] - _L_BB_delta[_qp])
                             + _dhepsilon_dphidelta[_qp]* (_L_BB_epsilon[_qp]- _L_BB_delta[_qp]))* _grad_B_diff_pot[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _phase_epsilon_var)
 {
   return (_grad_test[_i][_qp] *(_dhalpha_dphiepsilon[_qp] * (_L_BB_alpha[_qp] - _L_BB_epsilon[_qp]) 
                                +_dhbeta_dphiepsilon[_qp] *  (_L_BB_beta[_qp]  - _L_BB_epsilon[_qp])
                                +_dhgamma_dphiepsilon[_qp] * (_L_BB_gamma[_qp] - _L_BB_epsilon[_qp])
                             +_dhdelta_dphiepsilon[_qp]*  (_L_BB_delta[_qp]-_L_BB_epsilon[_qp]))* _grad_B_diff_pot[_qp] * _phi[_j][_qp]);
 }
 else
    return 0;
}
