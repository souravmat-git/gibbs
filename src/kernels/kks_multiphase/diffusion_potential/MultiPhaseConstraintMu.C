//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This code was written by S.Chatterjee

#include "MultiPhaseConstraintMu.h"
registerMooseObject("gibbsApp", MultiPhaseConstraintMu);

//*class template specialization
template <>
InputParameters
validParams<MultiPhaseConstraintMu>()
{    
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: c = c_{alpha}h_{alpha} + ...\
                              +c_{beta}h_{beta} + c_{gamma}h_{gamma}");
  //The non-variable that this kernel operates on is c_{gamma}                       
  params.addRequiredCoupledVar("phase_alpha", "phase field for alpha phase");
  params.addRequiredCoupledVar("phase_beta", "phase field for beta phase");
  params.addCoupledVar("phase_gamma",0.0, "phase field for gamma phase");
  params.addCoupledVar("phase_delta",0.0, "phase field for delta phase");
  params.addCoupledVar("phase_epsilon",0.0, "phase field for gamma phase");
  params.addParam<MaterialPropertyName>("xB_gamma", 0.0, "Phase composition in the gamma-phase");
  params.addParam<MaterialPropertyName>("xB_delta", 0.0, "Phase composition in the delta-phase");
  params.addParam<MaterialPropertyName>("xB_epsilon", 0.0, "Phase composition in the epsilon-phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_gamma", 0.0, "Thermodynamic factor in gamma phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_delta", 0.0, "Thermodynamic factor in delta phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_epsilon", 0.0, "Thermodynamic factor in epsilon phase");
  params.addParam<MaterialPropertyName>("h_gamma", 0.0, "Gamma");
  params.addParam<MaterialPropertyName>("h_delta", 0.0, "Delta");
  params.addParam<MaterialPropertyName>("h_epsilon", 0.0, "Epsilon");
  params.addRequiredCoupledVar("xB", "Component B mole fraction");   
  return params;
}

MultiPhaseConstraintMu::MultiPhaseConstraintMu(const InputParameters & parameters)
  : Kernel(parameters), 
    //alpha phase
    _phase_alpha(coupledValue("phase_alpha")),
    _phase_alpha_var(coupled("phase_alpha")),
    //beta phase
    _phase_beta(coupledValue("phase_beta")),
    _phase_beta_var(coupled("phase_beta")),
    //gamma phase
    _phase_gamma(coupledValue("phase_gamma")),
    _phase_gamma_var(coupled("phase_gamma")),
    //delta phase
    _phase_delta(coupledValue("phase_delta")),
    _phase_delta_var(coupled("phase_delta")),
     //epsilon phase
    _phase_epsilon(coupledValue("phase_epsilon")),
    _phase_epsilon_var(coupled("phase_epsilon")),
    //interpolation function
    _h_alpha(getMaterialProperty<Real>("h_alpha")),
    _h_beta(getMaterialProperty<Real>("h_beta")),
    _h_gamma(getMaterialProperty<Real>("h_gamma")),
    _h_delta(getMaterialProperty<Real>("h_delta")),
    _h_epsilon(getMaterialProperty<Real>("h_epsilon")),
    // Note: Only non-diagonal components of the interpolation marix
    // is required
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
    _dhbeta_dphigamma(getMaterialProperty<Real>("dhbeta_dphigamma")),
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
     //composition
    _xB(coupledValue("xB")),
    _xB_var(coupled("xB")),
    //Material property composition
    _xB_alpha(getMaterialProperty<Real>("xB_alpha")),
    _xB_beta(getMaterialProperty<Real>("xB_beta")),
    _xB_gamma(getMaterialProperty<Real>("xB_gamma")),
    _xB_delta(getMaterialProperty<Real>("xB_delta")),
    _xB_epsilon(getMaterialProperty<Real>("xB_epsilon")),
    //inverse of the thermodynamic factor
    _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
    _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
    _inv_B_tf_gamma(getMaterialProperty<Real>("inv_B_tf_gamma")),
    _inv_B_tf_delta(getMaterialProperty<Real>("inv_B_tf_delta")),
    _inv_B_tf_epsilon(getMaterialProperty<Real>("inv_B_tf_epsilon"))
{
}

Real
MultiPhaseConstraintMu::computeQpResidual()
{ 
  // Right-hand side of the equation
  
 const Real _weighted_sum = (_xB_alpha[_qp] * _h_alpha[_qp] 
                           + _xB_beta[_qp] * _h_beta[_qp] 
                           + _xB_gamma[_qp]* _h_gamma[_qp]
                           + _xB_delta[_qp]* _h_delta[_qp]
                           + _xB_epsilon[_qp]* _h_epsilon[_qp]);

  // w_i*(sum(c_{theta}* h_theta} - c)= 0
  return (_test[_i][_qp] * (_weighted_sum - _xB[_qp]));
}

Real
MultiPhaseConstraintMu::computeQpJacobian()
{
  return (_test[_i][_qp] * ( (_h_alpha[_qp]  *_inv_B_tf_alpha[_qp]) 
                           + (_h_beta[_qp]   *_inv_B_tf_beta[_qp])
                           + (_h_gamma[_qp]  *_inv_B_tf_gamma[_qp])
                           + (_h_delta[_qp]  *_inv_B_tf_delta[_qp]) 
                           + (_h_epsilon[_qp]*_inv_B_tf_epsilon[_qp])) *_phi[_j][_qp]);
}

Real
MultiPhaseConstraintMu::computeQpOffDiagJacobian(unsigned int jvar)
{  
  if (jvar == _phase_alpha_var)
  {  
    return (_test[_i][_qp] * ((_xB_beta[_qp]  - _xB_alpha[_qp])* _dhbeta_dphialpha[_qp] 
                             +(_xB_gamma[_qp] - _xB_alpha[_qp])* _dhgamma_dphialpha[_qp]
                             +(_xB_delta[_qp] - _xB_alpha[_qp])* _dhdelta_dphialpha[_qp]
                             +(_xB_epsilon[_qp] - _xB_alpha[_qp])* _dhepsilon_dphialpha[_qp]) * _phi[_j][_qp]);
  }     
  else if (jvar == _phase_beta_var)
  {
    return (_test[_i][_qp] * ((_xB_alpha[_qp] - _xB_beta[_qp])* _dhalpha_dphibeta[_qp] 
                             +(_xB_gamma[_qp] - _xB_beta[_qp])* _dhgamma_dphibeta[_qp] 
                             +(_xB_delta[_qp] - _xB_beta[_qp])* _dhdelta_dphibeta[_qp]
                             +(_xB_epsilon[_qp] - _xB_beta[_qp])*_dhepsilon_dphibeta[_qp])* _phi[_j][_qp]);
  }  
  else if (jvar == _phase_gamma_var)
  {
    return (_test[_i][_qp] *((_xB_alpha[_qp] - _xB_gamma[_qp])* _dhalpha_dphigamma[_qp] 
                            +(_xB_beta[_qp]  - _xB_gamma[_qp])* _dhbeta_dphigamma[_qp]
                            +(_xB_delta[_qp] - _xB_gamma[_qp])* _dhdelta_dphigamma[_qp]
                            +(_xB_epsilon[_qp] - _xB_gamma[_qp])* _dhepsilon_dphigamma[_qp])* _phi[_j][_qp]);
  } 
  else if (jvar == _phase_delta_var)
  {
    return (_test[_i][_qp] *((_xB_alpha[_qp] - _xB_delta[_qp])* _dhalpha_dphidelta[_qp] 
                            +(_xB_beta[_qp]  - _xB_delta[_qp])* _dhbeta_dphidelta[_qp]
                            +(_xB_gamma[_qp] - _xB_delta[_qp])* _dhgamma_dphidelta[_qp]
                            +(_xB_epsilon[_qp] - _xB_delta[_qp])* _dhepsilon_dphidelta[_qp])* _phi[_j][_qp]);
  }
  else if (jvar == _phase_epsilon_var)
  {
    return (_test[_i][_qp] *((_xB_alpha[_qp] - _xB_epsilon[_qp])* _dhalpha_dphiepsilon[_qp] 
                            +(_xB_beta[_qp]  - _xB_epsilon[_qp])* _dhbeta_dphiepsilon[_qp]
                            +(_xB_gamma[_qp] - _xB_epsilon[_qp])* _dhgamma_dphiepsilon[_qp]
                            +(_xB_delta[_qp] - _xB_epsilon[_qp])* _dhdelta_dphiepsilon[_qp])* _phi[_j][_qp]);
  }   
  else if (jvar == _xB_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp]);
  }  
  else  
    return 0.0;
}
