//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//*

#include "KKSPhaseConstraintMuD.h"

registerMooseObject("gibbsApp", KKSPhaseConstraintMuD);

template <>
InputParameters
validParams<KKSPhaseConstraintMuD>()
{
  InputParameters params = validParams<MultiCompMultiPhaseBase>();
  params.addClassDescription("Eqn: (1-h(eta))*xD_alpha + h(eta)*xD_beta - xD= 0."
                             "non-linear variable of this kernel is mu_D.");
  params.addRequiredCoupledVar("xD", "Component D mole fraction");
  params.addRequiredCoupledVar("C_diff_pot", "Component C diffusion potential");
  params.addRequiredCoupledVar("B_diff_pot", "Component B diffusion potential");
  params.addParam<MaterialPropertyName>("xD_gamma", 0.0 , "Mole fraction of D in gamma phase");
  params.addParam<MaterialPropertyName>("xD_delta", 0.0 , "Mole fraction of D in delta phase");
  params.addParam<MaterialPropertyName>("xD_epsilon", 0.0, "Mole fraction of D in epsilon phase");
  return params;
}

KKSPhaseConstraintMuD::KKSPhaseConstraintMuD(const InputParameters & parameters)
  :  MultiCompMultiPhaseBase(parameters),
    _xD(coupledValue("xD")),
    _xD_var(coupled("xD")),
     //Material property required by this kernel
    _xD_alpha(getMaterialProperty<Real>("xD_alpha")),
    _xD_beta(getMaterialProperty<Real>("xD_beta")),
    _xD_gamma(getMaterialProperty<Real>("xD_gamma")),
    _xD_delta(getMaterialProperty<Real>("xD_delta")),
    _xD_epsilon(getMaterialProperty<Real>("xD_epsilon")),
    //Coupled variables
    _C_diff_pot(coupledValue("C_diff_pot")),
    _C_diff_pot_var(coupled("C_diff_pot")),
    _B_diff_pot(coupledValue("B_diff_pot")),
    _B_diff_pot_var(coupled("B_diff_pot"))
{
}

Real
KKSPhaseConstraintMuD::computeQpResidual()
{
   const Real _weighted_sum = (_xD_alpha[_qp] * _h_alpha[_qp] 
                            + _xD_beta[_qp]   * _h_beta[_qp] 
                            + _xD_gamma[_qp]  * _h_gamma[_qp]
                            + _xD_delta[_qp]  * _h_delta[_qp]
                            + _xD_epsilon[_qp]* _h_epsilon[_qp]);
  // The kernel operates on the variable: Diffusion potential of comp D
  return (_test[_i][_qp] * (_weighted_sum - _xD[_qp]));     
}

Real
KKSPhaseConstraintMuD::computeQpJacobian()
{
  return (_test[_i][_qp] * MultiCompMultiPhaseBase::chi_DD() *_phi[_j][_qp]);  
}

Real
KKSPhaseConstraintMuD::computeQpOffDiagJacobian(unsigned int jvar)
{  

  if (jvar == _phase_alpha_var)
  {  
    return (_test[_i][_qp] * ((_xD_beta[_qp]  - _xD_alpha[_qp])* _dhbeta_dphialpha[_qp] 
                             +(_xD_gamma[_qp] - _xD_alpha[_qp])* _dhgamma_dphialpha[_qp]
                             +(_xD_delta[_qp] - _xD_alpha[_qp])* _dhdelta_dphialpha[_qp]
                             +(_xD_epsilon[_qp]- _xD_alpha[_qp])* _dhepsilon_dphialpha[_qp]) * _phi[_j][_qp]);
  }     
  else if (jvar == _phase_beta_var)
  {
    return (_test[_i][_qp] * ((_xD_alpha[_qp] - _xD_beta[_qp])* _dhalpha_dphibeta[_qp] 
                             +(_xD_gamma[_qp] - _xD_beta[_qp])* _dhgamma_dphibeta[_qp] 
                             +(_xD_delta[_qp] - _xD_beta[_qp])* _dhdelta_dphibeta[_qp]
                             +(_xD_epsilon[_qp] - _xD_beta[_qp])*_dhepsilon_dphibeta[_qp])* _phi[_j][_qp]);
  }  
  else if (jvar == _phase_gamma_var)
  {
    return (_test[_i][_qp] *((_xD_alpha[_qp] - _xD_gamma[_qp])*  _dhalpha_dphigamma[_qp] 
                            +(_xD_beta[_qp]  - _xD_gamma[_qp])*  _dhbeta_dphigamma[_qp]
                            +(_xD_delta[_qp] - _xD_gamma[_qp])*  _dhdelta_dphigamma[_qp]
                            +(_xD_epsilon[_qp] - _xD_gamma[_qp])*_dhepsilon_dphigamma[_qp])* _phi[_j][_qp]);
  } 
  else if (jvar == _phase_delta_var)
  {
    return (_test[_i][_qp] *((_xD_alpha[_qp] - _xD_delta[_qp])*   _dhalpha_dphidelta[_qp] 
                            +(_xD_beta[_qp]  - _xD_delta[_qp])*   _dhbeta_dphidelta[_qp]
                            +(_xD_gamma[_qp] - _xD_delta[_qp])*   _dhgamma_dphidelta[_qp]
                            +(_xD_epsilon[_qp] - _xD_delta[_qp])* _dhepsilon_dphidelta[_qp])* _phi[_j][_qp]);
  }
  else if (jvar == _phase_epsilon_var)
  {
    return (_test[_i][_qp] *((_xD_alpha[_qp] - _xD_epsilon[_qp])* _dhalpha_dphiepsilon[_qp] 
                            +(_xD_beta[_qp]  - _xD_epsilon[_qp])* _dhbeta_dphiepsilon[_qp]
                            +(_xD_gamma[_qp] - _xD_epsilon[_qp])* _dhgamma_dphiepsilon[_qp]
                            +(_xD_delta[_qp] - _xD_epsilon[_qp])* _dhdelta_dphiepsilon[_qp])* _phi[_j][_qp]);
  }   
  else if (jvar == _xD_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp]);
  }
  else if (jvar == _C_diff_pot_var)
  {
     return (_test[_i][_qp] * MultiCompMultiPhaseBase::chi_CD() *_phi[_j][_qp]);
  }
  else if (jvar == _B_diff_pot_var)
  {
    return (_test[_i][_qp] * MultiCompMultiPhaseBase::chi_BD()  *_phi[_j][_qp]);    
  }
  else
    return 0.0;
}
