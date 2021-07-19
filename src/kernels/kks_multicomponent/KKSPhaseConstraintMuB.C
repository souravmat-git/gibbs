//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This code modifies the KKSPhaseConcentration in MOOSE
//* by including the hand coded interpolation function
//* instead of using the Material property- S.Chatterjee

#include "KKSPhaseConstraintMuB.h"

registerMooseObject("gibbsApp", KKSPhaseConstraintMuB);

template <>
InputParameters
validParams<KKSPhaseConstraintMuB>()
{
  InputParameters params = validParams<MultiCompMultiPhaseBase>();
  params.addClassDescription("Eqn: (1-h(eta))*xB_alpha + h(eta)*xB_beta - xB = 0."
                             "non-linear variable of this kernel is muC");
  params.addRequiredCoupledVar("xB", "Component B mole fraction");
  params.addParam<MaterialPropertyName>("xB_gamma", 0.0, "Mole fraction of B in gamma phase");
  params.addParam<MaterialPropertyName>("xB_delta", 0.0, "Mole fraction of B in delta phase");
  params.addParam<MaterialPropertyName>("xB_epsilon", 0.0, "Mole fraction of B in epsilon phase");
  params.addRequiredCoupledVar("C_diff_pot", "Component C diffusion potential");
  params.addCoupledVar("D_diff_pot", 0.0, "Component D diffusion potential");
  return params;
}

KKSPhaseConstraintMuB::KKSPhaseConstraintMuB(const InputParameters & parameters)
  :  MultiCompMultiPhaseBase(parameters),
    _xB(coupledValue("xB")),
    _xB_var(coupled("xB")),
    //Material property: Mole fraction 
    _xB_alpha(getMaterialProperty<Real>("xB_alpha")),
    _xB_beta(getMaterialProperty<Real>("xB_beta")),
    _xB_gamma(getMaterialProperty<Real>("xB_gamma")),
    _xB_delta(getMaterialProperty<Real>("xB_delta")),
    _xB_epsilon(getMaterialProperty<Real>("xB_epsilon")),
    //For a ternary alloy A-B-C
    _C_diff_pot(coupledValue("C_diff_pot")),
    _C_diff_pot_var(coupled("C_diff_pot")),
    //For a quaternary alloy A-B-C-D
    _D_diff_pot(coupledValue("D_diff_pot")),
    _D_diff_pot_var(coupled("D_diff_pot"))
{
}

Real
KKSPhaseConstraintMuB::computeQpResidual()
{
   const Real _weighted_sum = (_xB_alpha[_qp] * _h_alpha[_qp] 
                            + _xB_beta[_qp] * _h_beta[_qp] 
                            + _xB_gamma[_qp]* _h_gamma[_qp]
                            + _xB_delta[_qp]* _h_delta[_qp]
                          + _xB_epsilon[_qp]* _h_epsilon[_qp]);
                          
  // The kernel operates on the variable: Diffusion potential of comp B  
   return (_test[_i][_qp] * (_weighted_sum - _xB[_qp]));
}

Real
KKSPhaseConstraintMuB::computeQpJacobian()
{
  return (_test[_i][_qp] * MultiCompMultiPhaseBase::chi_BB() *_phi[_j][_qp]);
}

Real
KKSPhaseConstraintMuB::computeQpOffDiagJacobian(unsigned int jvar)
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
    return (_test[_i][_qp] *((_xB_alpha[_qp] - _xB_delta[_qp])* _dhalpha_dphidelta[_qp] 
                            +(_xB_beta[_qp]  - _xB_delta[_qp])* _dhbeta_dphidelta[_qp]
                            +(_xB_gamma[_qp] - _xB_delta[_qp])* _dhgamma_dphidelta[_qp]
                            +(_xB_epsilon[_qp] - _xB_delta[_qp])* _dhepsilon_dphidelta[_qp])* _phi[_j][_qp]);
  }
  else if (jvar == _xB_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp]);
  }
  else if (jvar == _C_diff_pot_var)
  {
    return (_test[_i][_qp] * MultiCompMultiPhaseBase::chi_BC() *_phi[_j][_qp]);
  }
  else if (jvar == _D_diff_pot_var)
  {
    return (_test[_i][_qp] * MultiCompMultiPhaseBase::chi_BD() *_phi[_j][_qp]);  
  }
  else
    return 0.0;
}
