//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TCPhaseConstraintMuC.h"

registerMooseObject("gibbsApp", TCPhaseConstraintMuC);

template <>
InputParameters
validParams<TCPhaseConstraintMuC>()
{
  InputParameters params = validParams<MultiPhaseBase>();
  params.addClassDescription("Eqn: (1-h(eta))*xC_alpha + h(eta)*xC_beta - xC= 0."
                             "non-linear variable of this kernel is mu_C.");
  params.addRequiredCoupledVar("xC", "Overall Component C mole fraction");
  params.addParam<MaterialPropertyName>("xC_gamma", 0.0 , "Mole fraction of C in gamma phase");
  params.addParam<MaterialPropertyName>("xC_delta", 0.0 , "Mole fraction of C in delta phase");
  params.addParam<MaterialPropertyName>("xC_epsilon", 0.0, "Mole fraction of C in epsilon phase");
  //BC thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_BC_tf_gamma", 0.0, "Thermodynamic factor BC in gamma phase");
  params.addParam<MaterialPropertyName>("inv_BC_tf_delta", 0.0, "Thermodynamic factor BC in delta phase");
  params.addParam<MaterialPropertyName>("inv_BC_tf_epsilon", 0.0, "Thermodynamic factor BC in epsilon phase");
  //C thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_C_tf_gamma", 0.0, "Thermodynamic factor C in gamma phase");
  params.addParam<MaterialPropertyName>("inv_C_tf_delta", 0.0, "Thermodynamic factor C in delta phase");
  params.addParam<MaterialPropertyName>("inv_C_tf_epsilon", 0.0, "Thermodynamic factor C in epsilon phase");
  //Coupled variable
  params.addRequiredCoupledVar("B_diff_pot", "Component B diffusion potential");
  return params;
}

TCPhaseConstraintMuC::TCPhaseConstraintMuC(const InputParameters & parameters)
  : MultiPhaseBase(parameters),
    _xC(coupledValue("xC")),
    _xC_var(coupled("xC")),
     //Material property composition
    _xC_alpha(getMaterialProperty<Real>("xC_alpha")),
    _xC_beta(getMaterialProperty<Real>("xC_beta")),
    _xC_gamma(getMaterialProperty<Real>("xC_gamma")),
    _xC_delta(getMaterialProperty<Real>("xC_delta")),
    _xC_epsilon(getMaterialProperty<Real>("xC_epsilon")),
     //coeff BC of the susceptibility matrix
    _inv_BC_tf_alpha(getMaterialProperty<Real>("inv_BC_tf_alpha")),
    _inv_BC_tf_beta(getMaterialProperty<Real>("inv_BC_tf_beta")),
    _inv_BC_tf_gamma(getMaterialProperty<Real>("inv_BC_tf_gamma")),
    _inv_BC_tf_delta(getMaterialProperty<Real>("inv_BC_tf_delta")),
    _inv_BC_tf_epsilon(getMaterialProperty<Real>("inv_BC_tf_epsilon")),
    //coeff C of the susceptibility matrix
    _inv_C_tf_alpha(getMaterialProperty<Real>("inv_C_tf_alpha")),
    _inv_C_tf_beta(getMaterialProperty<Real>("inv_C_tf_beta")),
    _inv_C_tf_gamma(getMaterialProperty<Real>("inv_C_tf_gamma")),
    _inv_C_tf_delta(getMaterialProperty<Real>("inv_C_tf_delta")),
    _inv_C_tf_epsilon(getMaterialProperty<Real>("inv_C_tf_epsilon")),
    //For a ternary alloy A-B-C
    _B_diff_pot(coupledValue("B_diff_pot")),
    _B_diff_pot_var(coupled("B_diff_pot"))
{
}

Real
TCPhaseConstraintMuC::chi_BC() const
{
   return (_h_alpha[_qp]  * _inv_BC_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_BC_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_BC_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_BC_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_BC_tf_epsilon[_qp]);
}

Real
TCPhaseConstraintMuC::chi_CC() const
{
   return (_h_alpha[_qp]  * _inv_C_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_C_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_C_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_C_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_C_tf_epsilon[_qp]); 
}

Real
TCPhaseConstraintMuC::computeQpResidual()
{
  const Real _weighted_sum = (_xC_alpha[_qp] * _h_alpha[_qp] 
                            + _xC_beta[_qp]  * _h_beta[_qp] 
                            + _xC_gamma[_qp] * _h_gamma[_qp]
                            + _xC_delta[_qp] * _h_delta[_qp]
                            + _xC_epsilon[_qp]* _h_epsilon[_qp]);

  // w_i*(sum(c_{theta}* h_theta} - c)= 0
  return (_test[_i][_qp] * (_weighted_sum - _xC[_qp]));
}

Real
TCPhaseConstraintMuC::computeQpJacobian()
{
 return (_test[_i][_qp] * (TCPhaseConstraintMuC::chi_CC()) *_phi[_j][_qp]);
}

Real
TCPhaseConstraintMuC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_alpha_var)
  {  
    return (_test[_i][_qp] * ((_xC_beta[_qp]  - _xC_alpha[_qp])* _dhbeta_dphialpha[_qp] 
                             +(_xC_gamma[_qp] - _xC_alpha[_qp])* _dhgamma_dphialpha[_qp]
                             +(_xC_delta[_qp] - _xC_alpha[_qp])* _dhdelta_dphialpha[_qp]
                             +(_xC_epsilon[_qp]- _xC_alpha[_qp])* _dhepsilon_dphialpha[_qp]) * _phi[_j][_qp]);
  }     
  else if (jvar == _phase_beta_var)
  {
    return (_test[_i][_qp] * ((_xC_alpha[_qp] - _xC_beta[_qp])* _dhalpha_dphibeta[_qp] 
                             +(_xC_gamma[_qp] - _xC_beta[_qp])* _dhgamma_dphibeta[_qp] 
                             +(_xC_delta[_qp] - _xC_beta[_qp])* _dhdelta_dphibeta[_qp]
                             +(_xC_epsilon[_qp] - _xC_beta[_qp])*_dhepsilon_dphibeta[_qp])* _phi[_j][_qp]);
  }  
  else if (jvar == _phase_gamma_var)
  {
    return (_test[_i][_qp] *((_xC_alpha[_qp] - _xC_gamma[_qp])* _dhalpha_dphigamma[_qp] 
                            +(_xC_beta[_qp]  - _xC_gamma[_qp])* _dhbeta_dphigamma[_qp]
                            +(_xC_delta[_qp] - _xC_gamma[_qp])* _dhdelta_dphigamma[_qp]
                            +(_xC_epsilon[_qp] - _xC_gamma[_qp])* _dhepsilon_dphigamma[_qp])* _phi[_j][_qp]);
  } 
  else if (jvar == _phase_delta_var)
  {
    return (_test[_i][_qp] *((_xC_alpha[_qp] - _xC_delta[_qp])* _dhalpha_dphidelta[_qp] 
                            +(_xC_beta[_qp]  - _xC_delta[_qp])* _dhbeta_dphidelta[_qp]
                            +(_xC_gamma[_qp] - _xC_delta[_qp])* _dhgamma_dphidelta[_qp]
                            +(_xC_epsilon[_qp] - _xC_delta[_qp])* _dhepsilon_dphidelta[_qp])* _phi[_j][_qp]);
  }
  else if (jvar == _phase_epsilon_var)
  {
    return (_test[_i][_qp] *((_xC_alpha[_qp] - _xC_epsilon[_qp])* _dhalpha_dphiepsilon[_qp] 
                            +(_xC_beta[_qp]  - _xC_epsilon[_qp])* _dhbeta_dphiepsilon[_qp]
                            +(_xC_gamma[_qp] - _xC_epsilon[_qp])* _dhgamma_dphiepsilon[_qp]
                            +(_xC_delta[_qp] - _xC_epsilon[_qp])* _dhdelta_dphiepsilon[_qp])* _phi[_j][_qp]);
  }   
  else if (jvar == _xC_var)
  {
    return -(_test[_i][_qp] * _phi[_j][_qp]);
  }  
  else if (jvar == _B_diff_pot_var)
  {
    return (_test[_i][_qp] * TCPhaseConstraintMuC::chi_BC() *_phi[_j][_qp]);
  }
  else
    return 0.0;
}
