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

#include "MCPhaseConstraintMuC.h"

registerMooseObject("gibbsApp", MCPhaseConstraintMuC);

template <>
InputParameters
validParams<MCPhaseConstraintMuC>()
{
  InputParameters params = validParams<TCPhaseConstraintMuC>();
  params.addClassDescription("Eqn: (1-h(eta))*xC_alpha + h(eta)*xC_beta - xC = 0."
                             "non-linear variable of this kernel is muC");
  //CD thermodynamic factor
  params.addParam<MaterialPropertyName>("inv_CD_tf_gamma", 0.0, "Thermodynamic factor CD in gamma phase");
  params.addParam<MaterialPropertyName>("inv_CD_tf_delta", 0.0, "Thermodynamic factor CD in delta phase");
  params.addParam<MaterialPropertyName>("inv_CD_tf_epsilon", 0.0, "Thermodynamic factor CD in epsilon phase");
  params.addRequiredCoupledVar("D_diff_pot", "Component D diffusion potential");
  return params;
}

MCPhaseConstraintMuC::MCPhaseConstraintMuC(const InputParameters & parameters)
  : TCPhaseConstraintMuC(parameters),
    //coeff CD of the susceptibility matrix
    _inv_CD_tf_alpha(getMaterialProperty<Real>("inv_CD_tf_alpha")),
    _inv_CD_tf_beta(getMaterialProperty<Real>("inv_CD_tf_beta")),
    _inv_CD_tf_gamma(getMaterialProperty<Real>("inv_CD_tf_gamma")),
    _inv_CD_tf_delta(getMaterialProperty<Real>("inv_CD_tf_delta")),
    _inv_CD_tf_epsilon(getMaterialProperty<Real>("inv_CD_tf_epsilon")),
    //For a quaternary alloy A-B-C-D
    _D_diff_pot(coupledValue("D_diff_pot")),
    _D_diff_pot_var(coupled("D_diff_pot"))
{
}

Real
MCPhaseConstraintMuC::chi_CD() const
{
   return (_h_alpha[_qp]  * _inv_CD_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_CD_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_CD_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_CD_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_CD_tf_epsilon[_qp]);
}

Real
MCPhaseConstraintMuC::computeQpResidual()
{                         
  // The kernel operates on the variable: Diffusion potential of comp C  
   return (TCPhaseConstraintMuC::computeQpResidual());
}

Real
MCPhaseConstraintMuC::computeQpJacobian()
{
  return (TCPhaseConstraintMuC::computeQpJacobian());
}

Real
MCPhaseConstraintMuC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_alpha_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_beta_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_gamma_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_delta_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_epsilon_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _xC_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _B_diff_pot_var)
  {
    return (TCPhaseConstraintMuC::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _D_diff_pot_var)
  {
    return (_test[_i][_qp] * MCPhaseConstraintMuC::chi_CD() *_phi[_j][_qp]);  
  }
  else
    return 0.0;
}
