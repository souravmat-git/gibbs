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

#include "MCPhaseConstraintMuB.h"

registerMooseObject("gibbsApp", MCPhaseConstraintMuB);

template <>
InputParameters
validParams<MCPhaseConstraintMuB>()
{
  InputParameters params = validParams<TCPhaseConstraintMuB>();
  params.addClassDescription("Eqn: (1-h(eta))*xB_alpha + h(eta)*xB_beta - xB = 0."
                             "non-linear variable of this kernel is muB");
  //BD thermodynamic factor: Default it assumes two phases
  params.addParam<MaterialPropertyName>("inv_BD_tf_gamma", 0.0, "Thermodynamic factor BC in gamma phase");
  params.addParam<MaterialPropertyName>("inv_BD_tf_delta", 0.0, "Thermodynamic factor BC in delta phase");
  params.addParam<MaterialPropertyName>("inv_BD_tf_epsilon", 0.0, "Thermodynamic factor BC in epsilon phase");
  params.addRequiredCoupledVar("D_diff_pot", "Component D diffusion potential");
  return params;
}

MCPhaseConstraintMuB::MCPhaseConstraintMuB(const InputParameters & parameters)
  : TCPhaseConstraintMuB(parameters),
    //coeff BD of the susceptibility matrix
    _inv_BD_tf_alpha(getMaterialProperty<Real>("inv_BD_tf_alpha")),
    _inv_BD_tf_beta(getMaterialProperty<Real>("inv_BD_tf_beta")),
    _inv_BD_tf_gamma(getMaterialProperty<Real>("inv_BD_tf_gamma")),
    _inv_BD_tf_delta(getMaterialProperty<Real>("inv_BD_tf_delta")),
    _inv_BD_tf_epsilon(getMaterialProperty<Real>("inv_BD_tf_epsilon")),
    //For a quaternary alloy A-B-C-D
    _D_diff_pot(coupledValue("D_diff_pot")),
    _D_diff_pot_var(coupled("D_diff_pot"))
{
}

Real
MCPhaseConstraintMuB::chi_BD() const
{
   return (_h_alpha[_qp]  * _inv_BD_tf_alpha[_qp] 
         + _h_beta[_qp]   * _inv_BD_tf_beta[_qp]
         + _h_gamma[_qp]  * _inv_BD_tf_gamma[_qp]
         + _h_delta[_qp]  * _inv_BD_tf_delta[_qp]
         + _h_epsilon[_qp]* _inv_BD_tf_epsilon[_qp]);
 
}

Real
MCPhaseConstraintMuB::computeQpResidual()
{                         
  // The kernel operates on the variable: Diffusion potential of comp B  
   return (TCPhaseConstraintMuB::computeQpResidual());
}

Real
MCPhaseConstraintMuB::computeQpJacobian()
{
  return (TCPhaseConstraintMuB::computeQpJacobian());
}

Real
MCPhaseConstraintMuB::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_alpha_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_beta_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_gamma_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_delta_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_epsilon_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _xB_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _C_diff_pot_var)
  {
    return (TCPhaseConstraintMuB::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _D_diff_pot_var)
  {
    return (_test[_i][_qp] * MCPhaseConstraintMuB::chi_BD() *_phi[_j][_qp]);  
  }
  else
    return 0.0;
}
