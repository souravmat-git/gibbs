//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiPhaseChiTimeDerivative.h"

registerMooseObject("gibbsApp", MultiPhaseChiTimeDerivative);

template <>
InputParameters
validParams<MultiPhaseChiTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  //The non-variable that this kernel operates on is c_{gamma}                       
  params.addRequiredCoupledVar("phase_alpha", "phase field for alpha phase");
  params.addRequiredCoupledVar("phase_beta", "phase field for beta phase");
  params.addCoupledVar("phase_gamma",0.0, "phase field for gamma phase");
  params.addParam<MaterialPropertyName>("inv_B_tf_gamma", 0.0, "Phase composition in the gamma-phase");
  params.addParam<MaterialPropertyName>("inv_B_td_gamma", 0.0, "Thermodynamic factor in gamma phase");
  params.addParam<MaterialPropertyName>("h_gamma", 0.0, "Gamma");
  return params;
}

MultiPhaseChiTimeDerivative::MultiPhaseChiTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters),
  //alpha phase
  _phase_alpha(coupledValue("phase_alpha")),
  _phase_alpha_var(coupled("phase_alpha")),
  //beta phase
  _phase_beta(coupledValue("phase_beta")),
  _phase_beta_var(coupled("phase_beta")),
  //gamma phase
  _phase_gamma(coupledValue("phase_gamma")),
  _phase_gamma_var(coupled("phase_gamma")),
  //inverse of the thermodynamic factor
  _inv_B_tf_alpha(getMaterialProperty<Real>("inv_B_tf_alpha")),
  _inv_B_tf_beta(getMaterialProperty<Real>("inv_B_tf_beta")),
  _inv_B_tf_gamma(getMaterialProperty<Real>("inv_B_tf_gamma")),
  //third derivative
  _inv_B_td_alpha(getMaterialProperty<Real>("inv_B_td_alpha")),
  _inv_B_td_beta(getMaterialProperty<Real>("inv_B_td_beta")),
  _inv_B_td_gamma(getMaterialProperty<Real>("inv_B_td_gamma")),
   //interpolation function
  _h_alpha(getMaterialProperty<Real>("h_alpha")),
  _h_beta(getMaterialProperty<Real>("h_beta")),
  _h_gamma(getMaterialProperty<Real>("h_gamma")), 
  // Note: Only non-diagonal components of the interpolation marix
  _dhbeta_dphialpha(getMaterialProperty<Real>("dhbeta_dphialpha")),
  _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),  
    
   // for coupled variable phase_beta
   _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
   _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")), 
    
   // for coupled variable phase_gamma    
   _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
   _dhbeta_dphigamma(getMaterialProperty<Real>("dhbeta_dphigamma"))
{
}

Real
MultiPhaseChiTimeDerivative::interpolated_chi() const
{
  return ( _h_alpha[_qp]*_inv_B_tf_alpha[_qp] + _h_beta[_qp]*_inv_B_tf_beta[_qp]
         + _h_gamma[_qp]*_inv_B_tf_gamma[_qp]);
}

Real
MultiPhaseChiTimeDerivative::third_deriv() const
{
 //this is rate of change of the curvature of the free energy curve with composition
  return (_h_beta[_qp] * _inv_B_td_beta[_qp] + _h_alpha[_qp] * _inv_B_td_alpha[_qp]);
}

Real
MultiPhaseChiTimeDerivative::computeQpResidual()
{
  return (TimeDerivative::computeQpResidual() * MultiPhaseChiTimeDerivative::interpolated_chi());
}

Real
MultiPhaseChiTimeDerivative::computeQpJacobian()
{
  return (TimeDerivative::computeQpJacobian() * MultiPhaseChiTimeDerivative::interpolated_chi() +
         TimeDerivative::computeQpResidual() * MultiPhaseChiTimeDerivative::third_deriv() * _phi[_j][_qp]);
}

Real
MultiPhaseChiTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_alpha_var)
  {  
    return (TimeDerivative::computeQpOffDiagJacobian(jvar) * ( _inv_B_tf_beta[_qp]  - _inv_B_tf_alpha[_qp])* _dhbeta_dphialpha[_qp] 
                                                        +(_inv_B_tf_gamma[_qp] - _inv_B_tf_alpha[_qp])* _dhgamma_dphialpha[_qp]);
  }     
  else if (jvar == _phase_beta_var)
  {
    return (TimeDerivative::computeQpOffDiagJacobian(jvar) * ( _inv_B_tf_alpha[_qp] - _inv_B_tf_beta[_qp])* _dhalpha_dphibeta[_qp] 
                                                        +(_inv_B_tf_gamma[_qp]  - _inv_B_tf_beta[_qp])* _dhgamma_dphibeta[_qp]);
  }  
  else if (jvar == _phase_gamma_var)
  {
    return (TimeDerivative::computeQpOffDiagJacobian(jvar) * (_inv_B_tf_alpha[_qp] - _inv_B_tf_gamma[_qp])* _dhalpha_dphigamma[_qp] 
                                                        +( _inv_B_tf_beta[_qp] - _inv_B_tf_gamma[_qp])* _dhbeta_dphigamma[_qp]);
  }
  else
    return 0.0;  
}
