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

#include "MultiPhaseBase.h"

registerMooseObject("gibbsApp", MultiPhaseBase);

template <>
InputParameters
validParams<MultiPhaseBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("phase_alpha", "phase 1");
  params.addRequiredCoupledVar("phase_beta", "Phase 2");
  params.addCoupledVar("phase_gamma", 0.0, "Phase 3");
  params.addCoupledVar("phase_delta", 0.0, "Phase 4");
  params.addCoupledVar("phase_epsilon", 0.0, "Phase 5");    
  params.addRequiredParam<MaterialPropertyName>("h_alpha", "interpolation");
  params.addRequiredParam<MaterialPropertyName>("h_beta", "interpolation");
  params.addParam<MaterialPropertyName>("h_gamma",0.0, "interpolation");
  params.addParam<MaterialPropertyName>("h_delta",0.0, "interpolation");
  params.addParam<MaterialPropertyName>("h_epsilon",0.0, "interpolation");
  return params; 
}

MultiPhaseBase::MultiPhaseBase(const InputParameters & parameters)
  : Kernel(parameters),
  //can be extended to N...phases
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
  _dhdelta_dphiepsilon(getMaterialProperty<Real>("dhdelta_dphiepsilon"))
{
}

Real
MultiPhaseBase::computeQpResidual(){
  return 0;
}

Real
MultiPhaseBase::computeQpJacobian(){
  return 0;
}

Real
MultiPhaseBase::computeQpOffDiagJacobian(unsigned int /*jvar*/){
  return 0;
}
