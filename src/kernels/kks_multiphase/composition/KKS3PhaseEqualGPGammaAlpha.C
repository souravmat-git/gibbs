//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*

//*This Kernel implements the bulk part of the
//*Allen Cahn equation which tells us that: at
//*equilibrium the grandpotentials must be equal
//*In this code both the free energy can be 
//*extracted from the input file

#include "KKS3PhaseEqualGPGammaAlpha.h"

registerMooseObject("gibbsApp",KKS3PhaseEqualGPGammaAlpha);

template <>
InputParameters
validParams<KKS3PhaseEqualGPGammaAlpha>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: df/dphi_{alpha} = dh_{gamma}/dphi{_alpha}*(omega_gamma - omega_alpha)."
                             "omega_gamma, omega_alpha is the grand potential"
                             "This kernel operates on phi_alpha.");
  params.addRequiredCoupledVar("phase_comp_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("phase_comp_gamma", "Phase concentration in beta phase");
  params.addCoupledVar("phase_beta",0.0, "Beta phase field variable");
  params.addRequiredCoupledVar("phase_gamma", "Gamma phase field variable");
  params.addRequiredParam<MaterialPropertyName>("f_alpha", "free energy of the alpha phase");
  params.addRequiredParam<MaterialPropertyName>("f_gamma", "free_energy_beta");
  params.addRequiredParam<MaterialPropertyName>("df_dcalpha", "diffusion potential in alpha phase");
  params.addRequiredParam<MaterialPropertyName>("df_dcgamma","diffusion potential in beta phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2alpha", "Second derivative of alpha phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2gamma", "Second derivative in beta phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
 
  
  return params;
}

KKS3PhaseEqualGPGammaAlpha::KKS3PhaseEqualGPGammaAlpha(const InputParameters & parameters)
  : Kernel(parameters),
    _phase_comp_alpha(coupledValue("phase_comp_alpha")),
    _phase_comp_alpha_var(coupled("phase_comp_alpha")),
    _phase_comp_gamma(coupledValue("phase_comp_gamma")),
    _phase_comp_gamma_var(coupled("phase_comp_gamma")),
    _phase_beta(coupledValue("phase_beta")),
    _phase_beta_var(coupled("phase_beta")),
    _phase_gamma(coupledValue("phase_gamma")),
    _phase_gamma_var(coupled("phase_gamma")),
    _free_energy_alpha(getMaterialProperty<Real>("f_alpha")),
    _free_energy_gamma(getMaterialProperty<Real>("f_gamma")), 
    _df_alpha(getMaterialProperty<Real>("df_dcalpha")),
    _df_gamma(getMaterialProperty<Real>("df_dcgamma")),
    _d2f_alpha(getMaterialProperty<Real>("d2f_dc2alpha")),
    _d2f_gamma(getMaterialProperty<Real>("d2f_dc2gamma")), 
    // Interpolation functions
    _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),
    _d2hgamma_dphialpha2(getMaterialProperty<Real>("d2hgamma_dphialpha2")),
    _d2hgamma_dphialpha_dphibeta(getMaterialProperty<Real>("d2hgamma_dphialpha_dphibeta")),
    _d2hgamma_dphialpha_dphigamma(getMaterialProperty<Real>("d2hgamma_dphialpha_dphigamma")),
    _L(getMaterialProperty<Real>("mob_name"))
{

}    
        
Real
KKS3PhaseEqualGPGammaAlpha::computeQpResidual()
{
 
    
  const Real omega_gamma = _free_energy_gamma[_qp] - _df_gamma[_qp] * _phase_comp_gamma[_qp];
  const Real omega_alpha = _free_energy_alpha[_qp] - _df_alpha[_qp] * _phase_comp_alpha[_qp];
  
  return _test[_i][_qp] * _L[_qp] * _dhgamma_dphialpha[_qp] * (omega_gamma - omega_alpha);
}   
    
Real
KKS3PhaseEqualGPGammaAlpha::computeQpJacobian()
{
   //Second derivative of the interpolation function: _d2h
   //is required by this kernel
   //The is extracted from the Material class named InterpolationFunction.C
          
  const Real omega_gamma = _free_energy_gamma[_qp] - _df_gamma[_qp] * _phase_comp_gamma[_qp];
  const Real omega_alpha = _free_energy_alpha[_qp] -_df_alpha[_qp] * _phase_comp_alpha[_qp];
            
  return _test[_i][_qp] *_L[_qp] * _d2hgamma_dphialpha2[_qp] * (omega_gamma - omega_alpha) * _phi[_j][_qp];
}

Real
KKS3PhaseEqualGPGammaAlpha::computeQpOffDiagJacobian( unsigned int jvar)
{

  const Real omega_gamma = _free_energy_gamma[_qp] - _df_gamma[_qp] * _phase_comp_gamma[_qp];
  const Real omega_alpha = _free_energy_alpha[_qp] -_df_alpha[_qp] * _phase_comp_alpha[_qp];
    
  if (jvar == _phase_comp_gamma_var)
  {
    // Note the derivatives with respect to omega_{beta} is zero
    return  - _test[_i][_qp] * _L[_qp] * _dhgamma_dphialpha[_qp] * 
       ( _phase_comp_gamma[_qp] * _d2f_gamma[_qp] ) * _phi[_j][_qp]; 
  } 
  else if (jvar == _phase_comp_alpha_var)
  {
    return  _test[_i][_qp] * _L[_qp] * _dhgamma_dphialpha[_qp] * 
    ( _phase_comp_alpha[_qp] * _d2f_alpha[_qp] ) * _phi[_j][_qp];   
  }  
  else if (jvar == _phase_beta_var)
  {
    return  _test[_i][_qp] * _L[_qp] * _d2hgamma_dphialpha_dphibeta[_qp] * 
                                  (omega_gamma - omega_alpha) * _phi[_j][_qp];
  }
  else if (jvar == _phase_gamma_var)
  {
    return _test[_i][_qp] * _L[_qp] * _d2hgamma_dphialpha_dphigamma[_qp] *
                                  (omega_gamma - omega_alpha) * _phi[_j][_qp];                                  
  }
  
  else //anything else
   return 0.0 ;
}
