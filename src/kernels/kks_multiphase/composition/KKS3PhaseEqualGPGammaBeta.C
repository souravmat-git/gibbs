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

#include "KKS3PhaseEqualGPGammaBeta.h"

registerMooseObject("gibbsApp", KKS3PhaseEqualGPGammaBeta);

template <>
InputParameters
validParams<KKS3PhaseEqualGPGammaBeta>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: df/dphi = dh_{gamma}/dphi{_gamma}*(omega_gamma- omega_beta)."
                             "omega_beta, omega_gamma is the grand potential"
                             "This kernel operates on phi_beta.");
  params.addRequiredCoupledVar("phase_comp_gamma", "Phase concentration in gamma phase");
  params.addRequiredCoupledVar("phase_comp_beta", "Phase concentration in beta phase");
  params.addCoupledVar("phase_alpha",0.0, "Alpha phase field variable");
  params.addRequiredCoupledVar("phase_gamma", "gamma phase field variable");
  params.addRequiredParam<MaterialPropertyName>("f_beta", "free energy of the alpha phase");
  params.addRequiredParam<MaterialPropertyName>("f_gamma", "free_energy_beta");
  params.addRequiredParam<MaterialPropertyName>("df_dcbeta", "diffusion potential in alpha phase");
  params.addRequiredParam<MaterialPropertyName>("df_dcgamma","diffusion potential in beta phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2beta", "Second derivative of alpha phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2gamma", "Second derivative in beta phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
 
  
  return params;
}

KKS3PhaseEqualGPGammaBeta::KKS3PhaseEqualGPGammaBeta(const InputParameters & parameters)
  : Kernel(parameters),
    _phase_comp_gamma(coupledValue("phase_comp_gamma")),
    _phase_comp_gamma_var(coupled("phase_comp_gamma")),
    _phase_comp_beta(coupledValue("phase_comp_beta")),
    _phase_comp_beta_var(coupled("phase_comp_beta")),
    _phase_alpha(coupledValue("phase_alpha")),
    _phase_alpha_var(coupled("phase_alpha")),
    _phase_gamma(coupledValue("phase_gamma")),
    _phase_gamma_var(coupled("phase_gamma")),
    _free_energy_beta(getMaterialProperty<Real>("f_beta")),
    _free_energy_gamma(getMaterialProperty<Real>("f_gamma")), 
    _df_beta(getMaterialProperty<Real>("df_dcbeta")),
    _df_gamma(getMaterialProperty<Real>("df_dcgamma")),
    _d2f_beta(getMaterialProperty<Real>("d2f_dc2beta")),
    _d2f_gamma(getMaterialProperty<Real>("d2f_dc2gamma")), 
    // Interpolation functions
    _dhgamma_dphibeta(getMaterialProperty<Real>("dhbeta_dphigamma")),
    _d2hgamma_dphibeta2(getMaterialProperty<Real>("d2hbeta_dphigamma2")),
    _d2hgamma_dphibeta_dphialpha(getMaterialProperty<Real>("d2hbeta_dphigamma_dphialpha")),
    _d2hgamma_dphibeta_dphigamma(getMaterialProperty<Real>("d2hbeta_dphigamma_dphibeta")),
    _L(getMaterialProperty<Real>("mob_name"))
{
}    
        
Real
KKS3PhaseEqualGPGammaBeta::computeQpResidual()
{
 
    
  const Real omega_beta = _free_energy_beta[_qp] - _df_beta[_qp] * _phase_comp_beta[_qp];
  const Real omega_gamma = _free_energy_gamma[_qp] - _df_gamma[_qp] * _phase_comp_gamma[_qp];
  
  return _test[_i][_qp] * _L[_qp] * _dhgamma_dphibeta[_qp] * (omega_gamma - omega_beta);
}   
    
Real
KKS3PhaseEqualGPGammaBeta::computeQpJacobian()
{
   //Second derivative of the interpolation function: _d2h
   //is required by this kernel
   //The is extracted from the Material class named InterpolationFunction.C
          
  const Real omega_beta = _free_energy_beta[_qp] - _df_beta[_qp] * _phase_comp_beta[_qp];
  const Real omega_gamma = _free_energy_gamma[_qp] -_df_gamma[_qp] * _phase_comp_gamma[_qp];
            
  return _test[_i][_qp] *_L[_qp] * _d2hgamma_dphibeta2[_qp] * (omega_gamma - omega_beta) * _phi[_j][_qp];
}

Real
KKS3PhaseEqualGPGammaBeta::computeQpOffDiagJacobian( unsigned int jvar)
{

  const Real omega_beta = _free_energy_beta[_qp] - _df_beta[_qp] * _phase_comp_beta[_qp];
  const Real omega_gamma = _free_energy_gamma[_qp] -_df_gamma[_qp] * _phase_comp_gamma[_qp];
    
  if (jvar == _phase_comp_beta_var)
  {
    // Note the derivatives with respect to omega_{beta} is zero
    return  _test[_i][_qp] * _L[_qp] * _dhgamma_dphibeta[_qp] * 
       ( _phase_comp_beta[_qp] * _d2f_beta[_qp] ) * _phi[_j][_qp]; 
  } 
  else if (jvar == _phase_comp_gamma_var)
  {
    return - _test[_i][_qp] * _L[_qp] * _dhgamma_dphibeta[_qp] * 
    ( _phase_comp_gamma[_qp] * _d2f_gamma[_qp] ) * _phi[_j][_qp];   
  }  
  else if (jvar == _phase_alpha_var)
  {
    return  _test[_i][_qp] * _L[_qp] * _d2hgamma_dphibeta_dphialpha[_qp] * 
                                  (omega_gamma - omega_beta) * _phi[_j][_qp];
  }
  else if (jvar == _phase_gamma_var)
  {
    return _test[_i][_qp] * _L[_qp] * _d2hgamma_dphibeta_dphigamma[_qp] *
                                  (omega_gamma - omega_beta) * _phi[_j][_qp];                                  
  }
  
  else //anything else
   return 0.0 ;
}
