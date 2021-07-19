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
//*Allen Cahn equation which tells us that at
//*equilibrium the grandpotentials must be equal
//*In this code both the free energy can be 
//*extracted from the input file

#include "KKS3PhaseEqualGPAlphaBeta.h"

registerMooseObject("gibbsApp",KKS3PhaseEqualGPAlphaBeta);

template <>
InputParameters
validParams<KKS3PhaseEqualGPAlphaBeta>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: df/dphi_{beta} = dh_{alpha}/dphi{_beta}*(omega_alpha - omega_beta)."
                             "omega_beta, omega_alpha is the grand potential"
                             "This kernel operates on phi_beta.");
  params.addRequiredCoupledVar("phase_comp_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("phase_comp_beta", "Phase concentration in beta phase");
  params.addRequiredCoupledVar("phase_alpha", "Alpha phase field variable");
  params.addCoupledVar("phase_gamma", 0.0,"Gamma phase field variable");
  params.addRequiredParam<MaterialPropertyName>("f_alpha", "free energy of the alpha phase");
  params.addRequiredParam<MaterialPropertyName>("f_beta", "free_energy_beta");
  params.addRequiredParam<MaterialPropertyName>("df_dcalpha", "diffusion potential in alpha phase");
  params.addRequiredParam<MaterialPropertyName>("df_dcbeta","diffusion potential in beta phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2alpha", "Second derivative of alpha phase");
  params.addRequiredParam<MaterialPropertyName>("d2f_dc2beta", "Second derivative in beta phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
 
  
  return params;
}

KKS3PhaseEqualGPAlphaBeta::KKS3PhaseEqualGPAlphaBeta(const InputParameters & parameters)
  : Kernel(parameters),
    _phase_comp_alpha(coupledValue("phase_comp_alpha")),
    _phase_comp_alpha_var(coupled("phase_comp_alpha")),
    _phase_comp_beta(coupledValue("phase_comp_beta")),
    _phase_comp_beta_var(coupled("phase_comp_beta")),
    _phase_alpha(coupledValue("phase_alpha")),
    _phase_alpha_var(coupled("phase_alpha")),
    _phase_gamma(coupledValue("phase_gamma")),
    _phase_gamma_var(coupled("phase_gamma")),
    _free_energy_alpha(getMaterialProperty<Real>("f_alpha")),
    _free_energy_beta(getMaterialProperty<Real>("f_beta")), 
    _df_alpha(getMaterialProperty<Real>("df_dcalpha")),
    _df_beta(getMaterialProperty<Real>("df_dcbeta")),
    _d2f_alpha(getMaterialProperty<Real>("d2f_dc2alpha")),
    _d2f_beta(getMaterialProperty<Real>("d2f_dc2beta")), 
    // Interpolation functions
    _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
    _d2halpha_dphibeta2(getMaterialProperty<Real>("d2halpha_dphibeta2")),
    _d2halpha_dphibeta_dphialpha(getMaterialProperty<Real>("d2halpha_dphibeta_dphialpha")),
    _d2halpha_dphibeta_dphigamma(getMaterialProperty<Real>("d2halpha_dphibeta_dphigamma")),
    _L(getMaterialProperty<Real>("mob_name"))
{
}    
        
Real
KKS3PhaseEqualGPAlphaBeta::computeQpResidual()
{
 
  const Real omega_alpha = _free_energy_alpha[_qp] - _df_alpha[_qp] * _phase_comp_alpha[_qp];
  const Real omega_beta = _free_energy_beta[_qp] - _df_beta[_qp] * _phase_comp_beta[_qp];

  
  return _test[_i][_qp] * _L[_qp] * _dhalpha_dphibeta[_qp] * (omega_alpha - omega_beta);
}   
    
Real
KKS3PhaseEqualGPAlphaBeta::computeQpJacobian()
{
   //Second derivative of the interpolation function: _d2h
   //is required by this kernel
   //The is extracted from the Material class named InterpolationFunction.C
          
  const Real omega_alpha = _free_energy_alpha[_qp] - _df_alpha[_qp] * _phase_comp_alpha[_qp];
  const Real omega_beta = _free_energy_beta[_qp] - _df_beta[_qp] * _phase_comp_beta[_qp];
            
  return _test[_i][_qp] *_L[_qp] * _d2halpha_dphibeta2[_qp] * (omega_alpha - omega_beta) * _phi[_j][_qp];
}

Real
KKS3PhaseEqualGPAlphaBeta::computeQpOffDiagJacobian( unsigned int jvar)
{

  const Real omega_alpha = _free_energy_alpha[_qp] - _df_alpha[_qp] * _phase_comp_alpha[_qp];
  const Real omega_beta = _free_energy_beta[_qp] - _df_beta[_qp] * _phase_comp_beta[_qp];
    
  if (jvar == _phase_comp_beta_var)
  {
    // Note the negative sign here.
    return  _test[_i][_qp] * _L[_qp] * _dhalpha_dphibeta[_qp] * 
       ( _phase_comp_beta[_qp] * _d2f_beta[_qp] ) * _phi[_j][_qp]; 
  } 
  else if (jvar == _phase_comp_alpha_var)
  {
    return - _test[_i][_qp] * _L[_qp] * _dhalpha_dphibeta[_qp] * 
    ( _phase_comp_alpha[_qp] * _d2f_alpha[_qp] ) * _phi[_j][_qp];   
  }  
  else if (jvar == _phase_alpha_var)
  {
    return  _test[_i][_qp] * _L[_qp] * _d2halpha_dphibeta_dphialpha[_qp] * 
                                  (omega_alpha - omega_beta) * _phi[_j][_qp];
  }
  else if (jvar == _phase_gamma_var)
  {
    return _test[_i][_qp] * _L[_qp] * _d2halpha_dphibeta_dphigamma[_qp] *
                                  (omega_alpha - omega_beta) * _phi[_j][_qp];                                  
  }
  
  else //anything else
   return 0.0 ;
}
