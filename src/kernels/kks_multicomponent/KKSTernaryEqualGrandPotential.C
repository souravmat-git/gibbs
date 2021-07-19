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
//*This was written by S.Chatterjee

#include "KKSTernaryEqualGrandPotential.h"

registerMooseObject("gibbsApp",KKSTernaryEqualGrandPotential);

template <>
InputParameters
validParams<KKSTernaryEqualGrandPotential>()
{
  InputParameters params = validParams<KKSEqualGrandPotential>();
  params.addClassDescription("Eqn: df/dphi = dh/dphi*(omega_beta - omega_alpha)."
                             "omega_beta, omega_alpha is the grand potential"
                             "This kernel operates on eta.");
  params.addRequiredCoupledVar("xC_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("xC_beta", "Phase concentration in beta phase");
  return params;
}

KKSTernaryEqualGrandPotential::KKSTernaryEqualGrandPotential(const InputParameters & parameters)
  : KKSEqualGrandPotential(parameters),
   _xC_alpha(coupledValue("xC_alpha")),
   _xC_alpha_var(coupled("xC_alpha")),
   _xC_beta(coupledValue("xC_beta")),
   _xC_beta_var(coupled("xC_beta")), 
   _C_diff_pot_alpha(getMaterialProperty<Real>("C_diff_pot_alpha")),
   _C_diff_pot_beta(getMaterialProperty<Real>("C_diff_pot_beta")),
   _C_therm_factor_alpha(getMaterialProperty<Real>("C_therm_factor_alpha")),
   _C_therm_factor_beta(getMaterialProperty<Real>("C_therm_factor_beta")),
   _BC_therm_factor_alpha(getMaterialProperty<Real>("BC_therm_factor_alpha")),
   _BC_therm_factor_beta(getMaterialProperty<Real>("BC_therm_factor_beta"))
    //Note that we do not need to specify the material properties
   // that have already been defined by the class KKSEqualGrandPotential
{
}

Real
KKSTernaryEqualGrandPotential::TernaryGPDiff()
{
  //Defining a member function of type Real 
  //Note omega = f - \tilde{mu}_B*xB - \tilde{mu}_C*xC
  Real _GP_alpha = _f_alpha[_qp] - _B_diff_pot_alpha[_qp] * _xB_alpha[_qp] 
                                 - _C_diff_pot_alpha[_qp] * _xC_alpha[_qp];
  //Now, for phase beta
  Real _GP_beta =  _f_beta[_qp]  - _B_diff_pot_beta[_qp] * _xB_beta[_qp]
                                 - _C_diff_pot_beta[_qp] * _xC_beta[_qp];
                                 
  return (_GP_beta - _GP_alpha);
}
    
    
Real
KKSTernaryEqualGrandPotential::computeQpResidual()
{
  //Note that this kernel uses the residual from the kernel
  //KKSEqualGrandPotential for binary terms
  
  return (_test[_i][_qp] * _L[_qp] * _dh[_qp] * 
                                KKSTernaryEqualGrandPotential::TernaryGPDiff());
}   
    
Real
KKSTernaryEqualGrandPotential::computeQpJacobian()
{

  return (_test[_i][_qp] * _L[_qp] * _d2h[_qp] * 
                                KKSTernaryEqualGrandPotential::TernaryGPDiff());        
}

Real
KKSTernaryEqualGrandPotential::computeQpOffDiagJacobian(unsigned int jvar)
{ 
    
  if (jvar == _xB_alpha_var)
  {
    // Note the derivatives with respect to omega_{beta} is zero
    return   (_test[_i][_qp] * _L[_qp] * _dh[_qp] * ( _xB_alpha[_qp] * _B_therm_factor_alpha[_qp]
             + _BC_therm_factor_alpha[_qp] * _xC_alpha[_qp]) * _phi[_j][_qp]);  
  }
  
  else if (jvar == _xC_alpha_var)
  {
    // Note the derivatives with respect to omega_{beta} is zero
    return   (_test[_i][_qp] * _L[_qp] * _dh[_qp] * ( _xC_alpha[_qp] * _C_therm_factor_alpha[_qp] 
             + _BC_therm_factor_alpha[_qp] * _xB_alpha[_qp])* _phi[_j][_qp]); 
  } 
    
  else if (jvar == _xB_beta_var)
  {
    return - (_test[_i][_qp] * _L[_qp] * _dh[_qp] * ( _xB_beta[_qp] * _B_therm_factor_beta[_qp]
             + _BC_therm_factor_beta[_qp] * _xC_beta[_qp] ) * _phi[_j][_qp]);   
  } 
   
  else if (jvar == _xC_beta_var)
  {
    return - (_test[_i][_qp] * _L[_qp] * _dh[_qp] * ( _xC_beta[_qp] * _C_therm_factor_beta[_qp] 
              + _BC_therm_factor_beta[_qp] * _xB_beta[_qp]) * _phi[_j][_qp]);   
  } 
  else //anything else
   return 0.0 ;
}
