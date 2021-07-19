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
//*In this code both the free energy is obtained 
//* from the material class AnalyticalFreeEnergyMaterial

#include "KKSEqualGrandPotential.h"

registerMooseObject("gibbsApp",KKSEqualGrandPotential);

template <>
InputParameters
validParams<KKSEqualGrandPotential>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: df/dphi = dh/dphi*(omega_beta - omega_alpha)."
                             "omega_beta, omega_alpha is the grand potential"
                             "This kernel operates on eta.");
  params.addRequiredCoupledVar("xB_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("xB_beta", "Phase concentration in beta phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  return params;
}

KKSEqualGrandPotential::KKSEqualGrandPotential(const InputParameters & parameters)
  : Kernel(parameters),
   _xB_alpha(coupledValue("xB_alpha")),
   _xB_alpha_var(coupled("xB_alpha")),
   _xB_beta(coupledValue("xB_beta")),
   _xB_beta_var(coupled("xB_beta")),
   _f_alpha(getMaterialProperty<Real>("f_alpha")),
   _f_beta(getMaterialProperty<Real>("f_beta")), 
   _B_diff_pot_alpha(getMaterialProperty<Real>("B_diff_pot_alpha")),
   _B_diff_pot_beta(getMaterialProperty<Real>("B_diff_pot_beta")),
   _B_therm_factor_alpha(getMaterialProperty<Real>("B_therm_factor_alpha")),
   _B_therm_factor_beta(getMaterialProperty<Real>("B_therm_factor_beta")), 
   _dh(getMaterialProperty<Real>("dh")),
   _d2h(getMaterialProperty<Real>("d2h")),
   _L(getMaterialProperty<Real>("mob_name"))
{
} 
// This function returns the difference in the 
// grandpotential of the two phases
Real 
KKSEqualGrandPotential::grandPotentials_diff()
{
  const Real omega_beta =  _f_beta[_qp]  - _B_diff_pot_beta[_qp]  * _xB_beta[_qp];
  const Real omega_alpha = _f_alpha[_qp] - _B_diff_pot_alpha[_qp] * _xB_alpha[_qp];  

  return (omega_beta - omega_alpha);
}
       
Real
KKSEqualGrandPotential::computeQpResidual()
{ 
  return (_test[_i][_qp] * _L[_qp] * _dh[_qp] *
                              KKSEqualGrandPotential::grandPotentials_diff());
}   
    
Real
KKSEqualGrandPotential::computeQpJacobian()
{
   //Second derivative of the interpolation function: _d2h
   //is required by this kernel
   //The is extracted from the Material class named InterpolationFunction.C            
  return (_test[_i][_qp] *_L[_qp] * _d2h[_qp]* 
            KKSEqualGrandPotential::grandPotentials_diff() * _phi[_j][_qp]);
}

Real
KKSEqualGrandPotential::computeQpOffDiagJacobian(unsigned int jvar)
{   
  if (jvar == _xB_alpha_var)
  {
    return  (_test[_i][_qp] * _L[_qp] * _dh[_qp] * 
                ( _xB_alpha[_qp] * _B_therm_factor_alpha[_qp])* _phi[_j][_qp]); 
  } 
  else if (jvar == _xB_beta_var)
  {
    return -(_test[_i][_qp] * _L[_qp] * _dh[_qp] *
                 ( _xB_beta[_qp] * _B_therm_factor_beta[_qp]) * _phi[_j][_qp]);   
  }  
  else
   return 0.0 ;
}
