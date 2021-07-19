//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html 

//* This kernel implements the chemical potential
//* The equation this kernel mu - df/dc = 0
//* mu is the variable that this kernel operates on

#include "DiffusionPotentialKKSC.h"

registerMooseObject("gibbsApp", DiffusionPotentialKKSC);

template<>
InputParameters
validParams<DiffusionPotentialKKSC>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: mu_C - df_beta/dX_beta = 0"
                           "This kernel operates on X_C (mole_fraction)");
  params.addRequiredCoupledVar("xC_beta", "Component C in beta phase");
  params.addRequiredCoupledVar("xB_beta", "Component B in beta phase");
  params.addCoupledVar("xD_beta",0.0, "Component D in beta phase");
  params.addCoupledVar("xC_alpha",0.0, "phase composition in alpha phase"); 
  params.addCoupledVar("xC_gamma", 0.0, "phase_composition in gamma phase");
  params.addRequiredCoupledVar("diff_pot", "diffusion potential is the coupled variable");
  params.addParam<MaterialPropertyName>("CD_therm_factor_beta", 0.0 ,"mu(C).x(D)");
  return params;
}

DiffusionPotentialKKSC::DiffusionPotentialKKSC(const InputParameters & parameters)
  : Kernel(parameters),
  _xC_beta(coupledValue("xC_beta")),
  _xC_beta_var(coupled("xC_beta")),
  _xB_beta(coupledValue("xB_beta")),
  _xB_beta_var(coupled("xB_beta")),
  _xD_beta(coupledValue("xD_beta")),
  _xD_beta_var(coupled("xD_beta")),
  _xC_alpha(coupledValue("xC_alpha")),
  _xC_alpha_var(coupled("xC_alpha")),
  _xC_gamma(coupledValue("xC_gamma")),
  _xC_gamma_var(coupled("xC_gamma")),
  _diff_pot(coupledValue("diff_pot")),
  _diff_pot_var(coupled("diff_pot")),
  _C_diff_pot_beta(getMaterialProperty<Real>("C_diff_pot_beta")),
  _C_therm_factor_beta(getMaterialProperty<Real>("C_therm_factor_beta")),
  _BC_therm_factor_beta(getMaterialProperty<Real>("BC_therm_factor_beta")),
  _CD_therm_factor_beta(getMaterialProperty<Real>("CD_therm_factor_beta"))
{
}

Real
DiffusionPotentialKKSC::computeQpResidual()
{
  //This kernel takes df/dbeta Note: df/dc = df/dalpha = df/dbeta
  // diff_pot is with resepct to mole fraction
  //and diff_pot_beta is with respect to phase comp beta
  //Residual: N_{i} * (mu - df/dc_beta) = 0 ;
  // where N_{i} is the ith weight function
  //The non-linearvariable for this kernel is c1:
  
  return _test[_i][_qp] * (_diff_pot[_qp] - _C_diff_pot_beta[_qp]);
}

Real
DiffusionPotentialKKSC::computeQpJacobian()
{
  return  0.0; 
}

Real
DiffusionPotentialKKSC::computeQpOffDiagJacobian(unsigned int jvar)
{
 
 if (jvar == _xC_beta_var)
 {
    //Note: d2f/dc_beta = A_beta (only if parabola)
    return -(_test[_i][_qp] * _C_therm_factor_beta[_qp] * _phi[_j][_qp]);
 }
 else if  (jvar == _xB_beta_var)
 {
    return -(_test[_i][_qp] * _BC_therm_factor_beta[_qp]* _phi[_j][_qp]);
 }
 else if (jvar == _xD_beta_var)
 {
    return -(_test[_i][_qp] * _CD_therm_factor_beta[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _diff_pot_var)
 { 
  return  _test[_i][_qp] * _phi[_j][_qp] ;
 }
 else if (jvar == _xC_alpha_var)
 { 
    return 0.0;
    //Note : This will be zero for any solution model
  }
  
 else if (jvar == _xC_gamma_var) 
 {
    return 0.0;
 }
 else //anything else
 {
    return 0.0;
 }
 
}
