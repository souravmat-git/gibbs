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

#include "DiffusionPotentialWBM.h"

registerMooseObject("gibbsApp", DiffusionPotentialWBM);

template<>
InputParameters
validParams<DiffusionPotentialWBM>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: mu - df_beta/dc = 0"
                           "This kernel operates on c (mole_fraction)"),
  params.addRequiredCoupledVar("eta", "Phase field variable"),
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential of comp B is the coupled variable");
  return params;
}

DiffusionPotentialWBM::DiffusionPotentialWBM(const InputParameters & parameters)
  :Kernel(parameters),
  _eta(coupledValue("eta")),
  _eta_var(coupled("eta")),
  _B_diff_pot(coupledValue("B_diff_pot")),
  _B_diff_pot_var(coupled("B_diff_pot")),
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  _B_diff_pot_alpha(getMaterialProperty<Real>("B_diff_pot_alpha")),
  _B_therm_factor_alpha(getMaterialProperty<Real>("B_therm_factor_alpha")),
  _B_diff_pot_beta(getMaterialProperty<Real>("B_diff_pot_beta")),
  _B_therm_factor_beta(getMaterialProperty<Real>("B_therm_factor_beta"))
{
}

Real
DiffusionPotentialWBM::computeQpResidual()
{
  //Since the model is WBM
  //We have two diffusion potentials for each phase
  //and only one composition variable
  //Residual: N_{i} * (mu - df/dc_beta) = 0 ;
  // where N_{i} is the ith weight function
  //The non-linearvariable for this kernel is c1:
  _interpolated_mu = (_h[_qp] * _B_diff_pot_beta[_qp] + (1-_h[_qp]) * _B_diff_pot_alpha[_qp]);
  return _test[_i][_qp] * (_B_diff_pot[_qp] - _interpolated_mu);

}

Real
DiffusionPotentialWBM::computeQpJacobian()
{
  // This will be non-zero since the 
  // free energies of each phase is dependent on composition
  return  -_test[_i][_qp] * (_h[_qp] * _B_therm_factor_beta[_qp] + (1-_h[_qp]) * _B_therm_factor_alpha[_qp]) *_phi[_j][_qp];
}

Real
DiffusionPotentialWBM::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _B_diff_pot_var)
 {
    return _test[_i][_qp] *_phi[_j][_qp];  
 }
 else if (jvar == _eta_var)
 {
    return -_test[_i][_qp] * _dh[_qp] * (_B_diff_pot_beta[_qp] - _B_diff_pot_alpha[_qp])* _phi[_j][_qp];
 }  

 else
 {
    return 0.0;
 }
}
