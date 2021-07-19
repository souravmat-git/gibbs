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

#include "BinaryDiffusionPotentialKKS.h"

registerMooseObject("gibbsApp", BinaryDiffusionPotentialKKS);

template<>
InputParameters
validParams<BinaryDiffusionPotentialKKS>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: mu - df_beta/dc_beta = 0"
                           "This kernel operates on c (mole_fraction)");
  params.addRequiredCoupledVar("xB_beta", "Component B in beta phase");
  params.addRequiredCoupledVar("diff_pot", "diffusion potential of comp B");
  return params;
}

BinaryDiffusionPotentialKKS::BinaryDiffusionPotentialKKS(const InputParameters & parameters)
  : Kernel(parameters),
  _xB_beta(coupledValue("xB_beta")),
  _xB_beta_var(coupled("xB_beta")),
  _diff_pot(coupledValue("diff_pot")),
  _diff_pot_var(coupled("diff_pot")),
  _B_diff_pot_beta(getMaterialProperty<Real>("B_diff_pot_beta")),
  _B_therm_factor_beta(getMaterialProperty<Real>("B_therm_factor_beta"))
{
}

Real
BinaryDiffusionPotentialKKS::computeQpResidual()
{
  
  return (_test[_i][_qp] * (_diff_pot[_qp] - _B_diff_pot_beta[_qp]));
}

Real
BinaryDiffusionPotentialKKS::computeQpJacobian()
{
  return  (0.0); 
}

Real
BinaryDiffusionPotentialKKS::computeQpOffDiagJacobian(unsigned int jvar)
{
 
 if (jvar == _xB_beta_var)
 {
    return -(_test[_i][_qp] * _B_therm_factor_beta[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _diff_pot_var)
 { 
   return (_test[_i][_qp] * _phi[_j][_qp]) ;
 }
 else 
    return 0.0;
 
}
