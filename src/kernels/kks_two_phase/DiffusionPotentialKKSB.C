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

#include "DiffusionPotentialKKSB.h"

registerMooseObject("gibbsApp", DiffusionPotentialKKSB);

template<>
InputParameters
validParams<DiffusionPotentialKKSB>()
{
  InputParameters params = validParams<BinaryDiffusionPotentialKKS>();
  params.addClassDescription("Eqn: mu - df_beta/dc_beta = 0"
                           "This kernel operates on c (mole_fraction)");
  params.addRequiredCoupledVar("xC_beta", "Component C in beta phase");
  params.addCoupledVar("xD_beta",0.0, "Component D in beta phase");
  params.addRequiredParam<MaterialPropertyName>("BC_therm_factor_beta", "mu(B).x(C) - mu(B).x(C)");
  params.addParam<MaterialPropertyName>("BD_therm_factor_beta", 0.0,"mu(B).x(D) - mu(A).x(D)");
  return params;
}

DiffusionPotentialKKSB::DiffusionPotentialKKSB(const InputParameters & parameters)
  :BinaryDiffusionPotentialKKS(parameters),
  //For a ternary alloy A-B-C
  _xC_beta(coupledValue("xC_beta")),
  _xC_beta_var(coupled("xC_beta")),
  //For a quaternary alloy A-B-C-D
  _xD_beta(coupledValue("xD_beta")),
  _xD_beta_var(coupled("xD_beta")),
  _BC_therm_factor_beta(getMaterialProperty<Real>("BC_therm_factor_beta")),
  _BD_therm_factor_beta(getMaterialProperty<Real>("BD_therm_factor_beta"))
{
}

Real
DiffusionPotentialKKSB::computeQpResidual()
{
  
  return BinaryDiffusionPotentialKKS::computeQpResidual();
}

Real
DiffusionPotentialKKSB::computeQpJacobian()
{
  return  BinaryDiffusionPotentialKKS::computeQpJacobian();
}

Real
DiffusionPotentialKKSB::computeQpOffDiagJacobian(unsigned int jvar)
{
 
 if (jvar == _xB_beta_var)
 {
    return BinaryDiffusionPotentialKKS::computeQpOffDiagJacobian(jvar);
 }
 else if (jvar == _diff_pot_var)
 { 
    return  BinaryDiffusionPotentialKKS::computeQpOffDiagJacobian(jvar);
 }
 else if (jvar == _xC_beta_var)
 {
    return -(_test[_i][_qp] * _BC_therm_factor_beta[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _xD_beta_var)
 {
    return (_test[_i][_qp] * _BD_therm_factor_beta[_qp] * _phi[_j][_qp]);
 }
 else //anything else
 {
    return 0.0;
 }
 
}
