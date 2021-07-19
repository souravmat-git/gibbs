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

#include "DrivingForceSpinodal.h"

registerMooseObject("gibbsApp", DrivingForceSpinodal);

template <>
InputParameters
validParams<DrivingForceSpinodal>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn:Eqn: dh*(mu_A^{\beta} - mu_A^{\alpha})"
                             "This kernel operates on eta.");
  params.addRequiredCoupledVar("xB", "Mole fraction of component B");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  return params;
}

DrivingForceSpinodal::DrivingForceSpinodal(const InputParameters & parameters)
  : Kernel(parameters),
  _xB(coupledValue("xB")),
  _xB_var(coupled("xB")),
  _f_alpha(getMaterialProperty<Real>("f_alpha")),
  _B_diff_pot_alpha(getMaterialProperty<Real>("B_diff_pot_alpha")),
  _B_therm_factor_alpha(getMaterialProperty<Real>("B_therm_factor_alpha")),
  _A_chem_pot_eqm(getMaterialProperty<Real>("A_chem_pot_eqm")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  _L(getMaterialProperty<Real>("mob_name"))
{
} 
        
Real
DrivingForceSpinodal::computeQpResidual()
{ 
  return (_test[_i][_qp] * _L[_qp] * _dh[_qp] * (_f_alpha[_qp]- _xB[_qp]*_B_diff_pot_alpha[_qp] -_A_chem_pot_eqm[_qp]));
}   
    
Real
DrivingForceSpinodal::computeQpJacobian()
{
   return (_test[_i][_qp] * _L[_qp] * _d2h[_qp] *(_f_alpha[_qp]-_xB[_qp]*_B_diff_pot_alpha[_qp] -_A_chem_pot_eqm[_qp])* _phi[_j][_qp]);
}

Real
DrivingForceSpinodal::computeQpOffDiagJacobian( unsigned int jvar)
{
    
  if (jvar == _xB_var)
  {
    return - (_test[_i][_qp] * _L[_qp] * _dh[_qp] * _xB[_qp] * _B_therm_factor_alpha[_qp] *_phi[_j][_qp]); 
  } 
  else    
    return 0.0; 
}
