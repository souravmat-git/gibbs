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

#include "DrivingForceKKS.h"

registerMooseObject("gibbsApp", DrivingForceKKS);

template <>
InputParameters
validParams<DrivingForceKKS>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn:Eqn: dh*(mu_A^{\beta} - mu_A^{\alpha})"
                             "This kernel operates on eta.");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of component B");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "Non-dim. phase field mobility");
  params.addParam<MaterialPropertyName>("nd_factor", 1.0, "RT/Vm*barrier_height");
  return params;
}

DrivingForceKKS::DrivingForceKKS(const InputParameters & parameters)
  : Kernel(parameters),
   _B_diff_pot(coupledValue("B_diff_pot")),
   _B_diff_pot_var(coupled("B_diff_pot")),
   _A_chem_pot_alpha(getMaterialProperty<Real>("A_chem_pot_alpha")),
   _A_chem_pot_beta(getMaterialProperty<Real>("A_chem_pot_beta")), 
   _xB_alpha(getMaterialProperty<Real>("xB_alpha")),
   _xB_beta(getMaterialProperty<Real>("xB_beta")),
   _dh(getMaterialProperty<Real>("dh")),
   _d2h(getMaterialProperty<Real>("d2h")),
   _L(getMaterialProperty<Real>("mob_name")),
   _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
} 
        
Real
DrivingForceKKS::computeQpResidual()
{
  
  return (_nd_factor[_qp]*(_test[_i][_qp] * _L[_qp] * _dh[_qp] * (_A_chem_pot_beta[_qp] 
                                                  -_A_chem_pot_alpha[_qp])));
}   
    
Real
DrivingForceKKS::computeQpJacobian()
{
   return (_nd_factor[_qp]*(_test[_i][_qp] * _L[_qp] * _d2h[_qp] * (_A_chem_pot_beta[_qp] 
                                                  - _A_chem_pot_alpha[_qp])* _phi[_j][_qp]));
}

Real
DrivingForceKKS::computeQpOffDiagJacobian( unsigned int jvar)
{
    
  if (jvar == _B_diff_pot_var)
  {
    return (_nd_factor[_qp]*(_test[_i][_qp] * _L[_qp] * _dh[_qp] * (_xB_alpha[_qp] - _xB_beta[_qp])* _phi[_j][_qp])); 
  } 
  else    
      return 0.0;
 
}
