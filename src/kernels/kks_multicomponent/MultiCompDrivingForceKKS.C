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

#include "MultiCompDrivingForceKKS.h"

registerMooseObject("gibbsApp",  MultiCompDrivingForceKKS);

template <>
InputParameters
validParams< MultiCompDrivingForceKKS>()
{
  InputParameters params = validParams<DrivingForceKKS>();
  params.addClassDescription("Eqn: mu_{A}^{\beta} - mu_{A}^{\alpha}");
  params.addRequiredCoupledVar("C_diff_pot","Diffusion potential of comp C");
  params.addCoupledVar("D_diff_pot", 0.0, "Diffusion potential of comp C");
  params.addParam<MaterialPropertyName>("xD_alpha", 0.0,"mu(A).mu(D)");
  params.addParam<MaterialPropertyName>("xD_beta", 0.0,"mu(A).mu(D)");
  return params;
}

MultiCompDrivingForceKKS::MultiCompDrivingForceKKS(const InputParameters & parameters)
  : DrivingForceKKS(parameters),
   _C_diff_pot(coupledValue("C_diff_pot")),
   _C_diff_pot_var(coupled("C_diff_pot")),
   //Check if it is a quaternary alloy ? 
   _D_diff_pot(coupledValue("D_diff_pot")),
   _D_diff_pot_var(coupled("D_diff_pot")),
   _xC_alpha(getMaterialProperty<Real>("xC_alpha")),
   _xC_beta(getMaterialProperty<Real>("xC_beta")),
   _xD_alpha(getMaterialProperty<Real>("xD_alpha")),
   _xD_beta(getMaterialProperty<Real>("xD_beta"))
{
}
    
    
Real
MultiCompDrivingForceKKS::computeQpResidual()
{
  //Residual : L*dh*\mu_A^{\beta} -\mu_A^{\alpha}  
  return (DrivingForceKKS::computeQpResidual());
}   
    
Real
MultiCompDrivingForceKKS::computeQpJacobian()
{
 return (DrivingForceKKS::computeQpJacobian());        
}

Real
MultiCompDrivingForceKKS::computeQpOffDiagJacobian(unsigned int jvar)
{ 
    
  if (jvar == _B_diff_pot_var)
  {
    return  (DrivingForceKKS::computeQpOffDiagJacobian(jvar));  
  }  
  else if (jvar == _C_diff_pot_var)
  {
    return  (_test[_i][_qp] * _L[_qp] * _dh[_qp] * (- _xC_beta[_qp] + _xC_alpha[_qp]) * _phi[_j][_qp]); 
  }  
  else if (jvar == _D_diff_pot_var)
  {
    return  (_test[_i][_qp] * _L[_qp] * _dh[_qp] * (- _xD_beta[_qp] + _xD_alpha[_qp]) * _phi[_j][_qp]);
  } 
  else //anything else
   return 0.0 ;
}
