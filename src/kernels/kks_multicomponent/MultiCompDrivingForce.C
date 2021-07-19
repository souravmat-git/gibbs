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

#include "MultiCompDrivingForce.h"

registerMooseObject("gibbsApp",  MultiCompDrivingForce);

template <>
InputParameters
validParams< MultiCompDrivingForce>()
{
  InputParameters params = validParams<BinaryMultiPhaseDrivingForce>();
  params.addClassDescription("Eqn: mu_{A}^{\beta} - mu_{A}^{\alpha}");
  params.addRequiredCoupledVar("C_diff_pot","Diffusion potential of comp C");
  params.addCoupledVar("D_diff_pot", 0.0, "Diffusion potential of comp D");
  params.addRequiredParam<MaterialPropertyName>("xC_1", "Mole fraction of C in phase 1");
  params.addRequiredParam<MaterialPropertyName>("xC_2", "Mole fraction of C in phase 2");
  params.addParam<MaterialPropertyName>("xD_1", 0.0,"Mole fraction of D in phase 1");
  params.addParam<MaterialPropertyName>("xD_2", 0.0,"Mole fraction of D in phase 2");
  return params;
}

MultiCompDrivingForce::MultiCompDrivingForce(const InputParameters & parameters)
  : BinaryMultiPhaseDrivingForce(parameters),
   _C_diff_pot(coupledValue("C_diff_pot")),
   _C_diff_pot_var(coupled("C_diff_pot")),
   //Check if it is a quaternary alloy ? 
   _D_diff_pot(coupledValue("D_diff_pot")),
   _D_diff_pot_var(coupled("D_diff_pot")),
   _xC_1(getMaterialProperty<Real>("xC_1")),
   _xC_2(getMaterialProperty<Real>("xC_2")),
   _xD_1(getMaterialProperty<Real>("xD_1")),
   _xD_2(getMaterialProperty<Real>("xD_2"))
{
}
    
    
Real
MultiCompDrivingForce::computeQpResidual()
{
  //Residual : L*dh*\mu_A^{\beta} -\mu_A^{\alpha}  

  return (BinaryMultiPhaseDrivingForce::computeQpResidual());
}   
    
Real
MultiCompDrivingForce::computeQpJacobian()
{
 return (BinaryMultiPhaseDrivingForce::computeQpJacobian());        
}

Real
MultiCompDrivingForce::computeQpOffDiagJacobian(unsigned int jvar)
{ 
    
  if (jvar == _B_diff_pot_var)
  {
    return  (BinaryMultiPhaseDrivingForce::computeQpOffDiagJacobian(jvar));  
  } 
  else if (jvar == _phase_2_var)
  {
    return (BinaryMultiPhaseDrivingForce::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_3_var)
  {
    return (BinaryMultiPhaseDrivingForce::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_4_var)
  {
    return (BinaryMultiPhaseDrivingForce::computeQpOffDiagJacobian(jvar));
  }
  else if (jvar == _phase_5_var)
  {
    return (BinaryMultiPhaseDrivingForce::computeQpOffDiagJacobian(jvar));
  } 
  else if (jvar == _C_diff_pot_var)
  {
    return  (_test[_i][_qp] * _L[_qp] *_nd_factor[_qp]* _dh[_qp] * (_xC_1[_qp] - _xC_2[_qp]) * _phi[_j][_qp]); 
  }  
  else if (jvar == _D_diff_pot_var)
  {
    return  (_test[_i][_qp] * _L[_qp] * _nd_factor[_qp]*_dh[_qp] * (_xD_1[_qp] - _xD_2[_qp]) * _phi[_j][_qp]);
  }
  else //anything else
   return 0.0 ;
}
