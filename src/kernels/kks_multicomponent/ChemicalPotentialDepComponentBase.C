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

#include "ChemicalPotentialDepComponentBase.h"

registerMooseObject("gibbsApp",ChemicalPotentialDepComponentBase);

template <>
InputParameters
validParams<ChemicalPotentialDepComponentBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*(mu_A^{\beta} - mu_A^{\alpha})"
                             "This kernel operates on eta.");
  params.addRequiredCoupledVar("xB_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("xB_beta", "Phase concentration in beta phase");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  return params;
}

ChemicalPotentialDepComponentBase::ChemicalPotentialDepComponentBase(const InputParameters & parameters)
  : Kernel(parameters),
   _xB_alpha(coupledValue("xB_alpha")),
   _xB_alpha_var(coupled("xB_alpha")),
   _xB_beta(coupledValue("xB_beta")),
   _xB_beta_var(coupled("xB_beta")),
   _A_chem_pot_alpha(getMaterialProperty<Real>("A_chem_pot_alpha")),
   _A_chem_pot_beta(getMaterialProperty<Real>("A_chem_pot_beta")),
   _AB_therm_factor_alpha(getMaterialProperty<Real>("AB_therm_factor_alpha")),
   _AB_therm_factor_beta(getMaterialProperty<Real>("AB_therm_factor_beta")),
   _dh(getMaterialProperty<Real>("dh")),
   _d2h(getMaterialProperty<Real>("d2h")),
   _L(getMaterialProperty<Real>("mob_name"))
{
} 

Real
ChemicalPotentialDepComponentBase::computeQpResidual()
{
  //Residual : L*dh*\mu_A^{\beta} -\mu_A^{\alpha}
  
  return (_test[_i][_qp] * _L[_qp] * _dh[_qp] * (_A_chem_pot_beta[_qp] 
                                                -_A_chem_pot_alpha[_qp]));
}   
    
Real
ChemicalPotentialDepComponentBase::computeQpJacobian()
{
 return (_test[_i][_qp] * _L[_qp] * _d2h[_qp] * (_A_chem_pot_beta[_qp] 
                                                -_A_chem_pot_alpha[_qp]) * _phi[_j][_qp]);        
}

Real
ChemicalPotentialDepComponentBase::computeQpOffDiagJacobian(unsigned int jvar)
{ 
    
  if (jvar == _xB_alpha_var)
  {

    return  -(_test[_i][_qp] * _L[_qp] * _dh[_qp] *_AB_therm_factor_alpha[_qp] * _phi[_j][_qp]);  
  }
    
  else if (jvar == _xB_beta_var)
  {
    return (_test[_i][_qp] * _L[_qp] * _dh[_qp] * _AB_therm_factor_beta[_qp]* _phi[_j][_qp]);   
  } 
   
  else //anything else
   return 0.0 ;
}
