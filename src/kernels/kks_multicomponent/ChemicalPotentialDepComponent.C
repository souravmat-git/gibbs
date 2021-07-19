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

#include "ChemicalPotentialDepComponent.h"

registerMooseObject("gibbsApp", ChemicalPotentialDepComponent);

template <>
InputParameters
validParams<ChemicalPotentialDepComponent>()
{
  InputParameters params = validParams<ChemicalPotentialDepComponentBase>();
  params.addClassDescription("Eqn: mu_{A}^{\beta} - mu_{A}^{\alpha}");
  params.addRequiredCoupledVar("xC_alpha","Component C in alpha phase");
  params.addRequiredCoupledVar("xC_beta","Component C in beta phase");
  params.addCoupledVar("xD_alpha", 0.0, "Component D in alpha phase");
  params.addCoupledVar("xD_beta", 0.0, "Component D in beta phase");
  params.addParam<MaterialPropertyName>("AD_therm_factor_alpha", 0.0,"mu(A).x(D)");
  params.addParam<MaterialPropertyName>("AD_therm_factor_beta",0.0,"mu(A).x(D)_beta");
  return params;
}

ChemicalPotentialDepComponent::ChemicalPotentialDepComponent(const InputParameters & parameters)
  : ChemicalPotentialDepComponentBase(parameters),
   _xC_alpha(coupledValue("xC_alpha")),
   _xC_alpha_var(coupled("xC_alpha")),
   _xC_beta(coupledValue("xC_beta")),
   _xC_beta_var(coupled("xC_beta")),
   //Check if it is a quaternary alloy ? 
   _xD_alpha(coupledValue("xD_alpha")),
   _xD_alpha_var(coupled("xD_alpha")),
   _xD_beta(coupledValue("xD_beta")),
   _xD_beta_var(coupled("xD_beta")),
   _AC_therm_factor_alpha(getMaterialProperty<Real>("AC_therm_factor_alpha")),
   _AC_therm_factor_beta(getMaterialProperty<Real>("AC_therm_factor_beta")),
   _AD_therm_factor_alpha(getMaterialProperty<Real>("AD_therm_factor_alpha")),
   _AD_therm_factor_beta(getMaterialProperty<Real>("AD_therm_factor_beta"))
{
}
    
    
Real
ChemicalPotentialDepComponent::computeQpResidual()
{
  //Residual : L*dh*\mu_A^{\beta} -\mu_A^{\alpha}
  
  return (ChemicalPotentialDepComponentBase::computeQpResidual());
}   
    
Real
ChemicalPotentialDepComponent::computeQpJacobian()
{
 return (ChemicalPotentialDepComponentBase::computeQpJacobian());        
}

Real
ChemicalPotentialDepComponent::computeQpOffDiagJacobian(unsigned int jvar)
{    
  if (jvar == _xB_alpha_var)
  {
    return  ChemicalPotentialDepComponentBase::computeQpOffDiagJacobian(jvar);  
  }  
  else if (jvar == _xC_alpha_var)
  {
    return -(_test[_i][_qp] * _L[_qp] * _dh[_qp] * _AC_therm_factor_alpha[_qp]* _phi[_j][_qp]); 
  } 
  else if (jvar == _xD_alpha_var)
  {
    return -(_test[_i][_qp] * _L[_qp]* _dh[_qp] * _AD_therm_factor_alpha[_qp] * _phi[_j][_qp]);
  } 
  else if (jvar == _xB_beta_var)
  {
    return ChemicalPotentialDepComponentBase::computeQpOffDiagJacobian(jvar);   
  }     
  else if (jvar == _xC_beta_var)
  {
    return (_test[_i][_qp] * _L[_qp] * _dh[_qp] * _AC_therm_factor_beta[_qp]* _phi[_j][_qp]);   
  }   
  else if (jvar == _xD_beta_var)
  {
    return (_test[_i][_qp] * _L[_qp] * _dh[_qp] * _AD_therm_factor_beta[_qp] * _phi[_j][_qp]);
  }
  else //anything else
   return 0.0 ;
}
