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

#include "MCContinuityEquationD.h"

registerMooseObject("gibbsApp", MCContinuityEquationD);

template<>
InputParameters
validParams<MCContinuityEquationD>()
{
  InputParameters params = validParams<MultiCompMultiPhaseBase>();
  params.addClassDescription("Continuity equation of component D");
  params.addRequiredCoupledVar("xB", "Mole fraction of comp B");
  params.addRequiredCoupledVar("xC", "Mole fraction of comp C");
  params.addRequiredCoupledVar("D_diff_pot", "Diffusion potential of comp D");
  params.addRequiredCoupledVar("C_diff_pot", "Diffusion potential of comp B");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");
  return params;
}

MCContinuityEquationD::MCContinuityEquationD(const InputParameters & parameters)
  : MultiCompMultiPhaseBase(parameters),
    //Mole fraction of component C
   _grad_xC(coupledGradient("xC")),
   _xC_var(coupled("xC")),
   //Mole fraction of component B
   _grad_xB(coupledGradient("xB")),
   _xB_var(coupled("xB")),
  //Diffusion potential of component D
   _grad_D_diff_pot(coupledGradient("D_diff_pot")),
   _D_diff_pot_var(coupled("D_diff_pot")),
  //Diffusion potential of component C
   _grad_B_diff_pot(coupledGradient("B_diff_pot")),
   _B_diff_pot_var(coupled("B_diff_pot")),
   //Diffusion potential of component B
   _grad_C_diff_pot(coupledGradient("C_diff_pot")),
   _C_diff_pot_var(coupled("C_diff_pot"))
{
}

Real
MCContinuityEquationD::computeQpResidual()
{  
    //acts on Mole fraction of component D
   return (_grad_test[_i][_qp] *(MultiCompMultiPhaseBase::L_DD_interp() * _grad_D_diff_pot[_qp] 
                               + MultiCompMultiPhaseBase::L_CD_interp() * _grad_C_diff_pot[_qp]
                               + MultiCompMultiPhaseBase::L_BD_interp() * _grad_B_diff_pot[_qp]));
}

Real
MCContinuityEquationD::computeQpJacobian()
{ 
    //derivative wrt Mole fraction of component D
   return (_grad_test[_i][_qp] * MultiCompMultiPhaseBase::DC_DD_interp() * _grad_phi[_j][_qp]); 
}

Real
MCContinuityEquationD::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xC_var)
  {
     return (_grad_test[_i][_qp] * MultiCompMultiPhaseBase::DC_DC_interp() * _grad_phi[_j][_qp]); ;
  }
  else if (jvar == _xB_var)
  {
     return (_grad_test[_i][_qp] * MultiCompMultiPhaseBase::DC_DB_interp() * _grad_phi[_j][_qp]); ;
  }
  else if (jvar == _D_diff_pot_var)
  {
      return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_DD_interp() * _grad_phi[_j][_qp])
                                +( (MultiCompMultiPhaseBase::dL_DD_muD_interp() * _grad_D_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_CD_muD_interp() * _grad_C_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_BD_muD_interp() * _grad_B_diff_pot[_qp]) ) * _phi[_j][_qp] ));
  }
  else if (jvar == _C_diff_pot_var)
  {
      return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_CD_interp() * _grad_phi[_j][_qp])
                                 +((MultiCompMultiPhaseBase::dL_DD_muC_interp() * _grad_D_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_CD_muC_interp() * _grad_C_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_BD_muC_interp() * _grad_B_diff_pot[_qp]) ) * _phi[_j][_qp] ));
  }
  else if (jvar == _B_diff_pot_var)
  {
   return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_BD_interp() * _grad_phi[_j][_qp])
                                 +( (MultiCompMultiPhaseBase::dL_DD_muB_interp() * _grad_D_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_CD_muB_interp() * _grad_C_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_BD_muB_interp() * _grad_B_diff_pot[_qp]) ) * _phi[_j][_qp] ));
                               
  }
  else if (jvar == _phase_alpha_var)
  {
    return (_grad_test[_i][_qp] *((MultiCompMultiPhaseBase::sum_dh_L_DD_alpha()*  _grad_D_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_alpha()*  _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_alpha()*  _grad_B_diff_pot[_qp]) ) *_phi[_j][_qp]);
  } 
  else if (jvar == _phase_beta_var)
  {
    return (_grad_test[_i][_qp] *((MultiCompMultiPhaseBase::sum_dh_L_DD_beta()*  _grad_D_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_beta()*  _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_beta()*  _grad_B_diff_pot[_qp]) ) *_phi[_j][_qp]);
  }
  else if (jvar == _phase_gamma_var)
  {
    return (_grad_test[_i][_qp] *((MultiCompMultiPhaseBase::sum_dh_L_DD_gamma()*  _grad_D_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_gamma()*  _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_gamma()*  _grad_B_diff_pot[_qp]) ) *_phi[_j][_qp]);
  } 
  else if (jvar == _phase_delta_var)
  {
    return (_grad_test[_i][_qp] *((MultiCompMultiPhaseBase::sum_dh_L_DD_delta()*  _grad_D_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_delta()*  _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_delta()*  _grad_B_diff_pot[_qp]) ) *_phi[_j][_qp]);
  }
  else if (jvar == _phase_epsilon_var)
  {
    return (_grad_test[_i][_qp] *((MultiCompMultiPhaseBase::sum_dh_L_DD_epsilon()*  _grad_D_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_epsilon()*  _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_epsilon()*  _grad_B_diff_pot[_qp]) ) *_phi[_j][_qp]);
  }
 else //anything else 
    return 0.0;
}
