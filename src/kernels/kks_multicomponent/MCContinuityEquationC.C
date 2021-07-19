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

#include "MCContinuityEquationC.h"

registerMooseObject("gibbsApp", MCContinuityEquationC);

template<>
InputParameters
validParams<MCContinuityEquationC>()
{
  InputParameters params = validParams<MultiCompMultiPhaseBase>();
  params.addClassDescription("Continuity equation for comp C");
  params.addRequiredCoupledVar("xB", "Mole fraction of comp B");
  params.addRequiredCoupledVar("xD", "Mole fraction of comp D");
  params.addRequiredCoupledVar("C_diff_pot", "Diffusion potential of comp C");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");
  params.addRequiredCoupledVar("D_diff_pot", "Diffusion potential of comp D");
  return params;
}

MCContinuityEquationC::MCContinuityEquationC(const InputParameters & parameters)
  : MultiCompMultiPhaseBase(parameters),
   //Mole fraction of componnet B
   _grad_xB(coupledGradient("xB")),
   _xB_var(coupled("xB")),
   //Mole fraction of component D
   _grad_xD(coupledGradient("xD")),
   _xD_var(coupled("xD")),
   //Diffusion potential of componnet C
   _grad_C_diff_pot(coupledGradient("C_diff_pot")),
   _C_diff_pot_var(coupled("C_diff_pot")),
   //Diffusion potential of componnet B
   _grad_B_diff_pot(coupledGradient("B_diff_pot")),
   _B_diff_pot_var(coupled("B_diff_pot")),
   //Diffusion potential of component D
   _grad_D_diff_pot(coupledGradient("D_diff_pot")),
   _D_diff_pot_var(coupled("D_diff_pot"))
{
}

Real
MCContinuityEquationC::computeQpResidual()
{ 
   return (_grad_test[_i][_qp] *(MultiCompMultiPhaseBase::L_CC_interp() * _grad_C_diff_pot[_qp] 
                               + MultiCompMultiPhaseBase::L_BC_interp() * _grad_B_diff_pot[_qp]
                               + MultiCompMultiPhaseBase::L_CD_interp() * _grad_D_diff_pot[_qp]));
}

Real
MCContinuityEquationC::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * (MultiCompMultiPhaseBase::DC_CC_interp() * _grad_phi[_j][_qp]) );
                                //+(MultiCompMultiPhaseBase::dD_CC_xC_interp()* _grad_u[_qp]* _phi[_j][_qp])) ); 
}

Real
MCContinuityEquationC::computeQpOffDiagJacobian(unsigned int jvar)
{
 if (jvar == _xB_var)
 {
    return (_grad_test[_i][_qp] * MultiCompMultiPhaseBase::DC_CB_interp() * _grad_phi[_j][_qp]);
 }
 else if (jvar == _xD_var)
 {
    return (_grad_test[_i][_qp] * MultiCompMultiPhaseBase::DC_CD_interp() * _grad_phi[_j][_qp]);
 }
 else if (jvar == _C_diff_pot_var)
 { 
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_CC_interp()* _grad_phi[_j][_qp])
                                 +((MultiCompMultiPhaseBase::dL_CC_muC_interp() * _grad_C_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_BC_muC_interp() * _grad_B_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_CD_muC_interp() * _grad_D_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }
 else if (jvar == _B_diff_pot_var)
 {
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_BC_interp()* _grad_phi[_j][_qp])
                                 +( (MultiCompMultiPhaseBase::dL_CC_muB_interp() * _grad_C_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_BC_muB_interp() * _grad_B_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_CD_muB_interp() * _grad_D_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }
 else if (jvar == _D_diff_pot_var)
 {   
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_CD_interp()* _grad_phi[_j][_qp])
                                 +( (MultiCompMultiPhaseBase::dL_CC_muD_interp() * _grad_C_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_BC_muD_interp() * _grad_B_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_CD_muD_interp() * _grad_D_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }
 else if (jvar == _phase_alpha_var)
 {
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_CC_alpha()* _grad_C_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_alpha()* _grad_B_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_alpha()* _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_beta_var)
 {
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_CC_beta()* _grad_C_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_beta()*  _grad_B_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_beta()*  _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_gamma_var)
 {
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_CC_gamma()* _grad_C_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_gamma()*  _grad_B_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_gamma()*  _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_delta_var)
 {
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_CC_delta()* _grad_C_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_delta()*  _grad_B_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_delta()*  _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_epsilon_var)
 {
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_CC_epsilon()* _grad_C_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_epsilon()*  _grad_B_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_CD_epsilon()*  _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
                                
 }  
 else //anything else
 
    return 0.0;
}
