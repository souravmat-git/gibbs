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

#include "MCContinuityEquationB.h"

registerMooseObject("gibbsApp", MCContinuityEquationB);

template<>
InputParameters
validParams<MCContinuityEquationB>()
{
  InputParameters params = validParams<MultiCompMultiPhaseBase>();
  params.addClassDescription("Continuity equation for component B");
  params.addRequiredCoupledVar("xC", "Mole fraction of component C");
  params.addRequiredCoupledVar("xD", "Mole fraction of component D");  
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");    
  params.addRequiredCoupledVar("C_diff_pot", "diffusion potential of comp C");
  params.addRequiredCoupledVar("D_diff_pot", "diffusion potential of comp D");
  return params;
}

MCContinuityEquationB::MCContinuityEquationB(const InputParameters & parameters)
  : MultiCompMultiPhaseBase(parameters),
   //Mole fraction of component C
   _grad_xC(coupledGradient("xC")),
   _xC_var(coupled("xC")),
  //Mole fraction of component D
   _grad_xD(coupledGradient("xD")),
   _xD_var(coupled("xD")),  
   //Diffusion potential of component B
   _grad_B_diff_pot(coupledGradient("B_diff_pot")),
   _B_diff_pot_var(coupled("B_diff_pot")),
   //Diffusion potential of component C
   _grad_C_diff_pot(coupledGradient("C_diff_pot")),
   _C_diff_pot_var(coupled("C_diff_pot")),
   //Diffusion potential of component D
   _grad_D_diff_pot(coupledGradient("D_diff_pot")),
   _D_diff_pot_var(coupled("D_diff_pot"))
{
}


Real
MCContinuityEquationB::computeQpResidual()
{
   //the variable that this kernel acts on is xB
   return (_grad_test[_i][_qp] *(MultiCompMultiPhaseBase::L_BB_interp() * _grad_B_diff_pot[_qp] 
                               + MultiCompMultiPhaseBase::L_BC_interp() * _grad_C_diff_pot[_qp]
                               + MultiCompMultiPhaseBase::L_BD_interp() * _grad_D_diff_pot[_qp]) );
}

Real
MCContinuityEquationB::computeQpJacobian()
{
  // The jacobian is J_B = D_BB*XB + D_BC*xC + D_BD*XD
  return (_grad_test[_i][_qp] *(MultiCompMultiPhaseBase::DC_BB_interp() * _grad_phi[_j][_qp]) );
}

Real
MCContinuityEquationB::computeQpOffDiagJacobian(unsigned int jvar)
{
 if (jvar == _xC_var)
 {
    return (_grad_test[_i][_qp] * (MultiCompMultiPhaseBase::DC_BC_interp() * _grad_phi[_j][_qp]) );
 }
 else if (jvar == _xD_var)
 {    
     return (_grad_test[_i][_qp] *(MultiCompMultiPhaseBase::DC_BD_interp() * _grad_phi[_j][_qp]) );
 }
 else if (jvar == _B_diff_pot_var) 
 { 
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_BB_interp() * _grad_phi[_j][_qp])
                                 +((MultiCompMultiPhaseBase::dL_BB_muB_interp() * _grad_B_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_BC_muB_interp() * _grad_C_diff_pot[_qp])
                                 + (MultiCompMultiPhaseBase::dL_BD_muB_interp() * _grad_D_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }
 else if (jvar == _C_diff_pot_var) 
 {
 
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_BC_interp() * _grad_phi[_j][_qp])
                                 +( (MultiCompMultiPhaseBase::dL_BB_muC_interp() * _grad_B_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_BC_muC_interp() * _grad_C_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_BD_muC_interp() * _grad_D_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }
 else if (jvar == _D_diff_pot_var)
 {
    
    return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::L_BD_interp() * _grad_phi[_j][_qp])
                                 +( (MultiCompMultiPhaseBase::dL_BB_muD_interp() * _grad_B_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_BC_muD_interp() * _grad_C_diff_pot[_qp])
                                 +  (MultiCompMultiPhaseBase::dL_BD_muD_interp() * _grad_D_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }
 else if (jvar == _phase_alpha_var)
 {
   return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_BB_alpha()* _grad_B_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_alpha()* _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_alpha()* _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
 }
 else if (jvar == _phase_beta_var)
 {
   return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_BB_beta()* _grad_B_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_beta()* _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_beta()* _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
 } 
  else if (jvar == _phase_gamma_var)
 {
   return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_BB_gamma()* _grad_B_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_gamma()* _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_gamma()* _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
 }
  else if (jvar == _phase_delta_var)
 {
   return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_BB_delta()* _grad_B_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_delta()* _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_delta()* _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
 }
  else if (jvar == _phase_epsilon_var)
 {
   return (_grad_test[_i][_qp] *( (MultiCompMultiPhaseBase::sum_dh_L_BB_epsilon()* _grad_B_diff_pot[_qp])
                                + (MultiCompMultiPhaseBase::sum_dh_L_BC_epsilon()* _grad_C_diff_pot[_qp]) 
                                + (MultiCompMultiPhaseBase::sum_dh_L_BD_epsilon()* _grad_D_diff_pot[_qp]) ) *_phi[_j][_qp]);
 } 
 else
    return 0.0; 
}
