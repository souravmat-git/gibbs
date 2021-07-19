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

#include "TernaryContinuityEquationC.h"

registerMooseObject("gibbsApp", TernaryContinuityEquationC);

template<>
InputParameters
validParams<TernaryContinuityEquationC>()
{
  InputParameters params = validParams<TernaryMultiPhaseBase>();
  params.addClassDescription("Continuity equation for comp C");
  params.addRequiredCoupledVar("xB", "Mole fraction of comp B");
  params.addRequiredCoupledVar("C_diff_pot", "Diffusion potential of comp C");
  params.addRequiredCoupledVar("B_diff_pot", "Diffusion potential of comp B");
  return params;
}

TernaryContinuityEquationC::TernaryContinuityEquationC(const InputParameters & parameters)
  : TernaryMultiPhaseBase(parameters),
   //Mole fraction of componnet B
   _grad_xB(coupledGradient("xB")),
   _xB_var(coupled("xB")),
   //Diffusion potential of componnet C
   _grad_C_diff_pot(coupledGradient("C_diff_pot")),
   _C_diff_pot_var(coupled("C_diff_pot")),
   //Diffusion potential of componnet B
   _grad_B_diff_pot(coupledGradient("B_diff_pot")),
   _B_diff_pot_var(coupled("B_diff_pot"))
{
}

Real
TernaryContinuityEquationC::computeQpResidual()
{ 
   return (_grad_test[_i][_qp] *(TernaryMultiPhaseBase::L_CC_interp() * _grad_C_diff_pot[_qp] 
                               + TernaryMultiPhaseBase::L_BC_interp() * _grad_B_diff_pot[_qp]));
}

Real
TernaryContinuityEquationC::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * (TernaryMultiPhaseBase::DC_CC_interp() * _grad_phi[_j][_qp]) );
                                //+(TernaryMultiPhaseBase::dD_CC_xC_interp()* _grad_u[_qp]* _phi[_j][_qp])) ); 
}

Real
TernaryContinuityEquationC::computeQpOffDiagJacobian(unsigned int jvar)
{
 if (jvar == _xB_var)
 {
    return (_grad_test[_i][_qp] * TernaryMultiPhaseBase::DC_CB_interp() * _grad_phi[_j][_qp]);
 }
 else if (jvar == _C_diff_pot_var)
 { 
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::L_CC_interp()* _grad_phi[_j][_qp])
                                 +((TernaryMultiPhaseBase::dL_CC_muC_interp() * _grad_C_diff_pot[_qp])
                                 + (TernaryMultiPhaseBase::dL_BC_muC_interp() * _grad_B_diff_pot[_qp])
                                  ) * _phi[_j][_qp] ));
 }
 else if (jvar == _B_diff_pot_var)
 {
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::L_BC_interp()* _grad_phi[_j][_qp])
                                 +( (TernaryMultiPhaseBase::dL_CC_muB_interp() * _grad_C_diff_pot[_qp])
                                 +  (TernaryMultiPhaseBase::dL_BC_muB_interp() * _grad_B_diff_pot[_qp]) ) * _phi[_j][_qp] ));
 }

 else if (jvar == _phase_alpha_var)
 {
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::sum_dh_L_CC_alpha()* _grad_C_diff_pot[_qp])
                                + (TernaryMultiPhaseBase::sum_dh_L_BC_alpha()* _grad_B_diff_pot[_qp]) 
                                  ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_beta_var)
 {
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::sum_dh_L_CC_beta()* _grad_C_diff_pot[_qp])
                                + (TernaryMultiPhaseBase::sum_dh_L_BC_beta()*  _grad_B_diff_pot[_qp]) 
                                  ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_gamma_var)
 {
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::sum_dh_L_CC_gamma()* _grad_C_diff_pot[_qp])
                                + (TernaryMultiPhaseBase::sum_dh_L_BC_gamma()*  _grad_B_diff_pot[_qp]) 
                                  ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_delta_var)
 {
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::sum_dh_L_CC_delta()* _grad_C_diff_pot[_qp])
                                + (TernaryMultiPhaseBase::sum_dh_L_BC_delta()*  _grad_B_diff_pot[_qp]) 
                                 ) *_phi[_j][_qp]);
                                
 }
 else if (jvar == _phase_epsilon_var)
 {
    return (_grad_test[_i][_qp] *( (TernaryMultiPhaseBase::sum_dh_L_CC_epsilon()* _grad_C_diff_pot[_qp])
                                + (TernaryMultiPhaseBase::sum_dh_L_BC_epsilon()*  _grad_B_diff_pot[_qp]) 
                                 ) *_phi[_j][_qp]);
                                
 }  
 else //anything else
   return 0.0;
}
