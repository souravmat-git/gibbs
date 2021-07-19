//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef MCCONTINUITYEQUATIOND_H
#define MCCONTINUITYEQUATIOND_H

#include "MultiCompMultiPhaseBase.h"

class MCContinuityEquationD;

template<>
InputParameters validParams<MCContinuityEquationD>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *w* mu  = 0
  *mu - Diffusion potential df/dc = df/dc_alpha = df/dc_beta
  *df/dc_beta - first derivate of the free energy of the beta
  * phase. This is assumed to be a parabolic form
  * The kernel operates on the variable : c
  **/

class MCContinuityEquationD:  public MultiCompMultiPhaseBase
{
  public: 
     MCContinuityEquationD(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:
  
    //Mole fraction of component C
    const VariableGradient & _grad_xC;
    unsigned int _xC_var;
    
    //Mole fraction of component B
    const VariableGradient & _grad_xB;
    unsigned int _xB_var;
    
     //Diffusion potential of comp. D
    const VariableGradient & _grad_D_diff_pot;
    unsigned int _D_diff_pot_var; 
    
    //Diffusion potential of comp. B
    const VariableGradient & _grad_B_diff_pot;
    unsigned int _B_diff_pot_var;
    
    //Diffusion potential of comp. C
    const VariableGradient & _grad_C_diff_pot;
    unsigned int _C_diff_pot_var;
    
};
#endif //MCCONTINUITYEQUATIOND_H
