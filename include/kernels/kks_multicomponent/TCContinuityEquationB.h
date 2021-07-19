//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#pragma once
#include "TernaryMultiPhaseBase.h"

class TCContinuityEquationB;

template<>
InputParameters validParams<TCContinuityEquationB>();

/**
  *This class enforces the continuity equation for mass
  * The kernel operates on the variable : xB
  * for a ternary A-B-C alloy
  **/

class  TCContinuityEquationB :  public TernaryMultiPhaseBase
{
  public: 
     TCContinuityEquationB(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:  
  
    //Note all phase-field variables 
    //are inherited from the base class
    
    //Only chemical and mechanical variables 
    //need to be coupled
     
    //Mole fraction of C
    const VariableGradient & _grad_xC;
    unsigned int _xC_var;
      
    //Diffusion potential of comp. B
    const VariableGradient & _grad_B_diff_pot;
    unsigned int _B_diff_pot_var;
    
    //Diffusion potential of comp. C
    const VariableGradient & _grad_C_diff_pot;
    unsigned int _C_diff_pot_var;
       
};
