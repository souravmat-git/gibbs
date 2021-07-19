//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef TERNARYCONTINUITYEQUATIONC_H
#define TERNARYCONTINUITYEQUATIONC_H

#include "TernaryMultiPhaseBase.h"

class TernaryContinuityEquationC;

template<>
InputParameters validParams< TernaryContinuityEquationC>();

/**
  *This class enforces the continuity equation on component C
  * The kernel operates on the variable : Xc
  **/

class  TernaryContinuityEquationC:  public TernaryMultiPhaseBase
{
  public: 
     TernaryContinuityEquationC(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:
    
    //Mole fraction of component B
    const VariableGradient & _grad_xB;
    unsigned int _xB_var;
    
    //Diffusion potential of comp. C
    const VariableGradient & _grad_C_diff_pot;
    unsigned int _C_diff_pot_var;
  
    //Diffusion potential of comp. B
    const VariableGradient & _grad_B_diff_pot;
    unsigned int _B_diff_pot_var;
    
    
};
#endif //TERNARYCONTINUITYEQUATIONC_H