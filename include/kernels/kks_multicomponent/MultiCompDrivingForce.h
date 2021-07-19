//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef MULTICOMPDRIVINGFORCE_H
#define MULTICOMPDRIVINGFORCE_H

#include "BinaryMultiPhaseDrivingForce.h"

class MultiCompDrivingForce;

template<>
InputParameters validParams<MultiCompDrivingForce>();

/**
  *This class enforces the 
  *Equality of chemical potentials of the dependent component
  *The kernel is implemented keeping in mind a quaternary alloy
  * A-B-C-D
  **/

class MultiCompDrivingForce :  public BinaryMultiPhaseDrivingForce
{
  public: 
    MultiCompDrivingForce(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:
      
    const VariableValue & _C_diff_pot;
    unsigned int _C_diff_pot_var;
    
    const VariableValue & _D_diff_pot;
    unsigned int _D_diff_pot_var;
    
    //mu(A).(mu(C) - mu(A)) = x(C)
    const MaterialProperty<Real> & _xC_1;
    const MaterialProperty<Real> & _xC_2;
    
    //mu(A).(mu(D) - mu(A)) = x(D)
    const MaterialProperty<Real> & _xD_1;
    const MaterialProperty<Real> & _xD_2;
       
};
#endif // MULTICOMPDRIVINGFORCE_H
