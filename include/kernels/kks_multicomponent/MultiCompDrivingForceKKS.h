//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef MULTICOMPDRIVINGFORCEKKS_H
#define MULTICOMPDRIVINGFORCEKKS_H

#include "Kernel.h"
#include "DrivingForceKKS.h"

class MultiCompDrivingForceKKS;

template<>
InputParameters validParams<MultiCompDrivingForceKKS>();

/**
  *This class enforces the 
  *Equality of chemical potentials of the dependent component
  *The kernel is implemented keeping in mind a quaternary alloy
  * A-B-C-D
  **/


class MultiCompDrivingForceKKS :  public DrivingForceKKS
{
  public: 
    MultiCompDrivingForceKKS (const InputParameters & parameters);
  
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
    
    //mu(A).(mu(C) - mu(A)) = x(A)
    const MaterialProperty<Real> & _xC_alpha;
    const MaterialProperty<Real> & _xC_beta;
    
    //mu(A).(mu(D) - mu(A))
    const MaterialProperty<Real> & _xD_alpha;
    const MaterialProperty<Real> & _xD_beta;
       
};
#endif // MULTICOMPDRIVINGFORCEKKS_H
