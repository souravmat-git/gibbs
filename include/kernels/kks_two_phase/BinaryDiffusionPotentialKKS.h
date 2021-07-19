//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#ifndef BINARYDIFFUSIONPOTENTIALKKS_H
#define BINARYDIFFUSIONPOTENTIALKKS_H

#include "Kernel.h"

class BinaryDiffusionPotentialKKS;

template<>
InputParameters validParams<BinaryDiffusionPotentialKKS>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *mu - df/dc_beta = 0
  *mu - Diffusion potential df/dc = df/dc_alpha = df/dc_beta
  *df/dc_beta - first derivate of the free energy of the beta
  **/

class BinaryDiffusionPotentialKKS :  public Kernel
{
  public: 
    BinaryDiffusionPotentialKKS(const InputParameters & parameters);
  
  protected:
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _xB_beta;
  unsigned int _xB_beta_var;
    
  const VariableValue & _diff_pot;
  unsigned int _diff_pot_var;
    
  //Diffusion potential of component B
  const MaterialProperty<Real> & _B_diff_pot_beta;
    
   //Thermodynamic factor w.r.t  component B
  const MaterialProperty<Real> & _B_therm_factor_beta; 
};
#endif // BINARYDIFFUSIONPOTENTIALKKS_H
