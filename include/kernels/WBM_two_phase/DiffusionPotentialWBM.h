//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef DIFFUSIONPOTENTIALWBM_H
#define DIFFUSIONPOTENTIALWBM_H

#include "Kernel.h"

class DiffusionPotentialWBM;

template<>
InputParameters validParams<DiffusionPotentialWBM>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *mu - df/dc= 0
  *mu - is the global diffusion potential
  *df/dc = h*df_beta/dc + [1-h]*df_alpha/dc
  * The kernel operates on the variable : c
  **/

class DiffusionPotentialWBM :  public Kernel
{
  public: 
    DiffusionPotentialWBM(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
    
    const VariableValue & _eta;
    unsigned int _eta_var;
    
    const VariableValue & _B_diff_pot;
    unsigned int _B_diff_pot_var;
    
    Real _interpolated_mu;
    
    //Keep in mind the free energy are functions of 
    //composition and not phase compositions
    //Note the difference in the naming of the two variables 
    //dfalpha_dc and df_dcalpha
    const MaterialProperty<Real> & _h;
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _B_diff_pot_alpha;
    const MaterialProperty<Real> & _B_therm_factor_alpha;
    const MaterialProperty<Real> & _B_diff_pot_beta;
    const MaterialProperty<Real> & _B_therm_factor_beta;

};
#endif // DIFFUSIONPOTENTIALWBM_H
