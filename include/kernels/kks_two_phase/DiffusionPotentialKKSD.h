//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef DIFFUSIONPOTENTIALKKSD_H
#define DIFFUSIONPOTENTIALKKSD_H

#include "Kernel.h"

class DiffusionPotentialKKSD;

template<>
InputParameters validParams<DiffusionPotentialKKSD>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *mu - df/dc_beta = 0
  *mu - Diffusion potential df/dc = df/dc_alpha = df/dc_beta
  *df/dc_beta - first derivate of the free energy of the beta
  * phase. This is assumed to be a parabolic form
  * The kernel operates on the variable : c
  **/

class DiffusionPotentialKKSD :  public Kernel
{
  public: 
    DiffusionPotentialKKSD(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:  
    const VariableValue & _xD_beta;
    unsigned int _xD_beta_var;
    
    //This is the coupled required for 
    //ternary and above alloys
    const VariableValue & _xB_beta;
    unsigned int _xB_beta_var;
    
    const VariableValue & _xC_beta;
    unsigned int _xC_beta_var;
    
    const VariableValue & _xD_alpha;
    unsigned int _xD_alpha_var;
    
    const VariableValue & _xD_gamma;
    unsigned int _xD_gamma_var;

    const VariableValue & _diff_pot;
    unsigned int _diff_pot_var;
    
    //Diffusion potential of component C
    const MaterialProperty<Real> & _D_diff_pot_beta;
    
    //Therm factor w.r.t -> C mu(C).x(C)
    const MaterialProperty<Real> & _D_therm_factor_beta;
    
    //Therm factor w.r.t -> mu(C).x(B)
    const MaterialProperty<Real> & _BD_therm_factor_beta;
    
    //Therm factor w.r.t-> mu(C).x(D)
    const MaterialProperty<Real> & _CD_therm_factor_beta;
};
#endif // _DIFFUSIONPOTENTIALKKSD_H
