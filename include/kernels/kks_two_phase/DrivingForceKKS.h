//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef DRIVINGFORCEKKS_H
#define DRIVINGFORCEKKS_H

#include "Kernel.h"

class DrivingForceKKS;

template<>
InputParameters validParams<DrivingForceKKS>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *df/dphi = dh/dphi * (omega_beta - omega_alpha)
  *omega_beta - Grandpotential of the beta phase
  *omega_alpha - GrandPotential of the alpha phase
  **/

class DrivingForceKKS :  public Kernel
{
  public: 
    DrivingForceKKS(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    //Coupled variable diffusion potential of comp B
    const VariableValue & _B_diff_pot;
    unsigned int _B_diff_pot_var;
        
   //chemical potential of the dependent component in both phases
    const MaterialProperty<Real> & _A_chem_pot_alpha;
    const MaterialProperty<Real> & _A_chem_pot_beta;
    
   //mu(A).(mu(B) - mu(A)) = x_B
    const MaterialProperty<Real> & _xB_alpha;
    const MaterialProperty<Real> & _xB_beta; 
    
    //In addition this kernel requires the first and second derivatives
    // of the interpolation function
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _d2h;
    
    //The mobility here is assumed to be a constant
    const MaterialProperty<Real> & _L;
    
    //A non-dimensional parameter
    const MaterialProperty<Real> & _nd_factor;
    
};
#endif // DRIVINGFORCEKKS_H
