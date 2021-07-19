//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef DRIVINGFORCESPINODAL_H
#define DRIVINGFORCESPINODAL_H

#include "Kernel.h"

class DrivingForceSpinodal;

template<>
InputParameters validParams<DrivingForceSpinodal>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *df/dphi = dh/dphi * (omega_beta - omega_alpha)
  *omega_beta - Grandpotential of the beta phase
  *omega_alpha - GrandPotential of the alpha phase
  **/

class DrivingForceSpinodal :  public Kernel
{
  public: 
    DrivingForceSpinodal(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

   //Coupled variable is the composition of B
   const VariableValue & _xB;
   unsigned int _xB_var;
        
   //Free energy of the spinodal
   const MaterialProperty<Real> & _f_alpha;
    
   //Diffusion potential of the spinodal
   const MaterialProperty<Real> & _B_diff_pot_alpha;
     
   //mu(A).(mu(B) - mu(A)) = x_B
   const MaterialProperty<Real> & _B_therm_factor_alpha;
    
   //eqm_chemical potential
   const MaterialProperty<Real> & _A_chem_pot_eqm; 
    
   //In addition this kernel requires the first and second derivatives
   // of the interpolation function
   const MaterialProperty<Real> & _dh;
   const MaterialProperty<Real> & _d2h;
    
   //The mobility here is assumed to be a constant
   const MaterialProperty<Real> & _L;    
};
#endif // DRIVINGFORCESPINODAL_H
