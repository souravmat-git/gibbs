//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef KKSEQUALGRANDPOTENTIAL_H
#define KKSEQUALGRANDPOTENTIAL_H

#include "Kernel.h"

class KKSEqualGrandPotential;

template<>
InputParameters validParams<KKSEqualGrandPotential>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *df/dphi = dh/dphi * (omega_beta - omega_alpha)
  *omega_beta - Grandpotential of the beta phase
  *omega_alpha - GrandPotential of the alpha phase
  **/

class KKSEqualGrandPotential :  public Kernel
{
  public: 
    KKSEqualGrandPotential(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

//private: Member function defined within private cannot be used 
// for derived class

    Real grandPotentials_diff();

    const VariableValue & _xB_alpha;
    unsigned int _xB_alpha_var;
    
    const VariableValue & _xB_beta;
    unsigned int _xB_beta_var;
        
    const MaterialProperty<Real> & _f_alpha;
    const MaterialProperty<Real> & _f_beta;
    const MaterialProperty<Real> & _B_diff_pot_alpha;
    const MaterialProperty<Real> & _B_diff_pot_beta;
    const MaterialProperty<Real> & _B_therm_factor_alpha;
    const MaterialProperty<Real> & _B_therm_factor_beta;

    
    //In addition this kernel requires the first and second derivatives
    // of the interpolation function
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _d2h;
    
    //The mobility here is assumed to be a constant
    const MaterialProperty<Real> & _L;
    
};
#endif // KKSEQUALGRANDPOTENTIAL_H
