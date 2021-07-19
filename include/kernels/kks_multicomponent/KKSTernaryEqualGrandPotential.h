//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef KKSTERNARYEQUALGRANDPOTENTIAL_H
#define KKSTERNARYEQUALGRANDPOTENTIAL_H

#include "Kernel.h"
#include "KKSEqualGrandPotential.h"

class KKSTernaryEqualGrandPotential;

template<>
InputParameters validParams<KKSTernaryEqualGrandPotential>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *df/dphi = dh/dphi * (omega_beta - omega_alpha)
  *omega_beta - Grandpotential of the beta phase
  *omega_alpha - GrandPotential of the alpha phase
  **/

//This class is derived from the base class KKSEqualGrandPotential
class KKSTernaryEqualGrandPotential :  public KKSEqualGrandPotential
{
  public: 
    KKSTernaryEqualGrandPotential(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:
  
    //Declaring member function to calculate the diff in GP
    Real TernaryGPDiff();
      
    const VariableValue & _xC_alpha;
    unsigned int _xC_alpha_var;
    
    const VariableValue & _xC_beta;
    unsigned int _xC_beta_var;

    const MaterialProperty<Real> & _C_diff_pot_alpha;
    const MaterialProperty<Real> & _C_diff_pot_beta;
    const MaterialProperty<Real> & _C_therm_factor_alpha;
    const MaterialProperty<Real> & _C_therm_factor_beta;
    const MaterialProperty<Real> & _BC_therm_factor_alpha;
    const MaterialProperty<Real> & _BC_therm_factor_beta;    
};
#endif // KKSTERNARYEQUALGRANDPOTENTIAL_H
