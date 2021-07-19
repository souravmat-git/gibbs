//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef CHEMICALPOTENTIALDEPCOMPONENT_H
#define CHEMICALPOTENTIALDEPCOMPONENT_H

#include "Kernel.h"
#include "ChemicalPotentialDepComponentBase.h"

class ChemicalPotentialDepComponent;

template<>
InputParameters validParams<ChemicalPotentialDepComponent>();

/**
  *This class enforces the 
  *Equality of chemical potentials of the dependent component
  *The kernel is implemented keeping in mind a quaternary alloy
  * A-B-C-D
  **/


class ChemicalPotentialDepComponent :  public ChemicalPotentialDepComponentBase
{
  public: 
    ChemicalPotentialDepComponent(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:

    const VariableValue & _xC_alpha;
    unsigned int _xC_alpha_var;
    
    const VariableValue & _xC_beta;
    unsigned int _xC_beta_var;
    
    const VariableValue & _xD_alpha;
    unsigned int _xD_alpha_var;
    
    const VariableValue & _xD_beta;
    unsigned int _xD_beta_var;
    
    //Thermodynamic factor of comp A w.r.t C in two phases
    const MaterialProperty<Real> & _AC_therm_factor_alpha;
    const MaterialProperty<Real> & _AC_therm_factor_beta;
    
    //Thermodynamic factor of comp A w.r.t D in two phases
    const MaterialProperty<Real> & _AD_therm_factor_alpha;
    const MaterialProperty<Real> & _AD_therm_factor_beta;
       
};
#endif // CHEMICALPOTENTIALDEPCOMPONENT_H
