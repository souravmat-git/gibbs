//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef CHEMICALPOTENTIALDEPCOMPONENTBASE_H
#define CHEMICALPOTENTIALDEPCOMPONENTBASE_H

#include "Kernel.h"

class ChemicalPotentialDepComponentBase;

template<>
InputParameters validParams<ChemicalPotentialDepComponentBase>();

/**
  *This is a base kernel for binary alloys
  * for implementing the chemical potential diff 
  * of the dependent component
  **/

class ChemicalPotentialDepComponentBase : public Kernel
{
  public: 
    ChemicalPotentialDepComponentBase(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
 
    //private: Member function defined within private cannot be used 
    // for derived class

    const VariableValue & _xB_alpha;
    unsigned int _xB_alpha_var;
    
    const VariableValue & _xB_beta;
    unsigned int _xB_beta_var;
    
    //chemical potential of the dependent component in alpha phase
    const MaterialProperty<Real> & _A_chem_pot_alpha;
    const MaterialProperty<Real> & _A_chem_pot_beta;
    
    //Thermodynamic factor of comp A w.r.t B in two phases
    const MaterialProperty<Real> & _AB_therm_factor_alpha;
    const MaterialProperty<Real> & _AB_therm_factor_beta; 
    
    //In addition this kernel requires the first and second derivatives
    // of the interpolation function
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _d2h;
    
    //The mobility here is assumed to be a constant
    const MaterialProperty<Real> & _L;
    
};
#endif // CHEMICALPOTENTIALDEPCOMPONENTBASE_H
