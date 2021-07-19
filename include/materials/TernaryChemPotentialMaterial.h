//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This code is based on the code 
//*from TabulatedFluidProperties.C

//#ifndef TERNARYCHEMPOTENTIALMATERIAL_H
//#define TERNARYCHEMPOTENTIALMATERIAL_H

#pragma once

// Forward Declarations
class TernaryChemPotentialMaterial;

//MOOSE includes
#include "TabulatedPhaseMaterial.h"
#include "TernaryChemPotentialData.h"

template <>
InputParameters validParams<TernaryChemPotentialMaterial>();

class TernaryChemPotentialMaterial : public TabulatedPhaseMaterial
{
public:
  TernaryChemPotentialMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    
    //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string _A_chem_pot_name, _AB_therm_factor_name, _AC_therm_factor_name;
    
    //chemical potential of comp A that this material interpolates and reurns
   MaterialProperty<Real> & _A_chem_pot_val;
      
   //Thermodynamic factor mu(A).x(B) that this material interpolates and returns
    MaterialProperty<Real> & _AB_therm_factor_val;
    
    //Thermodynamic factor w.r.t AC that this material interpolates and returns
    MaterialProperty<Real> & _AC_therm_factor_val;
    
    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_B;
    
    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_C;
    
    //TernaryChemPotentialData is an userObject, which reads the 
    //tables for ternary alloys with the properties of the phase
    //and linearly interpolates the value .
    const TernaryChemPotentialData & _table_object;
 
};
//#endif // TERNARYCHEMPOTENTIALMATERIAL_H
