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

//#ifndef TERNARYCONJUGATEKINETICMATERIAL_H
//#define TERNARYCONJUGATEKINETICMATERIAL_H

#pragma once

// Forward Declarations
class TernaryConjugateKineticMaterial;

//MOOSE includes
#include "TabulatedKineticMaterial.h"
#include "TernaryConjugateMobilityData.h"

template <>
InputParameters validParams<TernaryConjugateKineticMaterial>();

class TernaryConjugateKineticMaterial : public TabulatedKineticMaterial
{
public:
  TernaryConjugateKineticMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    
    //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string _L_BB_name, _L_BC_name, _L_CC_name,
                _dL_BB_muB_name, _dL_BC_muB_name, _dL_CC_muB_name,
                _dL_BB_muC_name, _dL_BC_muC_name, _dL_CC_muC_name;
    
    //Onsager mobility L_BB
   MaterialProperty<Real> & _L_BB_val;
    
    //Onsager mobility L_BC
   MaterialProperty<Real> & _L_BC_val;
   
   //Onsager mobility L_CC
   MaterialProperty<Real> & _L_CC_val;
    
   //LBB with respect to muB, muC
   MaterialProperty<Real> & _dL_BB_muB_val;    
   MaterialProperty<Real> & _dL_BB_muC_val;
        
   //LBC with respect to muB, muC
   MaterialProperty<Real> & _dL_BC_muB_val;   
   MaterialProperty<Real> & _dL_BC_muC_val;
    
   //LBC with respect to muB, muC
   MaterialProperty<Real> & _dL_CC_muB_val;   
   MaterialProperty<Real> & _dL_CC_muC_val;    
    
    //Independent variable on which this property depends
   const VariableValue & _B_diff_pot;
    
    //Independent variable on which this property depends
   const VariableValue & _C_diff_pot;
    
    //TernaryPhaseData is an userObject, which reads the 
    //tables for ternary alloys with the properties of the phase
    //and linearly interpolates the value .
   const TernaryConjugateMobilityData & _table_object;
};
