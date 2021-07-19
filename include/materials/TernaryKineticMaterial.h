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

//#ifndef TERNARYKINETICMATERIAL_H
//#define TERNARYKINETICMATERIAL_H

#pragma once

// Forward Declarations
class TernaryKineticMaterial;

//MOOSE includes
#include "TabulatedKineticMaterial.h"
#include "TernaryMobilityData.h"

template <>
InputParameters validParams<TernaryKineticMaterial>();

class TernaryKineticMaterial : public TabulatedKineticMaterial
{
public:
  TernaryKineticMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    
    //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string _L_BB_name, _L_BC_name, _L_CC_name,
                _dL_BB_xB_name, _dL_BC_xB_name, _dL_CC_xB_name,
                _dL_BB_xC_name, _dL_BC_xC_name, _dL_CC_xC_name;
    
    //Onsager mobility matrix
   MaterialProperty<Real> & _L_BB_val;
   MaterialProperty<Real> & _L_BC_val;
   MaterialProperty<Real> & _L_CC_val;
    
   //LBB with respect to XB, xC
   MaterialProperty<Real> & _dL_BB_xB_val;
   MaterialProperty<Real> & _dL_BB_xC_val;
        
   //LBC with respect to XB, xC
   MaterialProperty<Real> & _dL_BC_xB_val;
   MaterialProperty<Real> & _dL_BC_xC_val;
    
   //LBC with respect to XB, xC
   MaterialProperty<Real> & _dL_CC_xB_val; 
   MaterialProperty<Real> & _dL_CC_xC_val;    
    
    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_B;
    
    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_C;
    
    //TernaryPhaseData is an userObject, which reads the 
    //tables for ternary alloys with the properties of the phase
    //and linearly interpolates the value .
    const TernaryMobilityData & _table_object;
 
};
//#endif // TERNARYKINETICMATERIAL_H
