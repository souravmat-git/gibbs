//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

class TernaryConstantKineticMaterial;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<TernaryConstantKineticMaterial>();

class TernaryConstantKineticMaterial : public Material
{
public:
  TernaryConstantKineticMaterial(const InputParameters & parameters);

private:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
 
  //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string _L_BB_name, _L_BC_name, _L_CC_name,
                _dL_BB_muB_name, _dL_BC_muB_name, _dL_CC_muB_name,
                _dL_BB_muC_name, _dL_BC_muC_name, _dL_CC_muC_name;
    
   //Coefficients of the Onsager mobility matrix
   MaterialProperty<Real> & _L_BB;   
   MaterialProperty<Real> & _L_BC;
   MaterialProperty<Real> & _L_CC;
      
   //LBB with respect to muB, muC
   MaterialProperty<Real> & _dL_BB_muB;   
   MaterialProperty<Real> & _dL_BB_muC;
        
   //LBC with respect to muB, muC
   MaterialProperty<Real> & _dL_BC_muB;   
   MaterialProperty<Real> & _dL_BC_muC;
    
   //LBC with respect to muB, muC
   MaterialProperty<Real> & _dL_CC_muB;   
   MaterialProperty<Real> & _dL_CC_muC;  
   
   //Value of the coefficients of the Onsager mobility matrix
   const Real & _L_BB_val;
   const Real & _L_BC_val;
   const Real & _L_CC_val;
};
