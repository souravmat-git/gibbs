//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Forward Declarations
class QuaternaryKineticMaterial;

//MOOSE includes
#include "TabulatedKineticMaterial.h"
#include "QuaternaryMobilityData.h"

template <>
InputParameters validParams<QuaternaryKineticMaterial>();

class QuaternaryKineticMaterial : public TabulatedKineticMaterial
{
public:
  QuaternaryKineticMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    
    //String variable to hold the ceff, of the mobility 
    // matrix obtained from the input file
    std::string _L_BB_name, _L_CC_name, _L_DD_name,
                _L_BC_name, _L_BD_name, _L_CD_name,
                _dL_BB_xB_name, _dL_CC_xB_name, _dL_DD_xB_name,
                _dL_BC_xB_name, _dL_BD_xB_name, _dL_CD_xB_name,
                _dL_BB_xC_name, _dL_CC_xC_name, _dL_DD_xC_name,
                _dL_BC_xC_name, _dL_BD_xC_name, _dL_CD_xC_name,
                _dL_BB_xD_name, _dL_CC_xD_name, _dL_DD_xD_name,
                _dL_BC_xD_name, _dL_BD_xD_name, _dL_CD_xD_name;
    
   //Onsager mobility coeff diagonal terms
   MaterialProperty<Real> & _L_BB_val;
   MaterialProperty<Real> & _L_CC_val;
   MaterialProperty<Real> & _L_DD_val;
   
   //Onsager mobility coeff off-diagonal terms
   MaterialProperty<Real> & _L_BC_val;
   MaterialProperty<Real> & _L_BD_val;
   MaterialProperty<Real> & _L_CD_val;
    
   //Diagonal components with respect to xB
   MaterialProperty<Real> & _dL_BB_xB_val;        
   MaterialProperty<Real> & _dL_CC_xB_val;
   MaterialProperty<Real> & _dL_DD_xB_val;
   
   MaterialProperty<Real> & _dL_BC_xB_val;        
   MaterialProperty<Real> & _dL_BD_xB_val;
   MaterialProperty<Real> & _dL_CD_xB_val;
    
   //Diagonal components with respect to xC
   MaterialProperty<Real> & _dL_BB_xC_val;        
   MaterialProperty<Real> & _dL_CC_xC_val;
   MaterialProperty<Real> & _dL_DD_xC_val;
   
   MaterialProperty<Real> & _dL_BC_xC_val;        
   MaterialProperty<Real> & _dL_BD_xC_val;
   MaterialProperty<Real> & _dL_CD_xC_val; 
   
   //Diagonal components with respect to xD
   MaterialProperty<Real> & _dL_BB_xD_val;        
   MaterialProperty<Real> & _dL_CC_xD_val;
   MaterialProperty<Real> & _dL_DD_xD_val;
   
   MaterialProperty<Real> & _dL_BC_xD_val;        
   MaterialProperty<Real> & _dL_BD_xD_val;
   MaterialProperty<Real> & _dL_CD_xD_val;   
    
   //Independent variable on which this property depends
   const VariableValue & _xB;
    
   //Independent variable on which this property depends
   const VariableValue & _xC;
    
   //Independent variable on which this property depends
   const VariableValue & _xD;
    
   //QuaternaryMobilityData is an userObject, which reads the 
   //tables for quaternary alloys with the properties of the phase
   //and trilinearly interpolates the value .
   const QuaternaryMobilityData & _table_object;
 
};
