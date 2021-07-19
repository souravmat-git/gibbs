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
class PlaneElasticityStressTensor;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<PlaneElasticityStressTensor>();

class PlaneElasticityStressTensor : public Material
{
public:
  PlaneElasticityStressTensor(const InputParameters & parameters);
  
  //This class calculates the strain from the displacement gradient
  //Using the strain-displacement relations.

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    //Independent variables on which this property depends
    // are the displacement gradient vector along x 
    // and the displacement gradient vector along y
    
    const VariableGradient & _grad_ux;
    const VariableGradient & _grad_uy;
    
    std::string _phase_name;
    
    //For orthorhombic materials and 2D 4 elastic constants are required 
    const MaterialProperty<Real> & _C11;
    const MaterialProperty<Real> & _C12;
    const MaterialProperty<Real> & _C22;
    const MaterialProperty<Real> & _C66;
        
    MaterialProperty<Real> & _sxx_val;   
    MaterialProperty<Real> & _syy_val;
    MaterialProperty<Real> & _sxy_val;  
};
