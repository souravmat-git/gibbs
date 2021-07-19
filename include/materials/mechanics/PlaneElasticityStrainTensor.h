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
class PlaneElasticityStrainTensor;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<PlaneElasticityStrainTensor>();

class PlaneElasticityStrainTensor : public Material
{
public:
  PlaneElasticityStrainTensor(const InputParameters & parameters);
  
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
        
    MaterialProperty<Real> & _exx_val;   
    MaterialProperty<Real> & _eyy_val;
    MaterialProperty<Real> & _exy_val;  
};
