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
class ElasticFieldMaterial;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<ElasticFieldMaterial>();

class ElasticFieldMaterial : public Material
{
public:
  ElasticFieldMaterial(const InputParameters & parameters);
  
  //This class calculates the strain, stress and the elastic energy
  //of a phase with or without a prescribed eigenstrain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    //Independent variable on which this property depends
    const VariableGradient & _grad_ux;
    
    //Input material properties
    //for this material class are youngs modulus, poissons ratio and transformation strain
    const MaterialProperty<Real> & _eT;
    const MaterialProperty<Real> & _mat_const;
    
    //Strain in the x-direction   
    MaterialProperty<Real> & _ex_val; 
    //Stress in the x-direction
    MaterialProperty<Real> & _sx_val;
    //Elastic strain energy
    MaterialProperty<Real> & _elastic_energy_val;    
};
