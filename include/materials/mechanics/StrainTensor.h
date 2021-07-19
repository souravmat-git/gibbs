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
class StrainTensor;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<StrainTensor>();

class StrainTensor : public Material
{
public:
  StrainTensor(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    //Independent variable on which this property depends
    const VariableGradient & _grad_ux;
    
    //Strain in the x-direction
    MaterialProperty<Real> & _ex_val;    
 
};
