//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class ProjectionTensor;

//MOOSE includes
#include "Material.h"
#include "RankTwoTensor.h"

template <>
InputParameters validParams<ProjectionTensor>();

class ProjectionTensor : public Material
{
public:
  ProjectionTensor(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    //Independent variable on which this property depends
    const VariableGradient & _grad_eta;
    
    //Declare the second rank projection tensor to calculate
    MaterialProperty<RankTwoTensor> &  _projection_tensor;    
};
