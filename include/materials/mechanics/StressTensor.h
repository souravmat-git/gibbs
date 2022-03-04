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
class StressTensor;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<StressTensor>();

class StressTensor : public Material
{
public:
  StressTensor(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Independent variable on which this property depends
    const VariableGradient & _grad_ux;

    //Youngs modulus
    const MaterialProperty<Real> & _E;

    //Stress in the x-direction
    MaterialProperty<Real> & _sx_val;
};
