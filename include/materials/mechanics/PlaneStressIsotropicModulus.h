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
class PlaneStressIsotropicModulus;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<PlaneStressIsotropicModulus>();

class PlaneStressIsotropicModulus : public Material
{
public:
  PlaneStressIsotropicModulus(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //std::string _C11_name, _C12_name, _C22_name, _C66_name;

    std::string _phase_name;

    const MaterialProperty<Real> & _E;
    const MaterialProperty<Real> & _nu;

    //Return the following quantities
    MaterialProperty<Real> & _C11_val;
    MaterialProperty<Real> & _C12_val;
    MaterialProperty<Real> & _C22_val;
    MaterialProperty<Real> & _C66_val;
};
