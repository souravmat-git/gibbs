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
class PlaneStrainIsotropicModulus;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<PlaneStrainIsotropicModulus>();

class PlaneStrainIsotropicModulus : public Material
{
public:
  PlaneStrainIsotropicModulus(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    std::string _phase_name;

    const MaterialProperty<Real> & _E;
    const MaterialProperty<Real> & _nu;

    //Return the following quantities
    MaterialProperty<Real> & _C11_val;
    MaterialProperty<Real> & _C12_val;
    MaterialProperty<Real> & _C22_val;
    MaterialProperty<Real> & _C66_val;
};
