//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class BinaryConstantKineticMaterial;

//MOOSE includes
#include "Material.h"

//template <>
//InputParameters validParams<BinaryConstantKineticMaterial>();

class BinaryConstantKineticMaterial : public Material
{
public:
  BinaryConstantKineticMaterial(const InputParameters & parameters);

  static InputParameters validParams();

private:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

  //String variable to hold L_BB and dL_BB_muB
  std::string  _L_BB_name, _dL_BB_muB_name;

  //L_BB value that this material interpolates and returns
  MaterialProperty<Real> & _L_BB;
  //dL_BB_muB value that this material interpolates and reurns
  MaterialProperty<Real> & _dL_BB_muB;

  const Real & _L_BB_val;
};
