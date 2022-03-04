//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once
class VariableUnitNormal_MP;

//MOOSE includes
#include "Material.h"

//template <>
//InputParameters validParams<VariableUnitNormal_MP>();

class VariableUnitNormal_MP : public Material
{
public:
  VariableUnitNormal_MP(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

  const VariableGradient & _grad_phialpha;
  const VariableGradient & _grad_phibeta;
  const VariableGradient & _grad_phigamma;

  //Define unit normals
  MaterialProperty<RealVectorValue> & _n1;
  MaterialProperty<RealVectorValue> & _n2;

};
