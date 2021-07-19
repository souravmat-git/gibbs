//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

class VariableUnitNormalBase;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<VariableUnitNormalBase>();

class VariableUnitNormalBase : public Material
{
public:
  VariableUnitNormalBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  //Declare three member functions to check the
  //gradient of alpha, beta and gamma phases
  Real  unalpha_norm_sq() const;
  Real  unbeta_norm_sq() const;
  Real  ungamma_norm_sq() const;

  const VariableGradient & _grad_phialpha;
  const VariableGradient & _grad_phibeta;
  const VariableGradient & _grad_phigamma;

};
