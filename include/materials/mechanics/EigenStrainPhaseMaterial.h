//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

class EigenStrainPhaseMaterial;

#include "Material.h"
#include "RankTwoTensor.h"

//template <>
//InputParameters validParams<EigenStrainPhaseMaterial>();

class EigenStrainPhaseMaterial : public Material
{
public:
  EigenStrainPhaseMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

  const std::string _base_name;

  //Input eigen strain
  const MaterialProperty<RankTwoTensor> & _eigen_tensor;

  //The first and second deriavtives of prefactor
  const MaterialProperty<Real> & _dprefactor;
  const MaterialProperty<Real> & _d2prefactor;

  //Declare the following properties
  MaterialProperty<RankTwoTensor> & _deigen_dmu;
  MaterialProperty<RankTwoTensor> & _d2eigen_dmu2;

};
